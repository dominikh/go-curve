package curve

import (
	"iter"
	"math"
	"sort"
)

const maxSplineSplit = 100

var _ Shape = CubicBez{}
var _ ParametricCurve = CubicBez{}

type CubicBez struct {
	P0 Point
	P1 Point
	P2 Point
	P3 Point
}

func (q CubicBez) IsInf() bool {
	return q.P0.IsInf() || q.P1.IsInf() || q.P2.IsInf() || q.P3.IsInf()
}

func (q CubicBez) IsNaN() bool {
	return q.P0.IsNaN() || q.P1.IsNaN() || q.P2.IsNaN() || q.P3.IsNaN()
}

// BoundingBox implements [Shape].
func (c CubicBez) BoundingBox() Rect {
	return BoundingBox(c)
}

// PathElements implements [Shape].
func (c CubicBez) PathElements(tolerance float64) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		_ = yield(MoveTo(c.P0)) &&
			yield(CubicTo(c.P1, c.P2, c.P3))
	}
}

// Perimeter implements [Shape].
func (c CubicBez) Perimeter(accuracy float64) float64 {
	return c.Arclen(accuracy)
}

// Winding implements [Shape].
func (c CubicBez) Winding(pt Point) int {
	return 0
}

// Arclen returns the arclength of a cubic Bézier segment.
//
// This is an adaptive subdivision approach using Legendre-Gauss quadrature
func (c CubicBez) Arclen(accuracy float64) float64 {
	return c.arclen(accuracy, 0)
}

func (c CubicBez) arclen(accuracy float64, depth int) float64 {
	d03 := c.P3.Sub(c.P0)
	d01 := c.P1.Sub(c.P0)
	d12 := c.P2.Sub(c.P1)
	d23 := c.P3.Sub(c.P2)
	lplc := d01.Hypot() + d12.Hypot() + d23.Hypot() - d03.Hypot()
	dd1 := d12.Sub(d01)
	dd2 := d23.Sub(d12)
	// It might be faster to do direct multiplies, the data dependencies would be shorter.
	// The following values don't have the factor of 3 for first deriv
	dm := d01.Add(d23).Mul(0.25).Add(d12.Mul(0.5)) // first derivative at midpoint
	dm1 := dd2.Add(dd1).Mul(0.5)                   // second derivative at midpoint
	dm2 := dd2.Sub(dd1).Mul(0.25)                  // 0.5 * (third derivative at midpoint)

	var est float64
	for _, coeff := range gaussLegendreCoeffs8 {
		wi, xi := coeff[0], coeff[1]
		dNorm2 := dm.Add(dm1.Mul(xi)).Add(dm2.Mul(xi * xi)).Hypot2()
		ddNorm2 := dm1.Add(dm2.Mul(2.0 * xi)).Hypot2()
		f := ddNorm2 / dNorm2
		est += wi * f
	}
	if math.IsNaN(est) {
		// dNorm2 will be 0 as c approaches a singularity
		est = 0
	}

	estGauss8Error := min(math.Pow(est, 3)*2.5e-6, 3e-2) * lplc
	if estGauss8Error < accuracy {
		return arclenQuadratureCore(gaussLegendreCoeffs8Half[:], dm, dm1, dm2)
	}
	estGauss16Error := min(math.Pow(est, 6)*1.5e-11, 9e-3) * lplc
	if estGauss16Error < accuracy {
		return arclenQuadratureCore(gaussLegendreCoeffs16Half[:], dm, dm1, dm2)
	}
	estGauss24Error := min(math.Pow(est, 9)*3.5e-16, 3.5e-3) * lplc
	if estGauss24Error < accuracy || depth >= 20 {
		return arclenQuadratureCore(gaussLegendreCoeffs24Half[:], dm, dm1, dm2)
	}
	c0, c1 := c.Subdivide()
	return c0.arclen(accuracy*0.5, depth+1) + c1.arclen(accuracy*0.5, depth+1)
}

func (cb CubicBez) Eval(t float64) Point {
	mt := 1.0 - t
	a := Vec2(cb.P0).Mul(mt * mt * mt)
	b := Vec2(cb.P1).Mul(mt * mt * 3.0)
	c := Vec2(cb.P2).Mul(mt * 3.0)
	d := Vec2(cb.P3)
	v := a.Add(b.Add(c.Add(d.Mul(t)).Mul(t)).Mul(t))
	return Point(v)
}

// Subdivide subdivides the cubic into halves, using de Casteljau.
func (c CubicBez) Subdivide() (CubicBez, CubicBez) {
	pm := c.Eval(0.5)
	return CubicBez{
			c.P0,
			c.P0.Midpoint(c.P1),
			Point(Vec2(c.P0).Add(Vec2(c.P1).Mul(2.0)).Add(Vec2(c.P2)).Mul(0.25)),
			pm,
		},
		CubicBez{
			pm,
			Point(Vec2(c.P1).Add(Vec2(c.P2).Mul(2.0)).Add(Vec2(c.P3)).Mul(0.25)),
			c.P2.Midpoint(c.P3),
			c.P3,
		}
}

// SubdivideCurve subdivides the cubic into halves, using de Casteljau.
func (c CubicBez) SubdivideCurve() (ParametricCurve, ParametricCurve) {
	return c.Subdivide()
}

func (c CubicBez) Start() Point {
	return c.P0
}

func (c CubicBez) End() Point {
	return c.P3
}

type CubicToQuadraticSegment struct {
	Start, End float64
	Segment    QuadBez
}

// Quadratics converts the cubic Béziers to quadratic Béziers.
//
// The iterator returns the start and end parameter in the cubic of each quadratic
// segment, along with the quadratic.
//
// Note that the resulting quadratic Béziers are not in general G1 continuous;
// they are optimized for minimizing distance error.
//
// This iterator will always produce at least one value.
func (c CubicBez) Quadratics(accuracy float64) iter.Seq[CubicToQuadraticSegment] {
	// The maximum error, as a vector from the cubic to the best approximating
	// quadratic, is proportional to the third derivative, which is constant
	// across the segment. Thus, the error scales down as the third power of
	// the number of subdivisions. Our strategy then is to subdivide t evenly.
	//
	// This is an overestimate of the error because only the component
	// perpendicular to the first derivative is important. But the simplicity is
	// appealing.

	return func(yield func(CubicToQuadraticSegment) bool) {
		// This magic number is the square of 36 / sqrt(3).
		// See: https://web.archive.org/web/20210108052742/http://caffeineowl.com/graphics/2d/vectorial/cubic2quad01.html
		maxHypot2 := 432.0 * accuracy * accuracy
		p1x2 := Vec2(c.P1).Mul(3).Sub(Vec2(c.P0))
		p2x2 := Vec2(c.P2).Mul(3).Sub(Vec2(c.P3))
		err := p2x2.Sub(p1x2).Hypot2()
		n := max(int(math.Ceil(math.Sqrt(math.Cbrt(err/maxHypot2)))), 1)

		for i := range n {
			t0 := float64(i) / float64(n)
			t1 := float64(i+1) / float64(n)
			seg := c.Subsegment(t0, t1)
			p1x2 := Vec2(seg.P1).Mul(3).Sub(Vec2(seg.P0))
			p2x2 := Vec2(seg.P2).Mul(3).Sub(Vec2(seg.P3))
			result := QuadBez{seg.P0, Point(p1x2.Add(p2x2).Mul(1.0 / 4.0)), seg.P3}
			if !yield(CubicToQuadraticSegment{t0, t1, result}) {
				return
			}
		}
	}
}

func (c CubicBez) Subsegment(t0, t1 float64) CubicBez {
	p0 := c.Eval(t0)
	p3 := c.Eval(t1)
	d := c.Differentiate()
	scale := (t1 - t0) * (1.0 / 3.0)
	p1 := p0.Translate(Vec2(d.Eval(t0)).Mul(scale))
	p2 := p3.Translate(Vec2(d.Eval(t1)).Mul(scale).Negate())
	return CubicBez{p0, p1, p2, p3}
}

func (c CubicBez) SubsegmentCurve(t0, t1 float64) ParametricCurve {
	return c.Subsegment(t0, t1)
}

func (c CubicBez) Differentiate() QuadBez {
	return QuadBez{
		Point(c.P1.Sub(c.P0).Mul(3)),
		Point(c.P2.Sub(c.P1).Mul(3)),
		Point(c.P3.Sub(c.P2).Mul(3)),
	}
}

func arclenQuadratureCore(coeffs [][2]float64, dm Vec2, dm1 Vec2, dm2 Vec2) float64 {
	var sum float64
	for _, coeff := range coeffs {
		wi, xi := coeff[0], coeff[1]
		d := dm.Add(dm2.Mul(xi * xi))
		dpx := d.Add(dm1.Mul(xi)).Hypot()
		dmx := d.Sub(dm1.Mul(xi)).Hypot()
		sum += math.Sqrt(2.25) * wi * (dpx + dmx)
	}
	return sum
}

func (c CubicBez) Extrema() ([MaxExtrema]float64, int) {
	// two calls to oneCoord, up to 2 roots per call, for a total of 4 possible values.
	var out [MaxExtrema]float64
	var outN int
	oneCoord := func(d0, d1, d2 float64) {
		a := d0 - 2*d1 + d2
		b := 2 * (d1 - d0)
		c := d0
		roots, n := SolveQuadratic(c, b, a)
		for _, t := range roots[:n] {
			if t > 0.0 && t < 1.0 {
				out[outN] = t
				outN++
			}
		}
	}

	d0 := c.P1.Sub(c.P0)
	d1 := c.P2.Sub(c.P1)
	d2 := c.P3.Sub(c.P2)
	oneCoord(d0.X, d1.X, d2.X)
	oneCoord(d0.Y, d1.Y, d2.Y)
	sort.Float64s(out[:outN])
	return out, outN
}

// regularize preprocesses a cubic Bézier to ease numerical robustness.
//
// If the cubic Bézier segment has zero or near-zero derivatives, perturb
// the control points to make it easier to process (especially offset and
// stroke), avoiding numerical robustness problems.
func (c CubicBez) regularize(dimension float64) CubicBez {
	out := c
	// First step: if control point is too near the endpoint, nudge it away
	// along the tangent.
	dim2 := dimension * dimension
	if out.P0.DistanceSquared(out.P1) < dim2 {
		d02 := out.P0.DistanceSquared(out.P2)
		if d02 >= dim2 {
			// TODO: moderate if this would move closer to p3
			out.P1 = out.P0.Lerp(out.P2, math.Sqrt(dim2/d02))
		} else {
			out.P1 = out.P0.Lerp(out.P3, 1.0/3.0)
			out.P2 = out.P3.Lerp(out.P0, 1.0/3.0)
			return out
		}
	}
	if out.P3.DistanceSquared(out.P2) < dim2 {
		d13 := out.P1.DistanceSquared(out.P2)
		if d13 >= dim2 {
			// TODO: moderate if this would move closer to p0
			out.P2 = out.P3.Lerp(out.P1, math.Sqrt(dim2/d13))
		} else {
			out.P1 = out.P0.Lerp(out.P3, 1.0/3.0)
			out.P2 = out.P3.Lerp(out.P0, 1.0/3.0)
			return out
		}
	}
	if cuspType, ok := c.detectCusp(dimension); ok {
		d01 := out.P1.Sub(out.P0)
		d01h := d01.Hypot()
		d23 := out.P3.Sub(out.P2)
		d23h := d23.Hypot()
		switch cuspType {
		case CuspTypeLoop:
			out.P1 = out.P1.Translate(d01.Mul(dimension / d01h))
			out.P2 = out.P2.Translate(d23.Mul(dimension / d23h).Negate())
		case CuspTypeDoubleInflection:
			// Avoid making control distance smaller than dimension
			if d01h > 2.0*dimension {
				out.P1 = out.P1.Translate(d01.Mul(dimension / d01h).Negate())
			}
			if d23h > 2.0*dimension {
				out.P2 = out.P2.Translate(d23.Mul(dimension / d23h))
			}
		}
	}
	return out
}

// CustType defines the classification result for cusp detection.
type CuspType int

const (
	// Cusp is a loop.
	CuspTypeLoop CuspType = iota + 1
	// Cusp has two closely spaced inflection points.
	CuspTypeDoubleInflection
)

// detectCusp detects whether there is a cusp.
//
// It returns a cusp classification if there is a cusp with curvature greater than
// the reciprocal of the given dimension.
func (c CubicBez) detectCusp(dimension float64) (CuspType, bool) {
	d01 := c.P1.Sub(c.P0)
	d02 := c.P2.Sub(c.P0)
	d03 := c.P3.Sub(c.P0)
	d12 := c.P2.Sub(c.P1)
	d23 := c.P3.Sub(c.P2)
	det012 := d01.Cross(d02)
	det123 := d12.Cross(d23)
	det013 := d01.Cross(d03)
	det023 := d02.Cross(d03)
	if det012*det123 > 0.0 && det012*det013 < 0.0 && det012*det023 < 0.0 {
		q := c.Differentiate()
		// accuracy isn't used for quadratic nearest
		nearestDist, nearestT := q.Nearest(Point{}, 1e-9)
		// detect whether curvature at minimum derivative exceeds 1/dimension,
		// without division.
		d := q.Eval(nearestT)
		d2 := q.Differentiate().Eval(nearestT)
		cross := Vec2(d).Cross(Vec2(d2))
		if nearestDist*nearestDist*nearestDist <= cross*dimension*cross*dimension {
			a := 3.0*det012 + det023 - 2.0*det013
			b := -3.0*det012 + det013
			c := det012
			d := b*b - 4.0*a*c
			if d > 0.0 {
				return CuspTypeDoubleInflection, true
			} else {
				return CuspTypeLoop, true
			}
		}
	}
	return 0, false
}

func (c CubicBez) Transform(aff Affine) CubicBez {
	return CubicBez{
		P0: c.P0.Transform(aff),
		P1: c.P1.Transform(aff),
		P2: c.P2.Transform(aff),
		P3: c.P3.Transform(aff),
	}
}

// Nearest finds the nearest point, using subdivision.
func (c CubicBez) Nearest(pt Point, accuracy float64) (distSq, t float64) {
	var bestR option[float64]
	bestT := 0.0
	for qq := range c.Quadratics(accuracy) {
		t0, t1, q := qq.Start, qq.End, qq.Segment
		qDistSq, qT := q.Nearest(pt, accuracy)
		if !bestR.isSet || qDistSq < bestR.value {
			bestT = t0 + qT*(t1-t0)
			bestR.set(qDistSq)
		}
	}
	return bestR.value, bestT
}

func (c CubicBez) SignedArea() float64 {
	v := c.P0.X*(6.0*c.P1.Y+3.0*c.P2.Y+c.P3.Y) +
		3.0*(c.P1.X*(-2.0*c.P0.Y+c.P2.Y+c.P3.Y)-c.P2.X*(c.P0.Y+c.P1.Y-2.0*c.P3.Y)) -
		c.P3.X*(c.P0.Y+3.0*c.P1.Y+6.0*c.P2.Y)
	return v * (1.0 / 20.0)
}

func (c CubicBez) Tangents() (Vec2, Vec2) {
	const epsilon = 1e-12
	d01 := c.P1.Sub(c.P0)
	var d0, d1 Vec2
	if d01.Hypot2() > epsilon {
		d0 = d01
	} else {
		d02 := c.P2.Sub(c.P0)
		if d02.Hypot2() > epsilon {
			d0 = d02
		} else {
			d0 = c.P3.Sub(c.P0)
		}
	}
	d23 := c.P3.Sub(c.P2)
	if d23.Hypot2() > epsilon {
		d1 = d23
	} else {
		d13 := c.P3.Sub(c.P1)
		if d13.Hypot2() > epsilon {
			d1 = d13
		} else {
			d1 = c.P3.Sub(c.P0)
		}
	}
	return d0, d1
}

func (c CubicBez) Seg() PathSegment {
	return PathSegment{Kind: CubicKind, P0: c.P0, P1: c.P1, P2: c.P2, P3: c.P3}
}

func (c CubicBez) IntersectLine(line Line) ([3]LineIntersection, int) {
	const epsilon = 1e-9
	p0 := line.P0
	p1 := line.P1
	dx := p1.X - p0.X
	dy := p1.Y - p0.Y

	// The basic technique here is to determine x and y as a cubic polynomial
	// as a function of t. Then plug those values into the line equation for the
	// probe line (giving a sort of signed distance from the probe line) and solve
	// that for t.
	px0, px1, px2, px3 := cubicBezCoefficients(c.P0.X, c.P1.X, c.P2.X, c.P3.X)
	py0, py1, py2, py3 := cubicBezCoefficients(c.P0.Y, c.P1.Y, c.P2.Y, c.P3.Y)
	c0 := dy*(px0-p0.X) - dx*(py0-p0.Y)
	c1 := dy*px1 - dx*py1
	c2 := dy*px2 - dx*py2
	c3 := dy*px3 - dx*py3
	invlen2 := 1.0 / (dx*dx + dy*dy)
	ts, n := SolveCubic(c0, c1, c2, c3)
	var ret [3]LineIntersection
	var retN int
	for _, t := range ts[:n] {
		if t >= -epsilon && t <= 1+epsilon {
			x := px0 + t*px1 + t*t*px2 + t*t*t*px3
			y := py0 + t*py1 + t*t*py2 + t*t*t*py3
			u := ((x-p0.X)*dx + (y-p0.Y)*dy) * invlen2
			if u >= 0.0 && u <= 1.0 {
				ret[retN] = LineIntersection{u, t}
				retN++
			}
		}
	}
	return ret, retN
}

// Return polynomial coefficients given cubic bezier coordinates.
func cubicBezCoefficients(x0, x1, x2, x3 float64) (_, _, _, _ float64) {
	p0 := x0
	p1 := 3.0*x1 - 3.0*x0
	p2 := 3.0*x2 - 6.0*x1 + 3.0*x0
	p3 := x3 - 3.0*x2 + 3.0*x1 - x0
	return p0, p1, p2, p3
}

// Inflections returns the inflection points.
//
// The function returns up to two inflection points in the first return
// parameter, with the second parameter specifying the number of points
// returned.
func (cb CubicBez) Inflections() ([2]float64, int) {
	a := cb.P1.Sub(cb.P0)
	b := cb.P2.Sub(cb.P1).Sub(a)
	c := cb.P3.Sub(cb.P0).Sub(cb.P2.Sub(cb.P1).Mul(3))
	nums, n := SolveQuadratic(a.Cross(b), a.Cross(c), b.Cross(c))
	var out [2]float64
	var outN int
	for _, num := range nums[:n] {
		if num >= 0 && num <= 1 {
			out[outN] = num
			outN++
		}
	}
	return out, outN
}

// ApproxQuadSpline returns a quadratic B-spline approximating this cubic
// Bézier.
//
// Returns false if no suitable approximation is found within the given
// tolerance.
func (cb CubicBez) ApproxQuadSpline(accuracy float64) (QuadBSpline, bool) {
	for i := range maxSplineSplit {
		if spline, ok := cb.approxQuadSplineN(i+1, accuracy); ok {
			return spline, true
		}
	}
	return nil, false
}

// Approximate a cubic curve with a quadratic spline of n curves
func (c CubicBez) approxQuadSplineN(n int, accuracy float64) (QuadBSpline, bool) {
	if n == 1 {
		qb, ok := c.tryApproxQuadratic(accuracy)
		if !ok {
			return nil, false
		}
		return QuadBSpline{qb.P0, qb.P1, qb.P2}, true
	}

	cubicsNext_, cubicsDone := iter.Pull(c.splitIntoN(n))
	defer cubicsDone()
	cubicsNext := func() CubicBez {
		v, ok := cubicsNext_()
		if !ok {
			panic("unreachable")
		}
		return v
	}

	nextCubic := cubicsNext()
	nextQ1 := nextCubic.approxQuadControl(0.0)
	q2 := c.P0
	var d1 Vec2
	spline := []Point{c.P0, nextQ1}
	for i := 1; i <= n; i++ {
		currentCubic := nextCubic
		q0 := q2
		q1 := nextQ1
		if i < n {
			nextCubic = cubicsNext()
			nextQ1 = nextCubic.approxQuadControl(float64(i) / float64(n-1))

			spline = append(spline, nextQ1)
			q2 = q1.Midpoint(nextQ1)
		} else {
			q2 = currentCubic.P3
		}
		d0 := d1
		d1 = Vec2(q2).Sub(Vec2(currentCubic.P3))

		if d1.Hypot() > accuracy ||
			!(CubicBez{
				Point(d0),
				// XXX is negate fine
				q0.Lerp(q1, 2.0/3.0).Translate(Vec2(currentCubic.P1).Negate()),
				q2.Lerp(q1, 2.0/3.0).Translate(Vec2(currentCubic.P2).Negate()),
				Point(d1),
			}.fitsInside(accuracy)) {
			return nil, false
		}
	}
	spline = append(spline, c.P3)
	return QuadBSpline(spline), true
}

func (c CubicBez) splitIntoN(n int) iter.Seq[CubicBez] {
	// We special case some values of n so that we produce the same results as Kurbo,
	// which in turn produces the same results as cu2qu.

	switch n {
	case 1:
		return func(yield func(CubicBez) bool) { yield(c) }
	case 2:
		return func(yield func(CubicBez) bool) {
			l, r := c.Subdivide()
			_ = yield(l) && yield(r)
		}
	case 3:
		return func(yield func(CubicBez) bool) {
			l, m, r := c.subdivide3()
			_ = yield(l) && yield(m) && yield(r)
		}
	case 4:
		return func(yield func(CubicBez) bool) {
			l, r := c.Subdivide()
			ll, lr := l.Subdivide()
			rl, rr := r.Subdivide()
			_ = yield(ll) &&
				yield(lr) &&
				yield(rl) &&
				yield(rr)
		}
	case 6:
		return func(yield func(CubicBez) bool) {
			l, r := c.Subdivide()
			l1, l2, l3 := l.subdivide3()
			r1, r2, r3 := r.subdivide3()
			_ = yield(l1) &&
				yield(l2) &&
				yield(l3) &&
				yield(r1) &&
				yield(r2) &&
				yield(r3)
		}
	}

	return func(yield func(CubicBez) bool) {
		a, b, c, d := c.parameters()
		dt := 1.0 / float64(n)
		delta2 := dt * dt
		delta3 := dt * delta2
		for i := range n {
			t1 := float64(i) * dt
			t1_2 := t1 * t1
			a1 := a.Mul(delta3)
			// Note that we don't simplify a * 3.0 * t1 to a * (3.0 * t1) because
			// multiplication isn't associative in floating point math, and we want our
			// results to match cu2qu's.
			b1 := a.Mul(3.0).Mul(t1).Add(b).Mul(delta2)
			c1 := b.Mul(2.0).Mul(t1).Add(c).Add(a.Mul(3.0).Mul(t1_2)).Mul(dt)
			d1 := a.Mul(t1).Mul(t1_2).Add(b.Mul(t1_2)).Add(c.Mul(t1)).Add(d)

			// Rust port of cu2qu [calc_cubic_points](https://github.com/fonttools/fonttools/blob/3b9a73ff8379ab49d3ce35aaaaf04b3a7d9d1655/Lib/fontTools/cu2qu/cu2qu.py#L63-L68)
			p0 := Point(d1)
			p1 := Point(c1.Div(3)).Translate(d1)
			p2 := Point(b1.Add(c1).Div(3)).Translate(Vec2(p1))
			p3 := Point(a1.Add(d1).Add(c1).Add(b1))
			result := CubicBez{p0, p1, p2, p3}
			if !yield(result) {
				break
			}
		}
	}
}

func (cb CubicBez) parameters() (Vec2, Vec2, Vec2, Vec2) {
	c := cb.P1.Sub(cb.P0).Mul(3.0)
	b := cb.P2.Sub(cb.P1).Mul(3.0).Sub(c)
	d := Vec2(cb.P0)
	a := Vec2(cb.P3).Sub(d).Sub(c).Sub(b)
	return a, b, c, d
}

func (c CubicBez) subdivide3() (CubicBez, CubicBez, CubicBez) {
	p0, p1, p2, p3 :=
		Vec2(c.P0),
		Vec2(c.P1),
		Vec2(c.P2),
		Vec2(c.P3)

	// The original Python cu2qu code does not use the division operator to divide by
	// 27 but instead uses multiplication by the reciprocal 1 / 27. We want to match
	// it exactly to avoid any floating point differences.
	//
	// Source: https://github.com/fonttools/fonttools/blob/85c80be/Lib/fontTools/cu2qu/cu2qu.py#L215-L218
	// See also: https://github.com/linebender/kurbo/issues/272
	mid1 := Point(p0.Mul(8).Add(p1.Mul(12)).Add(p2.Mul(6)).Add(p3).Mul(1.0 / 27.0))
	deriv1 := p3.Add(p2.Mul(3)).Sub(p0.Mul(4)).Mul(1.0 / 27.0)
	mid2 := Point(p0.Add(p1.Mul(6)).Add(p2.Mul(12)).Add(p3.Mul(8)).Mul(1.0 / 27.0))
	deriv2 := p3.Mul(4).Sub(p1.Mul(3)).Sub(p0).Mul(1.0 / 27.0)

	left := CubicBez{
		c.P0,
		Point(p0.Mul(2).Add(p1).Div(3)),
		// XXX is mid2 + (-deriv2) the same as mid2 - deriv2, under float math?
		mid1.Translate(deriv1.Negate()),
		mid1,
	}
	// XXX is mid2 + (-deriv2) the same as mid2 - deriv2, under float math?
	mid := CubicBez{mid1, mid1.Translate(deriv1), mid2.Translate(deriv2.Negate()), mid2}
	right := CubicBez{
		mid2,
		mid2.Translate(deriv2),
		Point(p2.Add(p3.Mul(2)).Div(3)),
		c.P3,
	}
	return left, mid, right
}

// Does this curve fit inside the given distance from the origin?
//
// Rust port of cu2qu [cubic_farthest_fit_inside](https://github.com/fonttools/fonttools/blob/3b9a73ff8379ab49d3ce35aaaaf04b3a7d9d1655/Lib/fontTools/cu2qu/cu2qu.py#L281)
func (c CubicBez) fitsInside(distance float64) bool {
	if Vec2(c.P2).Hypot() <= distance && Vec2(c.P1).Hypot() <= distance {
		return true
	}
	mid := Vec2(c.P0).Add(Vec2(c.P1).Add(Vec2(c.P2)).Mul(3)).Add(Vec2(c.P3)).Mul(0.125)
	if mid.Hypot() > distance {
		return false
	}
	// Split in two. Note that cu2qu here uses a 3/8 subdivision. I don't know why.
	left, right := c.Subdivide()
	return left.fitsInside(distance) && right.fitsInside(distance)
}

func (c CubicBez) approxQuadControl(t float64) Point {
	p1 := c.P0.Translate(c.P1.Sub(c.P0).Mul(1.5))
	p2 := c.P3.Translate(c.P2.Sub(c.P3).Mul(1.5))
	return p1.Lerp(p2, t)
}

// Approximate a cubic with a single quadratic
//
// Returns a quadratic approximating the given cubic that maintains
// endpoint tangents if that is within tolerance, or false otherwise.
func (c CubicBez) tryApproxQuadratic(accuracy float64) (QuadBez, bool) {
	q1, ok := Line{c.P0, c.P1}.CrossingPoint(Line{c.P2, c.P3})
	if !ok {
		return QuadBez{}, false
	}

	c1 := c.P0.Lerp(q1, 2.0/3.0)
	c2 := c.P3.Lerp(q1, 2.0/3.0)
	if !(CubicBez{
		Point{},
		// XXX negate ok?
		c1.Translate(Vec2(c.P1).Negate()),
		c2.Translate(Vec2(c.P2).Negate()),
		Point{},
	}.fitsInside(accuracy)) {
		return QuadBez{}, false
	}
	return QuadBez{c.P0, q1, c.P3}, true
}

// CubicsToQuadraticSplines converts multiple cubic Bézier curves to quadratic
// splines. It ensures that the resulting splines have the same number of
// control points.
func CubicsToQuadraticSplines(curves []CubicBez, accuracy float64) ([]QuadBSpline, bool) {
	result := make([]QuadBSpline, 0, len(curves))
outer:
	for i := 1; i <= maxSplineSplit; i++ {
		result = result[:0]

		for _, curve := range curves {
			spline, ok := curve.approxQuadSplineN(i, accuracy)
			if !ok {
				continue outer
			}
			result = append(result, spline)
		}

		return result, true
	}
	return nil, false
}
