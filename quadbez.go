package curve

import (
	"iter"
	"math"
)

var _ Shape = QuadBez{}
var _ ParametricCurve = QuadBez{}

type QuadBez struct {
	P0 Point
	P1 Point
	P2 Point
}

func (q QuadBez) BoundingBox() Rect {
	return BoundingBox(q)
}

func (q QuadBez) PathElements(tolerance float64) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		_ = yield(MoveTo(q.P0)) &&
			yield(QuadTo(q.P1, q.P2))
	}
}

func (q QuadBez) Perimeter(accuracy float64) float64 {
	return q.Arclen(accuracy)
}

// / Raise the order by 1.
// /
// / Returns a cubic Bézier segment that exactly represents this quadratic.
func (q QuadBez) Raise() CubicBez {
	return CubicBez{
		q.P0,
		q.P0.Translate(q.P1.Sub(q.P0).Mul(2.0 / 3.0)),
		q.P2.Translate(q.P1.Sub(q.P2).Mul(2.0 / 3.0)),
		q.P2,
	}
}

func (q QuadBez) IsInf() bool {
	return q.P0.IsInf() || q.P1.IsInf() || q.P2.IsInf()
}

func (q QuadBez) IsNaN() bool {
	return q.P0.IsNaN() || q.P1.IsNaN() || q.P2.IsNaN()
}

// Arclen returns the arclength of the quadratic Bézier segment.
//
// This computation is based on an analytical formula. Since that formula suffers
// from numerical instability when the curve is very close to a straight line, we
// detect that case and fall back to Legendre-Gauss quadrature.
//
// Overall accuracy should be better than 1e-13 over the entire range.
func (q QuadBez) Arclen(accuracy float64) float64 {
	d2 := Vec2(q.P0).Sub(Vec2(q.P1).Mul(2)).Add(Vec2(q.P2))
	a := d2.Hypot2()
	d1 := q.P1.Sub(q.P0)
	c := d1.Hypot2()
	if a < 5e-4*c {
		// This case happens for nearly straight Béziers.
		//
		// Calculate arclength using Legendre-Gauss quadrature using formula from Behdad
		// in https://github.com/Pomax/BezierInfo-2/issues/77
		v0 := Vec2(q.P0).Mul(-0.492943519233745).
			Add(Vec2(q.P1).Mul(0.430331482911935)).
			Add(Vec2(q.P2).Mul(0.0626120363218102)).
			Hypot()
		v1 := q.P2.Sub(q.P0).Mul(0.4444444444444444).Hypot()
		v2 := Vec2(q.P0).Mul(-0.0626120363218102).
			Sub(Vec2(q.P1).Mul(0.430331482911935)).
			Add(Vec2(q.P2).Mul(0.492943519233745)).
			Hypot()
		return v0 + v1 + v2
	}
	b := 2.0 * d2.Dot(d1)

	sabc := math.Sqrt(a + b + c)
	a2 := math.Pow(a, -0.5)
	a32 := a2 * a2 * a2
	c2 := 2.0 * math.Sqrt(c)
	baC2 := b*a2 + c2

	v0 := 0.25*a2*a2*b*(2.0*sabc-c2) + sabc
	// TODO: justify and fine-tune this exact constant.
	if baC2 < 1e-13 {
		// This case happens for Béziers with a sharp kink.
		return v0
	} else {
		return v0 + 0.25*a32*(4.0*c*a-b*b)*math.Log(((2.0*a+b)*a2+2.0*sabc)/baC2)
	}
}

func (q QuadBez) Eval(t float64) Point {
	mt := 1.0 - t
	a := Vec2(q.P0).Mul(mt * mt)
	b := Vec2(q.P1).Mul(mt * 2.0)
	c := Vec2(q.P2).Mul(t)
	d := b.Add(c)
	return Point(a.Add(d.Mul(t)))
}

func (q QuadBez) Subdivide() (QuadBez, QuadBez) {
	pm := q.Eval(0.5)
	return QuadBez{q.P0, q.P0.Midpoint(q.P1), pm},
		QuadBez{pm, q.P1.Midpoint(q.P2), q.P2}
}

func (q QuadBez) SubdivideCurve() (ParametricCurve, ParametricCurve) {
	return q.Subdivide()
}

func (q QuadBez) Subsegment(t0 float64, t1 float64) QuadBez {
	p0 := q.Eval(t0)
	p2 := q.Eval(t1)
	p1 := p0.Translate(q.P1.Sub(q.P0).Lerp(q.P2.Sub(q.P1), t0).Mul(t1 - t0))
	return QuadBez{p0, p1, p2}
}

func (q QuadBez) SubsegmentCurve(t0 float64, t1 float64) ParametricCurve {
	return q.Subsegment(t0, t1)
}

func (q QuadBez) Differentiate() Line {
	return Line{
		Point(q.P1.Sub(q.P0).Mul(2)),
		Point(q.P2.Sub(q.P1).Mul(2)),
	}
}

func (q QuadBez) Start() Point {
	return q.P0
}

func (q QuadBez) End() Point {
	return q.P2
}

func (q QuadBez) Extrema() ([MaxExtrema]float64, int) {
	// Finding the extrema of a quadratic bezier means finding the roots in the
	// quadratic's first derivative, which is a line.

	var out [MaxExtrema]float64
	var outN int
	d0 := q.P1.Sub(q.P0)
	d1 := q.P2.Sub(q.P1)
	dd := d1.Sub(d0)
	if dd.X != 0.0 {
		t := -d0.X / dd.X
		if t > 0.0 && t < 1.0 {
			out[outN] = t
			outN++
		}
	}
	if dd.Y != 0 {
		t := -d0.Y / dd.Y
		if t > 0.0 && t < 1.0 {
			out[outN] = t
			outN++
			if outN == 2 && out[0] > t {
				out[0], out[1] = out[1], out[0]
			}
		}
	}
	return out, outN
}

func (q QuadBez) Nearest(pt Point, accuracy float64) (distSq, outT float64) {
	/// Find the nearest point, using analytical algorithm based on cubic root finding.
	evalT := func(pt Point, tBest *float64, rBest *option[float64], t float64, p0 Point) {
		r := p0.Sub(pt).Hypot2()
		if !rBest.isSet || r < rBest.value {
			rBest.set(r)
			*tBest = t
		}
	}
	tryT := func(
		q *QuadBez,
		pt Point,
		tBest *float64,
		rBest *option[float64],
		t float64,
	) bool {
		if !(t >= 0.0 && t <= 1.0) {
			return true
		}
		evalT(pt, tBest, rBest, t, q.Eval(t))
		return false
	}
	d0 := q.P1.Sub(q.P0)
	d1 := Vec2(q.P0).Add(Vec2(q.P2)).Sub(Vec2(q.P1).Mul(2.0))
	d := q.P0.Sub(pt)
	c0 := d.Dot(d0)
	c1 := 2.0*d0.Hypot2() + d.Dot(d1)
	c2 := 3.0 * d1.Dot(d0)
	c3 := d1.Hypot2()
	roots, n := SolveCubic(c0, c1, c2, c3)
	var rBest option[float64]
	tBest := 0.0
	needEnds := n == 0

	for _, t := range roots[:n] {
		b := tryT(&q, pt, &tBest, &rBest, t)
		if b {
			needEnds = true
		}
	}
	if needEnds {
		evalT(pt, &tBest, &rBest, 0.0, q.P0)
		evalT(pt, &tBest, &rBest, 1.0, q.P2)
	}

	return rBest.value, tBest
}

func (q QuadBez) Transform(aff Affine) QuadBez {
	return QuadBez{
		P0: q.P0.Transform(aff),
		P1: q.P1.Transform(aff),
		P2: q.P2.Transform(aff),
	}
}

func (q QuadBez) SignedArea() float64 {
	v := q.P0.X*(2.0*q.P1.Y+q.P2.Y) +
		2.0*(q.P1.X*(q.P2.Y-q.P0.Y)) -
		q.P2.X*(q.P0.Y+2.0*q.P1.Y)
	return v * (1.0 / 6.0)
}

func (q QuadBez) Tangents() (Vec2, Vec2) {
	const epsilon = 1e-12
	d01 := q.P1.Sub(q.P0)
	var d0, d1 Vec2
	if d01.Hypot2() > epsilon {
		d0 = d01
	} else {
		d0 = q.P2.Sub(q.P0)
	}
	d12 := q.P2.Sub(q.P1)
	if d12.Hypot2() > epsilon {
		d1 = d12
	} else {
		d1 = q.P2.Sub(q.P0)
	}
	return d0, d1
}

// An approximation to $\int (1 + 4x^2) ^ -0.25 dx$
//
// This is used for flattening curves.
func approxParabolaIntegral(x float64) float64 {
	const d = 0.67
	return x / (1.0 - d + math.Sqrt(math.Sqrt(math.Pow(d, 4)+0.25*x*x)))
}

// An approximation to the inverse parabola integral.
func approxParabolaInvIntegral(x float64) float64 {
	const b = 0.39
	return x * (1.0 - b + math.Sqrt(b*b+0.25*x*x))
}

// Maps a value from 0..1 to 0..1.
func (q QuadBez) determineSubdivT(params *flattenParams, x float64) float64 {
	a := params.a0 + (params.a2-params.a0)*x
	u := approxParabolaInvIntegral(a)
	return (u - params.u0) * params.uscale
}

// / Estimate the number of subdivisions for flattening.
func (q QuadBez) estimateSubdiv(sqrtTol float64) flattenParams {
	// Determine transformation to $y = x^2$ parabola.
	d01 := q.P1.Sub(q.P0)
	d12 := q.P2.Sub(q.P1)
	dd := d01.Sub(d12)
	cross := q.P2.Sub(q.P0).Cross(dd)
	x0 := d01.Dot(dd) * (1.0 / cross)
	x2 := d12.Dot(dd) * (1.0 / cross)
	scale := math.Abs(cross / (dd.Hypot() * (x2 - x0)))

	// Compute number of subdivisions needed.
	a0 := approxParabolaIntegral(x0)
	a2 := approxParabolaIntegral(x2)
	var val float64
	if !math.IsInf(scale, 0) {
		da := math.Abs(a2 - a0)
		sqrtScale := math.Sqrt(scale)
		if math.Signbit(x0) == math.Signbit(x2) {
			val = da * sqrtScale
		} else {
			// Handle cusp case (segment contains curvature maximum)
			xmin := sqrtTol / sqrtScale
			val = sqrtTol * da / approxParabolaIntegral(xmin)
		}
	}
	u0 := approxParabolaInvIntegral(a0)
	u2 := approxParabolaInvIntegral(a2)
	uscale := 1.0 / (u2 - u0)
	return flattenParams{
		a0,
		a2,
		u0,
		uscale,
		val,
	}
}

type flattenParams struct {
	a0     float64
	a2     float64
	u0     float64
	uscale float64
	/// The number of subdivisions * 2 * sqrtTol.
	val float64
}

func (q QuadBez) Seg() PathSegment {
	return PathSegment{Kind: QuadKind, P0: q.P0, P1: q.P1, P2: q.P2}
}

func (q QuadBez) IntersectLine(line Line) ([3]LineIntersection, int) {
	const epsilon = 1e-9
	p0 := line.P0
	p1 := line.P1
	dx := p1.X - p0.X
	dy := p1.Y - p0.Y

	// The basic technique here is to determine x and y as a quadratic polynomial
	// as a function of t. Then plug those values into the line equation for the
	// probe line (giving a sort of signed distance from the probe line) and solve
	// that for t.
	px0, px1, px2 := quadBezCoefficients(q.P0.X, q.P1.X, q.P2.X)
	py0, py1, py2 := quadBezCoefficients(q.P0.Y, q.P1.Y, q.P2.Y)
	c0 := dy*(px0-p0.X) - dx*(py0-p0.Y)
	c1 := dy*px1 - dx*py1
	c2 := dy*px2 - dx*py2
	invlen2 := 1.0 / (dx*dx + dy*dy)
	ts, n := SolveQuadratic(c0, c1, c2)
	var ret [3]LineIntersection
	var retN int
	for _, t := range ts[:n] {
		if t >= -epsilon && t <= 1+epsilon {
			x := px0 + t*px1 + t*t*px2
			y := py0 + t*py1 + t*t*py2
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
func quadBezCoefficients(x0, x1, x2 float64) (_, _, _ float64) {
	p0 := x0
	p1 := 2.0*x1 - 2.0*x0
	p2 := x2 - 2.0*x1 + x0
	return p0, p1, p2
}
