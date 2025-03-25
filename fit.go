package curve

import (
	"fmt"
	"iter"
	"math"
)

// As described in [Simplifying Bézier paths], strictly optimizing for
// Fréchet distance can create bumps. The problem is curves with long
// control arms (distance from the control point to the corresponding
// endpoint). We mitigate that by applying a penalty as a multiplier to
// the measured error (approximate Fréchet distance). This is ReLU-like,
// with a value of 1.0 below the elbow, and a given slope above it. The
// values here have been determined empirically to give good results.
//
// [Simplifying Bézier paths]: https://raphlinus.github.io/curves/2023/04/18/bezpath-simplify.html
const (
	dPenaltyElbow = 0.65
	dPenaltySlope = 2.0
)

const numSamples = 20

// FittableCurve describes curves in a way useful for curve fitting. Curves that
// implement this interface can be used as source curves in [FitToBezPath] and
// [FitToBezPathOpt].
//
// The interface can represent source curves with cusps and corners, though if
// the corners are known in advance, it may be better to run curve fitting on
// subcurves bounded by the corners.
//
// The interface primarily works by sampling the source curve and computing the
// position and derivative at each sample. Those derivatives are then used for
// multiple sub-tasks, including ensuring G1 continuity at subdivision points,
// computing the area and moment of the curve for curve fitting, and casting
// rays for evaluation of a distance metric to test accuracy.
//
// A major motivation is computation of offset curves, which often have cusps,
// but the presence and location of those cusps is not generally known. It is
// also intended for conversion between curve types (for example, piecewise
// Euler spiral or NURBS), and distortion effects such as perspective transform.
//
// Note general similarities to [ParametricCurve] but also important
// differences. Instead of separate methods for evaluating the curve and its
// derivative, there is a single SamplePtDeriv method, which can be more
// efficient and also handles cusps more robustly. Also, there is no method for
// subsegmenting, as that is not needed and would be annoying to implement.
type FittableCurve interface {
	// SamplePtTangent evaluates the curve and its tangent at parameter t.
	//
	// For a regular curve (one not containing a cusp or corner), the
	// derivative is a good choice for the tangent vector and the sign
	// parameter can be ignored. Otherwise, the sign parameter selects which
	// side of the discontinuity the tangent will be sampled from.
	//
	// Generally, t is in the range [0, 1].
	SamplePtTangent(t float64, sign float64) CurveFitSample
	// SamplePtDeriv evaluates the point and derivative at parameter t.
	//
	// In curves with cusps, the derivative can go to zero.
	SamplePtDeriv(t float64) (Point, Vec2)
	// BreakCusp findf a cusp or corner within the given range.
	//
	// If the range contains a corner or cusp, return it. If there is more
	// than one such discontinuity, any can be reported, as the function will
	// be called repeatedly after subdivision of the range.
	//
	// Do not report cusps at the endpoints of the range, as this may cause
	// potentially infinite subdivision. In particular, when a cusp is reported
	// and this method is called on a subdivided range bounded by the reported
	// cusp, then the subsequent call should not report a cusp there.
	//
	// The definition of what exactly constitutes a cusp is somewhat loose.
	// If a cusp is missed, then the curve fitting algorithm will attempt to
	// fit the curve with a smooth curve, which is generally not a disaster but
	// will usually result in more subdivision. Conversely, it might be useful
	// to report near-cusps, specifically points of curvature maxima where the
	// curvature is large but still mathematically finite.
	BreakCusp(start, end float64) (float64, bool)
}

type MomentIntegraler interface {
	MomentIntegrals(start, end float64) (float64, float64, float64)
}

// MomentIntegrals computes moment integrals.
//
// This function computes the integrals of y dx, x y dx, and y^2 dx over the
// length of this curve. From these integrals it is fairly straightforward
// to derive the moments needed for curve fitting.
//
// By default it uses quadrature integration with Green's theorem, in terms of
// samples evaluated with [FittableCurve.SamplePtDeriv]. If pcf implements
// [MomentIntegraler], then this function defers to it.
func MomentIntegrals(pcf FittableCurve, start, end float64) (float64, float64, float64) {
	if pcf, ok := pcf.(MomentIntegraler); ok {
		return pcf.MomentIntegrals(start, end)
	}
	t0 := 0.5 * (start + end)
	dt := 0.5 * (end - start)

	var a, x, y float64
	for _, v := range gaussLegendreCoeffs16 {
		wi, xi := v[0], v[1]
		t := t0 + xi*dt
		p, d := pcf.SamplePtDeriv(t)
		a_ := wi * d.X * p.Y
		x_ := p.X * a_
		y_ := p.Y * a_
		a += a_
		x += x_
		y += y_
	}
	return a * dt, x * dt, y * dt
}

// FitToCubic fits a single cubic to a range of the source curve.
//
// Returns the cubic segment and the square of the error. Returns false if no
// fitting cubic could be computed.
func FitToCubic(
	source FittableCurve,
	rangeStart float64,
	rangeEnd float64,
	accuracy float64,
) (CubicBez, float64, bool) {
	start := source.SamplePtTangent(rangeStart, 1.0)
	end := source.SamplePtTangent(rangeEnd, -1.0)
	d := end.Point.Sub(start.Point)
	chord2 := d.Hypot2()
	acc2 := accuracy * accuracy
	if chord2 <= acc2 {
		// Special case very short chords; try to fit a line.
		return tryFitLine(source, accuracy, rangeStart, rangeEnd, start.Point, end.Point)
	}
	th := d.Angle()
	mod2Pi := func(th float64) float64 {
		thScaled := th * (1.0 / math.Pi) * 0.5
		return math.Pi * 2.0 * (thScaled - math.Round(thScaled))
	}
	th0 := mod2Pi(start.Tangent.Angle() - th)
	th1 := mod2Pi(th - end.Tangent.Angle())

	area, x, y := MomentIntegrals(source, rangeStart, rangeEnd)
	x0, y0 := start.Point.X, start.Point.Y
	dx, dy := d.X, d.Y
	// Subtract off area of chord
	area -= dx * (y0 + 0.5*dy)
	// area is signed area of closed curve segment.
	// This quantity is invariant to translation and rotation.

	// Subtract off moment of chord
	dy3 := dy * (1.0 / 3.0)
	x -= dx * (x0*y0 + 0.5*(x0*dy+y0*dx) + dy3*dx)
	y -= dx * (y0*y0 + y0*dy + dy3*dy)
	// Translate start point to origin; convert raw integrals to moments.
	x -= x0 * area
	y = 0.5*y - y0*area
	// Rotate into place (this also scales up by chordlength for efficiency).
	moment := d.X*x + d.Y*y
	// moment is the chordlength times the x moment of the curve translated
	// so its start point is on the origin, and rotated so its end point is on the
	// x axis.

	chord2Inv := 1.0 / chord2
	unitArea := area * chord2Inv
	mx := moment * (chord2Inv * chord2Inv)
	// unitArea is signed area scaled to unit chord; mx is scaled x moment

	chord := math.Sqrt(chord2)
	aff := Translate(Vec2(start.Point)).Mul(Rotate(th)).Mul(Scale(chord, chord))
	curveDist := curveDistFromCurve(source, rangeStart, rangeEnd)
	var bestC option[CubicBez]
	var bestErr2 option[float64]
	fits, fitsN := cubicFit(th0, th1, unitArea, mx)
	for _, cfit := range fits[:fitsN] {
		cand := cfit.cbez
		d0 := cfit.d0
		d1 := cfit.d1
		c := cand.Transform(aff)
		if err2, ok := curveDist.evalDist(source, c, acc2); ok {
			scaleF := func(d float64) float64 {
				return 1.0 + max(d-dPenaltyElbow, 0.0)*dPenaltySlope
			}
			scale := math.Pow(max(scaleF(d0), scaleF(d1)), 2)
			err2 := err2 * scale
			if err2 < acc2 && (!bestErr2.isSet || err2 < bestErr2.value) {
				bestC.set(c)
				bestErr2.set(err2)
			}
		}
	}
	if bestC.isSet && bestErr2.isSet {
		return bestC.value, bestErr2.value, true
	} else {
		return CubicBez{}, 0, false
	}
}

// cubicFit returns curves matching area and moment, given unit chord.
func cubicFit(th0 float64, th1 float64, area float64, mx float64) ([4]struct {
	cbez CubicBez
	d0   float64
	d1   float64
}, int) {
	// Note: maybe we want to take unit vectors instead of angle? Shouldn't
	// matter much either way though.
	s0, c0 := math.Sincos(th0)
	s1, c1 := math.Sincos(th1)
	a4 := -9.0 *
		c0 *
		(((2.0*s1*c1*c0+s0*(2.0*c1*c1-1.0))*c0-2.0*s1*c1)*c0 -
			c1*c1*s0)
	a3 := 12.0 *
		((((c1*(30.0*area*c1-s1)-15.0*area)*c0+2.0*s0-
			c1*s0*(c1+30.0*area*s1))*
			c0+
			c1*(s1-15.0*area*c1))*
			c0 -
			s0*c1*c1)
	a2 := 12.0 *
		((((70.0*mx+15.0*area)*s1*s1+c1*(9.0*s1-70.0*c1*mx-5.0*c1*area))*
			c0-
			5.0*s0*s1*(3.0*s1-4.0*c1*(7.0*mx+area)))*
			c0 -
			c1*(9.0*s1-70.0*c1*mx-5.0*c1*area))
	a1 := 16.0 *
		(((12.0*s0-5.0*c0*(42.0*mx-17.0*area))*s1-
			70.0*c1*(3.0*mx-area)*s0-
			75.0*c0*c1*area*area)*
			s1 -
			75.0*c1*c1*area*area*s0)
	a0 := 80.0 * s1 * (42.0*s1*mx - 25.0*area*(s1-c1*area))
	// TODO: "roots" is not a good name for this variable, as it also contains
	// the real part of complex conjugate pairs.
	roots := make([]float64, 0, 4)
	const epsilon = 1e-12
	if math.Abs(a4) > epsilon {
		a := a3 / a4
		b := a2 / a4
		c := a1 / a4
		d := a0 / a4
		if quads, ok := FactorQuarticInner(a, b, c, d, false); ok {
			for _, quad := range quads {
				qc1, qc0 := quad[0], quad[1]
				qroots, n := SolveQuadratic(qc0, qc1, 1.0)
				if n == 0 {
					// Real part of pair of complex roots
					roots = append(roots, -0.5*qc1)
				} else {
					roots = append(roots, qroots[:n]...)
				}
			}
		}
	} else if math.Abs(a3) > epsilon {
		qroots, n := SolveCubic(a0, a1, a2, a3)
		roots = append(roots, qroots[:n]...)
	} else if math.Abs(a2) > epsilon || math.Abs(a1) > epsilon || math.Abs(a0) > epsilon {
		qroots, n := SolveQuadratic(a0, a1, a2)
		roots = append(roots, qroots[:n]...)
	} else {
		return [4]struct {
			cbez CubicBez
			d0   float64
			d1   float64
		}{{
			CubicBez{
				Pt(0.0, 0.0),
				Pt(1.0/3.0, 0.0),
				Pt(2.0/3.0, 0.0),
				Pt(1.0, 0.0),
			},
			1.0 / 3.0,
			1.0 / 3.0,
		}}, 1
	}

	s01 := s0*c1 + s1*c0
	var outN int
	var out [4]struct {
		cbez CubicBez
		d0   float64
		d1   float64
	}
	for _, root := range roots {
		var d0, d1 float64
		if root > 0.0 {
			d := (root*s0 - area*(10.0/3.0)) / (0.5*root*s01 - s1)
			if d > 0.0 {
				d0, d1 = root, d
			} else {
				d0, d1 = s1/s01, 0.0
			}
		} else {
			d0, d1 = 0.0, s0/s01
		}
		// We could implement a maximum d value here.
		if d0 >= 0.0 && d1 >= 0.0 {
			out[outN] = struct {
				cbez CubicBez
				d0   float64
				d1   float64
			}{
				CubicBez{
					Pt(0.0, 0.0),
					Pt(d0*c0, d0*s0),
					Pt(1.0-d1*c1, d1*s1),
					Pt(1.0, 0.0),
				},
				d0,
				d1,
			}
			outN++
		}
	}
	return out, outN
}

// CurveFitSample describes a sample point of a curve for fitting.
type CurveFitSample struct {
	// A point on the curve at the sample location.
	Point Point
	// A vector tangent to the curve at the sample location.
	Tangent Vec2
}

// Intersect intersects a ray orthogonal to the tangent with the given cubic.
//
// Returns a vector of t values on the cubic.
func (cfs CurveFitSample) Intersect(c CubicBez) ([3]float64, int) {
	p1 := c.P1.Sub(c.P0).Mul(3)
	p2 := Vec2(c.P2).Mul(3).Sub(Vec2(c.P1).Mul(6)).Add(Vec2(c.P0).Mul(3))
	p3 := c.P3.Sub(c.P0).Sub(c.P2.Sub(c.P1).Mul(3))
	c0 := c.P0.Sub(cfs.Point).Dot(cfs.Tangent)
	c1 := p1.Dot(cfs.Tangent)
	c2 := p2.Dot(cfs.Tangent)
	c3 := p3.Dot(cfs.Tangent)

	roots, n := SolveCubic(c0, c1, c2, c3)
	var out [3]float64
	nn := 0
	for _, t := range roots[:n] {
		if t >= 0.0 && t <= 1.0 {
			out[nn] = t
			nn++
		}
	}
	return out, nn
}

// FitToBezPath generates a Bézier path that fits the source curve.
//
// This function recursively subdivides the curve in half by the parameter when
// the accuracy is not met. That gives a reasonably optimized result but not
// necessarily the minimum number of segments.
//
// In general, the resulting Bézier path should have a [Fréchet distance] less
// than the provided accuracy parameter. However, this is not a rigorous
// guarantee, as the error metric is computed approximately.
//
// This function is intended for use when the source curve is piecewise
// continuous, with the discontinuities reported by the BreakCusp method. In
// applications (such as stroke expansion) where this property may not hold, it
// is up to the client to detect and handle such cases. Even so, best effort is
// made to avoid infinite subdivision.
//
// When a higher degree of optimization is desired (at considerably more runtime
// cost), consider using [FitToBezPathOpt] instead.
//
// See [Fitting cubic Bézier curves] and [Parallel curves of cubic Béziers] for
// an explanation of the approach used.
//
// [Fréchet distance]: https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance
// [Fitting cubic Bézier curves]: https://raphlinus.github.io/curves/2021/03/11/bezier-fitting.html
// [Parallel curves of cubic Béziers]: https://raphlinus.github.io/curves/2022/09/09/parallel-beziers.html
func FitToBezPath(source FittableCurve, accuracy float64) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		fitToBezPathRec(source, 0.0, 1.0, accuracy, yield, true)
	}
}

func fitToBezPathRec(
	source FittableCurve,
	start float64,
	end float64,
	accuracy float64,
	yield func(PathElement) bool,
	needMove bool,
) (cont, needMoveRet bool) {
	// Discussion question: possibly should take endpoint samples, to avoid
	// duplication of that work.

	startPt := source.SamplePtTangent(start, 1.0).Point
	endPt := source.SamplePtTangent(end, -1.0).Point
	if startPt.DistanceSquared(endPt) <= accuracy*accuracy {
		if c, _, ok := tryFitLine(source, accuracy, start, end, startPt, endPt); ok {
			if needMove {
				if !yield(MoveTo(c.P0)) {
					return false, true
				}
			}
			return yield(CubicTo(c.P1, c.P2, c.P3)), false
		}
	}
	var t float64
	if cusp, ok := source.BreakCusp(start, end); ok {
		t = cusp
	} else if c, _, ok := FitToCubic(source, start, end, accuracy); ok {
		if needMove {
			if !yield(MoveTo(c.P0)) {
				return false, true
			}
		}
		return yield(CubicTo(c.P1, c.P2, c.P3)), false
	} else {
		// A smarter approach is possible than midpoint subdivision, but would be
		// a significant increase in complexity.
		t = 0.5 * (start + end)
	}
	if t == start || t == end {
		// infinite recursion, just draw a line
		p1 := startPt.Lerp(endPt, 1.0/3.0)
		p2 := endPt.Lerp(startPt, 1.0/3.0)
		if needMove {
			if !yield(MoveTo(startPt)) {
				return false, true
			}
		}
		return yield(CubicTo(p1, p2, endPt)), false
	}
	cont, needMove = fitToBezPathRec(source, start, t, accuracy, yield, needMove)
	if !cont {
		return false, needMove
	}
	cont, needMove = fitToBezPathRec(source, t, end, accuracy, yield, needMove)
	if !cont {
		return false, needMove
	}
	return true, needMove
}

// tryFitLine tries fitting a line.
//
// This is especially useful for very short chords, in which the standard
// cubic fit is not numerically stable. The tangents are not considered, so
// it's useful in cusp and near-cusp situations where the tangents are not
// reliable, as well.
//
// Returns the line raised to a cubic and the error, if within tolerance.
func tryFitLine(
	source FittableCurve,
	accuracy float64,
	rangeStart float64,
	rangeEnd float64,
	start Point,
	end Point,
) (CubicBez, float64, bool) {
	acc2 := accuracy * accuracy
	chordLine := Line{start, end}
	const shortN = 7
	maxErr2 := 0.0
	dt := (rangeEnd - rangeStart) / float64(shortN+1)
	for i := range shortN {
		t := rangeStart + float64(i+1)*dt
		p, _ := source.SamplePtDeriv(t)
		err2, _ := chordLine.Nearest(p, accuracy)
		if err2 > acc2 {
			// Not in tolerance; likely subdivision will help.
			return CubicBez{}, 0, false
		}
		maxErr2 = max(err2, maxErr2)
	}
	p1 := start.Lerp(end, 1.0/3.0)
	p2 := end.Lerp(start, 1.0/3.0)
	c := CubicBez{start, p1, p2, end}
	return c, maxErr2, true
}

// FitToBezPathOpt generates a highly optimized Bézier path that fits the source
// curve.
//
// This function is considerably slower than [FitToBezPath], as it computes
// optimal subdivision points. Its result is expected to be very close to the
// optimum possible Bézier path for the source curve, in that it has a minimal
// number of curve segments, and a minimal error over all paths with that number
// of segments.
//
// In general, the resulting Bézier path should have a [Fréchet distance] less
// than the provided accuracy parameter. However, this is not a rigorous
// guarantee, as the error metric is computed approximately.
//
// [Fréchet distance]: https://en.wikipedia.org/wiki/Fr%C3%A9chet_distance
func FitToBezPathOpt(source FittableCurve, accuracy float64) BezPath {
	var cusps []float64
	path := BezPath{}
	t0 := 0.0
	for {
		t1 := 1.0
		if len(cusps) > 0 {
			t1 = cusps[len(cusps)-1]
		}

		if t, ok := fitToBezPathOptInner(source, accuracy, t0, t1, &path); ok {
			cusps = append(cusps, t)
		} else {
			if len(cusps) > 0 {
				t := cusps[len(cusps)-1]
				t0 = t
				cusps = cusps[:len(cusps)-1]
			} else {
				break
			}
		}
	}
	return path
}

// fitToBezPathOptInner fits a range without cusps.
//
// Returns true and a cusp's location if a cusp was found.
func fitToBezPathOptInner(
	source FittableCurve,
	accuracy float64,
	rangeStart float64,
	rangeEnd float64,
	path *BezPath,
) (float64, bool) {
	if t, ok := source.BreakCusp(rangeStart, rangeEnd); ok {
		return t, true
	}
	var err float64
	if c, err2, ok := FitToCubic(source, rangeStart, rangeEnd, accuracy); ok {
		err = math.Sqrt(err2)
		if err < accuracy {
			if rangeStart == 0.0 {
				path.MoveTo(c.P0)
			}
			path.CubicTo(c.P1, c.P2, c.P3)
			return 0, false
		}
	} else {
		err = 2.0 * accuracy
	}
	t0, t1 := rangeStart, rangeEnd
	n := uint(0)
	var lastErr float64
loop:
	for {
		n++
		switch kind, v := fitOptSegment(source, accuracy, t0, t1); kind {
		case FitResultParamVal:
			t0 = v
		case FitResultSegmentError:
			lastErr = v
			break loop
		case FitResultCuspFound:
			return v, true
		}
	}
	t0 = rangeStart
	const epsilon = 1e-9
	f := func(x float64) result[float64, float64] {
		return fitOptErrDelta(source, x, accuracy, t0, t1, n)
	}
	k1 := 0.2 / accuracy
	ya := -err
	yb := accuracy - lastErr
	var x float64
	if k := solveITPFallible(f, 0.0, accuracy, epsilon, 1, k1, ya, yb); k.isOK {
		x = k.ok[1]
	} else {
		return k.err, true
	}
	pathLen := len(*path)
	for i := range n {
		var t1 float64
		if i < n-1 {
			switch kind, v := fitOptSegment(source, x, t0, rangeEnd); kind {
			case FitResultParamVal:
				t1 = v
			case FitResultSegmentError:
				t1 = rangeEnd
			case FitResultCuspFound:
				path.Truncate(pathLen)
				return v, true
			}
		} else {
			t1 = rangeEnd
		}
		c, _, ok := FitToCubic(source, t0, t1, accuracy)
		if !ok {
			panic("unreachable")
		}
		if i == 0 && rangeStart == 0.0 {
			path.MoveTo(c.P0)
		}
		path.CubicTo(c.P1, c.P2, c.P3)
		t0 = t1
		if t0 == rangeEnd {
			// This is unlikely but could happen when not monotonic.
			break
		}
	}
	return 0, false
}

func measureOneSeg(source FittableCurve, rangeStart, rangeEnd float64, limit float64) (float64, bool) {
	_, err2, ok := FitToCubic(source, rangeStart, rangeEnd, limit)
	return math.Sqrt(err2), ok
}

type FitResultKind int

const (
	// The parameter value that meets the desired accuracy.
	FitResultParamVal FitResultKind = iota + 1
	// Error of the measured segment.
	FitResultSegmentError
	// The parameter value where a cusp was found.
	FitResultCuspFound
)

func fitOptSegment(source FittableCurve, accuracy float64, rangeStart, rangeEnd float64) (FitResultKind, float64) {
	if t, ok := source.BreakCusp(rangeStart, rangeEnd); ok {
		return FitResultCuspFound, t
	}
	missingErr := accuracy * 2.0
	var err float64
	if v, ok := measureOneSeg(source, rangeStart, rangeEnd, accuracy); ok {
		err = v
	} else {
		err = missingErr
	}
	if err <= accuracy {
		return FitResultSegmentError, err
	}
	t0, t1 := rangeStart, rangeEnd
	f := func(x float64) result[float64, float64] {
		if t, ok := source.BreakCusp(rangeStart, rangeEnd); ok {
			return result[float64, float64]{err: t}
		}
		var err float64
		f, ok := measureOneSeg(source, t0, x, accuracy)
		if ok {
			err = f
		} else {
			err = missingErr
		}
		return result[float64, float64]{isOK: true, ok: err - accuracy}
	}
	const epsilon = 1e-9
	k1 := 2.0 / (t1 - t0)
	if k := solveITPFallible(f, t0, t1, epsilon, 1, k1, -accuracy, err-accuracy); k.isOK {
		return FitResultParamVal, k.ok[0]
	} else {
		return FitResultCuspFound, k.err
	}
}

// fitOptErrDelta returns the delta error (accuracy - error of last segment) or
// a cusp.
func fitOptErrDelta(
	source FittableCurve,
	accuracy float64,
	limit float64,
	rangeStart, rangeEnd float64,
	n uint,
) result[float64, float64] {
	t0, t1 := rangeStart, rangeEnd
	for range n - 1 {
		switch kind, v := fitOptSegment(source, accuracy, t0, t1); kind {
		case FitResultParamVal:
			t0 = v
			// In this case, n - 1 will work, which of course means the error is highly
			// non-monotonic. We should probably harvest that solution.
		case FitResultSegmentError:
			return result[float64, float64]{isOK: true, ok: accuracy}
		case FitResultCuspFound:
			return result[float64, float64]{err: v}
		default:
			panic(fmt.Sprintf("invalid kind; %v", kind))
		}
	}
	var err float64
	v, ok := measureOneSeg(source, t0, t1, limit)
	if ok {
		err = v
	} else {
		err = accuracy * 2
	}
	return result[float64, float64]{isOK: true, ok: accuracy - err}
}

// curveDist is an acceleration structure for estimating curve distance.
type curveDist struct {
	samples    [numSamples]CurveFitSample
	arcparams  []float64
	rangeStart float64
	rangeEnd   float64
	// A "spicy" curve is one with potentially extreme curvature variation,
	// so use arc length measurement for better accuracy.
	spicy bool
}

func curveDistFromCurve(source FittableCurve, rangeStart, rangeEnd float64) curveDist {
	step := (rangeEnd - rangeStart) * (1.0 / float64(numSamples+1))
	var lastTan option[Vec2]
	spicy := false
	const spicyThresh = 0.2
	var samples [numSamples]CurveFitSample
	for i := range numSamples + 2 {
		sample := source.SamplePtTangent(rangeStart+float64(i)*step, 1.0)
		if lastTan.isSet {
			cross := sample.Tangent.Cross(lastTan.value)
			dot := sample.Tangent.Dot(lastTan.value)
			if math.Abs(cross) > spicyThresh*math.Abs(dot) {
				spicy = true
			}
		}
		lastTan.set(sample.Tangent)
		if i > 0 && i < numSamples+1 {
			samples[i-1] = sample
		}
	}
	return curveDist{
		samples:    samples,
		arcparams:  nil,
		rangeStart: rangeStart,
		rangeEnd:   rangeEnd,
		spicy:      spicy,
	}
}

func (cd *curveDist) computeArcParams(source FittableCurve) {
	const nSubsamples = 10
	start, end := cd.rangeStart, cd.rangeEnd
	dt := (end - start) * (1.0 / float64((numSamples+1)*nSubsamples))
	arclen := 0.0
	for i := range numSamples + 1 {
		for j := range nSubsamples {
			t := start + dt*(float64(i*nSubsamples+j)+0.5)
			_, deriv := source.SamplePtDeriv(t)
			arclen += deriv.Hypot()
		}
		if i < numSamples {
			cd.arcparams = append(cd.arcparams, arclen)
		}
	}
	arclenInv := 1.0 / arclen
	for i := range cd.arcparams {
		cd.arcparams[i] *= arclenInv
	}
}

// evalArc evaluates distance based on arc length parametrization.
func (cd *curveDist) evalArc(c CubicBez, acc2 float64) (float64, bool) {
	// TODO: this could perhaps be tuned.
	const epsilon = 1e-9
	cArclen := c.Arclen(epsilon)
	maxErr2 := 0.0
	for i := range cd.samples {
		sample := cd.samples[i]
		s := cd.arcparams[i]
		t := SolveForArclen(c, cArclen*s, epsilon)
		err := sample.Point.DistanceSquared(c.Eval(t))
		maxErr2 = max(err, maxErr2)
		if maxErr2 > acc2 {
			return 0, false
		}
	}
	return maxErr2, true
}

// evalRay evaluates distance to a cubic approximation.
//
// If distance exceeds stated accuracy, can return false. Note that
// acc2 is the square of the target.
//
// Returns the square of the error, which is intended to be a good
// approximation of the Fréchet distance.
func (cd *curveDist) evalRay(c CubicBez, acc2 float64) (float64, bool) {
	maxErr2 := 0.0
	for _, sample := range cd.samples {
		best := acc2 + 1.0
		roots, n := sample.Intersect(c)
		for _, t := range roots[:n] {
			err := sample.Point.DistanceSquared(c.Eval(t))
			best = min(best, err)
		}
		maxErr2 = max(best, maxErr2)
		if maxErr2 > acc2 {
			return 0, false
		}
	}
	return maxErr2, true
}

func (cd *curveDist) evalDist(source FittableCurve, c CubicBez, acc2 float64) (float64, bool) {
	// Always compute cheaper distance, hoping for early-out.
	rayDist, ok := cd.evalRay(c, acc2)
	if !ok {
		return 0, false
	}
	if !cd.spicy {
		return rayDist, true
	}
	if len(cd.arcparams) == 0 {
		cd.computeArcParams(source)
	}
	return cd.evalArc(c, acc2)
}
