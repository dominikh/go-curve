package curve

import (
	"math"
)

// CubicOffset represents the [offset curve] of a cubic Bézier. It can be used
// for curve fitting, to produce a composite cubic Bézier representing a cubic
// Bézier's offset curve.
//
// See the [Parallel curves of cubic Béziers] blog post for a discussion of how
// this algorithm works and what kind of results can be expected. In general, it
// is expected to perform much better than most published algorithms. The number
// of curve segments needed to attain a given accuracy scales as O(n^6) with
// accuracy.
//
// [Parallel curves of cubic Béziers]: https://raphlinus.github.io/curves/2022/09/09/parallel-beziers.html
// [offset curve]: https://en.wikipedia.org/wiki/Parallel_curve
type CubicOffset struct {
	c  CubicBez
	q  QuadBez
	d  float64
	c0 float64
	c1 float64
	c2 float64
}

var _ FittableCurve = (*CubicOffset)(nil)

// / Create a new curve from Bézier segment and offset, with numerical robustness tweaks.
// /
// / The dimension represents a minimum feature size; the regularization is allowed to
// / perturb the curve by this amount in order to improve the robustness.
func NewCubicOffset(c CubicBez, d float64, dimension float64) CubicOffset {
	c = c.regularize(dimension)

	q := c.Differentiate()
	d0 := Vec2(q.P0)
	d1 := q.P1.Sub(q.P0).Mul(2)
	d2 := Vec2(q.P0).Sub(Vec2(q.P1).Mul(2)).Add(Vec2(q.P2))
	return CubicOffset{
		c:  c,
		q:  q,
		d:  d,
		c0: d * d1.Cross(d0),
		c1: d * 2.0 * d2.Cross(d0),
		c2: d * d2.Cross(d1),
	}
}

func (co *CubicOffset) evalOffset(t float64) Vec2 {
	dp := Vec2(co.q.Eval(t))
	norm := Vec(-dp.Y, dp.X)
	// TODO: deal with hypot = 0
	return norm.Mul(co.d).Mul(1.0 / dp.Hypot())
}

func (co *CubicOffset) Eval(t float64) Point {
	// Point on source curve.
	return co.c.Eval(t).Translate(co.evalOffset(t))
}

// / Evaluate derivative of curve.
func (co *CubicOffset) evalDeriv(t float64) Vec2 {
	return Vec2(co.q.Eval(t)).Mul(co.cuspSign(t))
}

// Compute a function which has a zero-crossing at cusps, and is positive at low
// curvatures on the source curve.
func (co *CubicOffset) cuspSign(t float64) float64 {
	ds2 := Vec2(co.q.Eval(t)).Hypot2()
	return ((co.c2*t+co.c1)*t+co.c0)/(ds2*math.Sqrt(ds2)) + 1.0
}

func (co *CubicOffset) SamplePtTangent(t float64, sign float64) CurveFitSample {
	p := co.Eval(t)
	const cuspEpsilon = 1e-8
	cusp := co.cuspSign(t)
	if math.Abs(cusp) < cuspEpsilon {
		// This is a numerical derivative, which is probably good enough for all
		// practical purposes, but an analytical derivative would be more
		// elegant.
		//
		// Also, we're not dealing with second or higher order cusps.
		cusp = sign * (co.cuspSign(t+cuspEpsilon) - co.cuspSign(t-cuspEpsilon))
	}
	tangent := Vec2(co.q.Eval(t))
	if math.Signbit(cusp) {
		tangent = tangent.Negate()
	}
	return CurveFitSample{p, tangent}
}

func (co *CubicOffset) SamplePtDeriv(t float64) (Point, Vec2) {
	return co.Eval(t), co.evalDeriv(t)
}

func (co *CubicOffset) BreakCusp(start, end float64) (float64, bool) {
	const cuspEpsilon = 1e-8
	// When an endpoint is on (or very near) a cusp, move just far enough
	// away from the cusp that we're confident we have the right sign.
	breakCuspHelper := func(x, d float64) (float64, float64) {
		cusp := co.cuspSign(x)
		for math.Abs(cusp) < cuspEpsilon && d < 1.0 {
			x += d
			oldCusp := cusp
			cusp = co.cuspSign(x)
			if math.Abs(cusp) > math.Abs(oldCusp) {
				break
			}
			d *= 2.0
		}
		return x, cusp
	}
	a, cusp0 := breakCuspHelper(start, 1e-12)
	b, cusp1 := breakCuspHelper(end, -1e-12)
	if a >= b || cusp0*cusp1 >= 0.0 {
		// Discussion point: maybe we should search for double cusps in the interior
		// of the range.
		return 0, false
	}
	s := sign(cusp1)
	f := func(t float64) float64 {
		return s * co.cuspSign(t)
	}
	k1 := 0.2 / (b - a)
	const itpEpsilon = 1e-12
	x := SolveITP(f, a, b, itpEpsilon, 1, k1, s*cusp0, s*cusp1)
	return x, true
}

func sign(x float64) float64 {
	if math.Signbit(x) {
		return -1
	} else {
		return 1
	}
}
