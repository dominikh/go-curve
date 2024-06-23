package curve

import (
	"iter"
	"math"
)

//! Simplification of a Bézier path.
//!
//! This module is currently experimental.
//!
//! The methods in this module create a `SimplifyBezPath` object, which can then
//! be fed to [`fit_to_bezpath`] or [`fit_to_bezpath_opt`] depending on the degree
//! of optimization desired.
//!
//! The implementation uses a number of techniques to achieve high performance and
//! accuracy. The parameter (generally written `t`) evenly divides the curve segments
//! in the original, so sampling can be done in constant time. The derivatives are
//! computed analytically, as that is straightforward with Béziers.
//!
//! The areas and moments are computed analytically (using Green's theorem), and
//! the prefix sum is stored. Thus, it is possible to analytically compute the area
//! and moment of any subdivision of the curve, also in constant time, by taking
//! the difference of two stored prefix sum values, then fixing up the subsegments.
//!
//! A current limitation (hoped to be addressed in the future) is that non-regular
//! cubic segments may have tangents computed incorrectly. This can easily happen,
//! for example when setting a control point equal to an endpoint.
//!
//! In addition, this method does not report corners (adjoining segments where the
//! tangents are not continuous). It is not clear whether it's best to handle such
//! cases here, or in client code.
//!
//! [`fit_to_bezpath`]: crate::fit_to_bezpath
//! [`fit_to_bezpath_opt`]: crate::fit_to_bezpath_opt

type simplifyBezPath struct {
	ss []simplifyCubic
}

// / Set up a new Bézier path for simplification.
// /
// / Currently this is not dealing with discontinuities at all, but it
// / could be extended to do so.
func newSimplifyBezPath(seq iter.Seq[PathSegment]) simplifyBezPath {
	var a, x, y float64
	var ss []simplifyCubic
	for seg := range seq {
		c := seg.Cubic()
		ai, xi, yi := momentIntegrals2(c)
		a += ai
		x += xi
		y += yi
		ss = append(ss, simplifyCubic{c, [3]float64{a, x, y}})
	}
	return simplifyBezPath{ss}
}

// BreakCusp implements ParamCurveFit.
func (simp simplifyBezPath) BreakCusp(start float64, end float64) (float64, bool) {
	return 0, false
}

// SamplePtDeriv implements ParamCurveFit.
func (simp simplifyBezPath) SamplePtDeriv(t float64) (Point, Vec2) {
	i, t0 := simp.scale(t)
	n := len(simp.ss)
	if i == n {
		i -= 1
		t0 = 1.0
	}
	c := simp.ss[i].c
	return c.Eval(t0), Vec2(c.Differentiate().Eval(t0)).Mul(float64(n))
}

// SamplePtTangent implements ParamCurveFit.
func (simp simplifyBezPath) SamplePtTangent(t float64, sign float64) CurveFitSample {
	i, t0 := simp.scale(t)
	if i == len(simp.ss) {
		i -= 1
		t0 = 1.0
	}
	c := simp.ss[i].c
	p := c.Eval(t0)
	tangent := Vec2(c.Differentiate().Eval(t0))
	return CurveFitSample{p, tangent}
}

func (simp simplifyBezPath) MomentIntegrals(start, end float64) (float64, float64, float64) {
	// We could use the default implementation, but provide our own, mostly
	// because it is possible to efficiently provide an analytically accurate
	// answer.
	i0, t0 := simp.scale(start)
	i1, t1 := simp.scale(end)
	if i0 == i1 {
		return simp.momentIntegrals(i0, t0, t1)
	} else {
		a0, x0, y0 := simp.momentIntegrals(i0, t0, 1.0)
		a1, x1, y1 := simp.momentIntegrals(i1, 0.0, t1)
		a, x, y := a0+a1, x0+x1, y0+y1
		if i1 > i0+1 {
			axy2 := simp.ss[i0].moments
			axy3 := simp.ss[i1-1].moments
			a2, x2, y2 := axy2[0], axy2[1], axy2[2]
			a3, x3, y3 := axy3[0], axy3[1], axy3[2]
			a += a3 - a2
			x += x3 - x2
			y += y3 - y2
		}
		return a, x, y
	}
}

// / Resolve a `t` value to a cubic.
// /
// / Also return the resulting `t` value for the selected cubic.
func (simp simplifyBezPath) scale(t float64) (int, float64) {
	tScale := t * float64(len(simp.ss))
	tFloor := math.Floor(tScale)
	return int(tFloor), tScale - tFloor
}

func (simp simplifyBezPath) momentIntegrals(i int, start, end float64) (float64, float64, float64) {
	if end == start {
		return 0, 0, 0
	} else {
		return momentIntegrals2(simp.ss[i].c.Subsegment(start, end))
	}
}

var _ FittableCurve = simplifyBezPath{}

type SimplifyOptions struct {
	AngleThresh float64
	OptLevel    OptLevel
}

type simplifyCubic struct {
	c       CubicBez
	moments [3]float64
}

var DefaultSimplifyOptions = SimplifyOptions{1e-3, Subdivide}

// XXX find a better name
func momentIntegrals2(c CubicBez) (float64, float64, float64) {
	x0, y0 := c.P0.X, c.P0.Y
	x1, y1 := c.P1.X-x0, c.P1.Y-y0
	x2, y2 := c.P2.X-x0, c.P2.Y-y0
	x3, y3 := c.P3.X-x0, c.P3.Y-y0

	r0 := 3.0 * x1
	r1 := 3.0 * y1
	r2 := x2 * y3
	r3 := x3 * y2
	r4 := x3 * y3
	r5 := 27.0 * y1
	r6 := x1 * x2
	r7 := 27.0 * y2
	r8 := 45.0 * r2
	r9 := 18.0 * x3
	r10 := x1 * y1
	r11 := 30.0 * x1
	r12 := 45.0 * x3
	r13 := x2 * y1
	r14 := 45.0 * r3
	r15 := x1 * x1
	r16 := 18.0 * y3
	r17 := x2 * x2
	r18 := 45.0 * y3
	r19 := x3 * x3
	r20 := 30.0 * y1
	r21 := y2 * y2
	r22 := y3 * y3
	r23 := y1 * y1
	a := -r0*y2 - r0*y3 + r1*x2 + r1*x3 - 6.0*r2 + 6.0*r3 + 10.0*r4

	// Scale and add chord
	lift := x3 * y0
	area := a*0.05 + lift
	x := r10*r9 - r11*r4 + r12*r13 + r14*x2 - r15*r16 - r15*r7 - r17*r18 +
		r17*r5 +
		r19*r20 +
		105.0*r19*y2 +
		280.0*r19*y3 -
		105.0*r2*x3 +
		r5*r6 -
		r6*r7 -
		r8*x1
	y := -r10*r16 - r10*r7 - r11*r22 + r12*r21 + r13*r7 + r14*y1 - r18*x1*y2 +
		r20*r4 -
		27.0*r21*x1 -
		105.0*r22*x2 +
		140.0*r22*x3 +
		r23*r9 +
		27.0*r23*x2 +
		105.0*r3*y3 -
		r8*y2

	mx := x*(1.0/840.0) + x0*area + 0.5*x3*lift
	my := y*(1.0/420.0) + y0*a*0.1 + y0*lift

	return area, mx, my
}

type simplifyState struct {
	queue       BezPath
	needsMoveTo bool

	accuracy float64
	options  SimplifyOptions
	yield    func(PathElement) bool
}

func (ss *simplifyState) addSegment(seg PathSegment) {
	if !ss.queue.HasSegments() {
		ss.queue.MoveTo(seg.Start())
	}
	switch seg.Kind {
	case LineKind:
		ss.queue.LineTo(seg.P1)
	case QuadKind:
		ss.queue.QuadTo(seg.P1, seg.P2)
	case CubicKind:
		ss.queue.CubicTo(seg.P1, seg.P2, seg.P3)
	}
}

func (ss *simplifyState) flush() bool {
	if !ss.queue.HasSegments() {
		return true
	}

	if len(ss.queue) == 2 {
		// Queue is just one segment (count is moveto + primitive)
		// Just output the segment, no simplification is possible.
		els := ss.queue
		if !ss.needsMoveTo {
			els = els[1:]
		}
		for _, el := range els {
			if !ss.yield(el) {
				break
			}
		}
	} else {
		s := newSimplifyBezPath(Segments(ss.queue.PathElements(0)))
		var b BezPath
		switch ss.options.OptLevel {
		case Subdivide:
			first := true
			for el := range FitToBezPath(s, ss.accuracy) {
				if first {
					first = false
					if !ss.needsMoveTo {
						continue
					}
				}
				if !ss.yield(el) {
					return false
				}
			}
		case Optimized:
			b = FitToBezPathOpt(s, ss.accuracy)
			els := b
			if !ss.needsMoveTo {
				els = els[1:]
			}
			for _, el := range els {
				if !ss.yield(el) {
					return false
				}
			}
		}

	}
	ss.needsMoveTo = false
	ss.queue.Truncate(0)
	return true
}

// / Simplify a Bézier path.
// /
// / This function simplifies an arbitrary Bézier path; it is designed to handle
// / multiple subpaths and also corners.
// /
// / The underlying curve-fitting approach works best if the source path is very
// / smooth. If it contains higher frequency noise, then results may be poor, as
// / the resulting curve matches the original with G1 continuity at each subdivision
// / point, and also preserves the area. For such inputs, consider some form of
// / smoothing or low-pass filtering before simplification. In particular, if the
// / input is derived from a sequence of points, consider fitting a smooth spline.
// /
// / We may add such capabilities in the future, possibly as opt-in smoothing
// / specified through the options.
func Simplify(
	path iter.Seq[PathElement],
	accuracy float64,
	options SimplifyOptions,
) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		var lastPt option[Point]
		var lastSeg option[PathSegment]
		state := simplifyState{
			accuracy: accuracy,
			options:  options,
			yield:    yield,
		}

		for el := range path {
			var thisSeg option[PathSegment]
			switch el.Kind {
			case MoveToKind:
				state.flush()
				state.needsMoveTo = true
				lastPt.set(el.P0)
			case LineToKind:
				last := lastPt.unwrap()
				if last == el.P0 {
					continue
				}
				thisSeg.set(Line{last, el.P0}.Seg())
			case QuadToKind:
				last := lastPt.unwrap()
				if last == el.P0 && last == el.P1 {
					continue
				}
				thisSeg.set(QuadBez{last, el.P0, el.P1}.Seg())
			case CubicToKind:
				last := lastPt.unwrap()
				if last == el.P0 && last == el.P1 && last == el.P2 {
					continue
				}
				thisSeg.set(CubicBez{last, el.P0, el.P1, el.P2}.Seg())
			case ClosePathKind:
				state.flush()
				if !yield(ClosePath()) {
					return
				}
				state.needsMoveTo = true
				lastSeg.clear()
				continue
			}
			if seg := thisSeg.value; thisSeg.isSet {
				if last := lastSeg.value; lastSeg.isSet {
					_, lastTan := last.Tangents()
					thisTan, _ := seg.Tangents()
					if math.Abs(lastTan.Cross(thisTan)) > math.Abs(lastTan.Dot(thisTan))*options.AngleThresh {
						state.flush()
					}
				}
				lastPt.set(seg.End())
				state.addSegment(seg)
			}
			lastSeg = thisSeg
		}
		state.flush()
	}
}
