package curve

import (
	"iter"
	"math"
)

// Line represents a line segment. It is both a [Shape] and a [ParametricCurve].
type Line struct {
	/// The line's start point.
	P0 Point
	/// The line's end point.
	P1 Point
}

var _ Shape = Line{}
var _ ParametricCurve = Line{}
var _ ArclenSolver = Line{}

// Length returns the length of the line.
func (l Line) Length() float64 {
	return l.P1.Sub(l.P0).Hypot()
}

// Arclen returns the length of the line
func (l Line) Arclen(accuracy float64) float64 {
	return l.Length()
}

func (l Line) SolveForArclen(arclen float64, accuracy float64) float64 {
	return arclen / l.P1.Sub(l.P0).Hypot()
}

// / Computes the point where two lines, if extended to infinity, would cross.
func (l Line) CrossingPoint(o Line) (Point, bool) {
	ab := l.P1.Sub(l.P0)
	cd := o.P1.Sub(o.P0)
	pcd := ab.Cross(cd)
	if pcd == 0 {
		return Point{}, false
	}
	h := ab.Cross(l.P0.Sub(o.P0)) / pcd
	return o.P0.Translate(cd.Mul(h)), true
}

func (l Line) IsInf() bool {
	return l.P0.IsInf() || l.P1.IsInf()
}

func (l Line) IsNaN() bool {
	return l.P0.IsNaN() || l.P1.IsNaN()
}

func (l Line) Translate(v Vec2) Line {
	return Line{
		P0: l.P0.Translate(v),
		P1: l.P1.Translate(v),
	}
}

func (l Line) BoundingBox() Rect {
	return Rect{
		X0: l.P0.X,
		Y0: l.P0.Y,
		X1: l.P1.X,
		Y1: l.P1.Y,
	}
}

func (l Line) Perimeter(accuracy float64) float64 {
	return l.Length()
}

func (l Line) PathElements(tolerance float64) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		_ = yield(MoveTo(l.P0)) &&
			yield(LineTo(l.P1))
	}
}

func (l Line) Eval(t float64) Point {
	return l.P0.Lerp(l.P1, t)
}

func (l Line) Nearest(pt Point, accuracy float64) (distSq, t float64) {
	d := l.P1.Sub(l.P0)
	dotp := d.Dot(pt.Sub(l.P0))
	dSquared := d.Dot(d)
	if dotp <= 0.0 {
		return pt.Sub(l.P0).Hypot2(), 0.0
	} else if dotp >= dSquared {
		return pt.Sub(l.P1).Hypot2(), 1.0
	} else {
		t := dotp / dSquared
		dist := pt.Sub(l.Eval(t)).Hypot2()
		return dist, t
	}
}

func (l Line) Transform(aff Affine) Line {
	return Line{
		P0: l.P0.Transform(aff),
		P1: l.P1.Transform(aff),
	}
}

func (l Line) Start() Point { return l.P0 }
func (l Line) End() Point   { return l.P1 }

func (l Line) Subsegment(start, end float64) Line {
	return Line{l.Eval(start), l.Eval(end)}
}

func (l Line) SubsegmentCurve(start, end float64) ParametricCurve {
	return l.Subsegment(start, end)
}

func (l Line) Subdivide() (Line, Line) {
	return l.Subsegment(0.0, 0.5), l.Subsegment(0.5, 1.0)
}

func (l Line) SubdivideCurve() (ParametricCurve, ParametricCurve) {
	return l.Subdivide()
}

func (l Line) Extrema() ([MaxExtrema]float64, int) {
	return [MaxExtrema]float64{}, 0
}

func (l Line) SignedArea() float64 {
	return Vec2(l.P0).Cross(Vec2(l.P1)) * 0.5
}

func (l Line) Tangents() (Vec2, Vec2) {
	d := l.P1.Sub(l.P0)
	return d, d
}

func (l Line) Seg() PathSegment {
	return PathSegment{Kind: LineKind, P0: l.P0, P1: l.P1}
}

func (l Line) IntersectLine(o Line) ([3]LineIntersection, int) {
	const epsilon = 1e-9
	p0 := o.P0
	p1 := o.P1
	dx := p1.X - p0.X
	dy := p1.Y - p0.Y

	det := dx*(l.P1.Y-l.P0.Y) - dy*(l.P1.X-l.P0.X)
	if math.Abs(det) < epsilon {
		// Lines are coincident (or nearly so).
		return [3]LineIntersection{}, 0
	}
	t := dx*(p0.Y-l.P0.Y) - dy*(p0.X-l.P0.X)
	// t = position on self
	t /= det
	if t >= -epsilon && t <= 1+epsilon {
		// u = position on probe line
		u :=
			(l.P0.X-p0.X)*(l.P1.Y-l.P0.Y) - (l.P0.Y-p0.Y)*(l.P1.X-l.P0.X)
		u /= det
		if u >= 0.0 && u <= 1.0 {
			return [3]LineIntersection{{u, t}}, 1
		}
	}
	return [3]LineIntersection{}, 0
}
