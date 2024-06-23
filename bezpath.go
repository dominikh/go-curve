package curve

import (
	"fmt"
	"io"
	"iter"
	"math"
	"slices"
)

type PathElementKind int

const (
	/// Move directly to the point without drawing anything, starting a new
	/// subpath.
	MoveToKind PathElementKind = iota + 1
	/// Draw a line from the current location to the point.
	LineToKind
	/// Draw a quadratic bezier using the current location and the two points.
	QuadToKind
	/// Draw a cubic bezier using the current location and the three points.
	CubicToKind
	/// Close off the path.
	ClosePathKind
)

// / The element of a Bézier path.
// /
// / A valid path has `MoveTo` at the beginning of each subpath.
type PathElement struct {
	Kind PathElementKind
	P0   Point
	P1   Point
	P2   Point
}

func (el PathElement) String() string {
	var kind string
	switch el.Kind {
	case MoveToKind:
		kind = "MoveTo"
	case LineToKind:
		kind = "LineTo"
	case QuadToKind:
		kind = "QuadTo"
	case CubicToKind:
		kind = "CubicTo"
	case ClosePathKind:
		kind = "ClosePath"
	default:
		kind = "InvalidPathElement"
	}
	return fmt.Sprintf("%s(%s, %s, %s)", kind, el.P0, el.P1, el.P2)
}

func (el PathElement) Transform(aff Affine) PathElement {
	switch el.Kind {
	case MoveToKind:
		return MoveTo(el.P0.Transform(aff))
	case LineToKind:
		return LineTo(el.P0.Transform(aff))
	case QuadToKind:
		return QuadTo(el.P0.Transform(aff), el.P1.Transform(aff))
	case CubicToKind:
		return CubicTo(el.P0.Transform(aff), el.P1.Transform(aff), el.P2.Transform(aff))
	case ClosePathKind:
		return ClosePath()
	default:
		return PathElement{}
	}
}

func (el PathElement) IsInf() bool {
	return el.P0.IsInf() ||
		el.P1.IsInf() ||
		el.P2.IsInf()
}

func (el PathElement) IsNaN() bool {
	return el.P0.IsNaN() ||
		el.P1.IsNaN() ||
		el.P2.IsNaN()
}

func MoveTo(pt Point) PathElement {
	return PathElement{Kind: MoveToKind, P0: pt}
}

func LineTo(pt Point) PathElement {
	return PathElement{Kind: LineToKind, P0: pt}
}

func QuadTo(p0, p1 Point) PathElement {
	return PathElement{Kind: QuadToKind, P0: p0, P1: p1}
}

func CubicTo(p0, p1, p2 Point) PathElement {
	return PathElement{Kind: CubicToKind, P0: p0, P1: p1, P2: p2}
}

func ClosePath() PathElement {
	return PathElement{Kind: ClosePathKind}
}

type PathSegmentKind int

const (
	// A line segment.
	LineKind PathSegmentKind = iota + 1
	// A quadratic Bézier segment.
	QuadKind
	// A cubic Bézier segment.
	CubicKind
)

// PathSegment represents a segment of a Bézier path. This type acts as a sort of tagged
// union representing all possible path segments ([Path], [QuadBez], and [CubicBez]).
type PathSegment struct {
	// We don't use an interface for PathSegment because we want {Line, Quad,
	// Cubic}.Transform to return their respective types, not PathSegment. But we cannot
	// encode that in Go interfaces.
	//
	// This also avoids having to allocate for path segments.

	Kind PathSegmentKind
	P0   Point
	P1   Point
	P2   Point
	P3   Point
}

// BoundingBox implements Shape.
func (seg PathSegment) BoundingBox() Rect {
	switch seg.Kind {
	case LineKind:
		return seg.Line().BoundingBox()
	case QuadKind:
		return seg.Quad().BoundingBox()
	case CubicKind:
		return seg.Cubic().BoundingBox()
	default:
		return Rect{}
	}
}

// PathElements implements Shape.
func (seg PathSegment) PathElements(tolerance float64) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		if !yield(PathElement{Kind: MoveToKind, P0: seg.P0}) {
			return
		}
		switch seg.Kind {
		case LineKind:
			yield(PathElement{Kind: LineToKind, P0: seg.P1})
		case QuadKind:
			yield(PathElement{Kind: QuadToKind, P0: seg.P1, P1: seg.P2})
		case CubicKind:
			yield(PathElement{Kind: CubicToKind, P0: seg.P1, P1: seg.P2, P2: seg.P3})
		}
	}
}

// Perimeter implements Shape.
func (seg PathSegment) Perimeter(accuracy float64) float64 {
	return seg.Arclen(accuracy)
}

var _ Shape = PathSegment{}
var _ ParametricCurve = PathSegment{}

// Line returns the line represented by this segment. This is only valid when Kind ==
// LineKind.
func (seg PathSegment) Line() Line { return Line{seg.P0, seg.P1} }

// Quad returns the quadratic Bézier represented by this segment. This is only valid when Kind ==
// QuadKind.
func (seg PathSegment) Quad() QuadBez { return QuadBez{seg.P0, seg.P1, seg.P2} }

// Cubic converts seg to a cubic Bézier. This is valid for any Kind.
func (seg PathSegment) Cubic() CubicBez {
	switch seg.Kind {
	case LineKind:
		p0 := seg.P0
		p1 := seg.P1
		return CubicBez{p0, p0, p1, p1}
	case QuadKind:
		return seg.Quad().Raise()
	case CubicKind:
		return CubicBez{seg.P0, seg.P1, seg.P2, seg.P3}
	default:
		return CubicBez{}
	}
}

func (seg PathSegment) Transform(aff Affine) PathSegment {
	return PathSegment{
		Kind: seg.Kind,
		P0:   seg.P0.Transform(aff),
		P1:   seg.P0.Transform(aff),
		P2:   seg.P0.Transform(aff),
		P3:   seg.P0.Transform(aff),
	}
}

func (seg PathSegment) IsInf() bool {
	return seg.P0.IsInf() || seg.P1.IsInf() || seg.P2.IsInf() || seg.P3.IsInf()
}

func (seg PathSegment) IsNaN() bool {
	return seg.P0.IsNaN() || seg.P1.IsNaN() || seg.P2.IsNaN() || seg.P3.IsNaN()
}

func (seg PathSegment) Eval(t float64) Point {
	switch seg.Kind {
	case LineKind:
		return seg.Line().Eval(t)
	case QuadKind:
		return seg.Quad().Eval(t)
	case CubicKind:
		return seg.Cubic().Eval(t)
	default:
		return Point{}
	}
}

func (seg PathSegment) Subsegment(start, end float64) PathSegment {
	switch seg.Kind {
	case LineKind:
		return seg.Line().Subsegment(start, end).Seg()
	case QuadKind:
		return seg.Quad().Subsegment(start, end).Seg()
	case CubicKind:
		return seg.Cubic().Subsegment(start, end).Seg()
	default:
		return PathSegment{}
	}
}

func (seg PathSegment) Start() Point {
	return seg.Eval(0)
}

func (seg PathSegment) End() Point {
	return seg.Eval(1)
}

func (seg PathSegment) SubsegmentCurve(start, end float64) ParametricCurve {
	return seg.Subsegment(start, end)
}

func (seg PathSegment) Subdivide() (PathSegment, PathSegment) {
	return seg.Subsegment(0.0, 0.5), seg.Subsegment(0.5, 1.0)
}

func (seg PathSegment) SubdivideCurve() (ParametricCurve, ParametricCurve) {
	return seg.Subdivide()
}

func (seg PathSegment) Arclen(accuracy float64) float64 {
	switch seg.Kind {
	case LineKind:
		return seg.Line().Arclen(accuracy)
	case QuadKind:
		return seg.Quad().Arclen(accuracy)
	case CubicKind:
		return seg.Cubic().Arclen(accuracy)
	default:
		return 0
	}
}

func (seg PathSegment) SolveForArclen(arclen, accuracy float64) float64 {
	switch seg.Kind {
	case LineKind:
		return SolveForArclen(seg.Line(), arclen, accuracy)
	case QuadKind:
		return SolveForArclen(seg.Quad(), arclen, accuracy)
	case CubicKind:
		return SolveForArclen(seg.Cubic(), arclen, accuracy)
	default:
		return 0
	}
}

func (seg PathSegment) SignedArea() float64 {
	switch seg.Kind {
	case LineKind:
		return seg.Line().SignedArea()
	case QuadKind:
		return seg.Quad().SignedArea()
	case CubicKind:
		return seg.Cubic().SignedArea()
	default:
		return 0
	}
}

func (seg PathSegment) Nearest(pt Point, accuracy float64) (distSq, t float64) {
	switch seg.Kind {
	case LineKind:
		return seg.Line().Nearest(pt, accuracy)
	case QuadKind:
		return seg.Quad().Nearest(pt, accuracy)
	case CubicKind:
		return seg.Cubic().Nearest(pt, accuracy)
	default:
		return 0, 0
	}
}

func (seg PathSegment) Extrema() ([MaxExtrema]float64, int) {
	switch seg.Kind {
	case LineKind:
		return seg.Line().Extrema()
	case QuadKind:
		return seg.Quad().Extrema()
	case CubicKind:
		return seg.Cubic().Extrema()
	default:
		return [MaxExtrema]float64{}, 0
	}
}

// PathElement returns the PathElement corresponding to the segment, discarding the
// segment's starting point.
func (seg PathSegment) PathElement() PathElement {
	switch seg.Kind {
	case LineKind:
		return LineTo(seg.Line().P1)
	case QuadKind:
		return QuadTo(seg.Quad().P1, seg.Quad().P2)
	case CubicKind:
		return CubicTo(seg.Cubic().P1, seg.Cubic().P2, seg.Cubic().P3)
	default:
		return PathElement{}
	}
}

// Reverse returns a new PathSegment describing the same path as this one, but with the
// points reversed.
func (seg PathSegment) Reverse() PathSegment {
	switch seg.Kind {
	case LineKind:
		seg.P0, seg.P1 = seg.P1, seg.P0
		return seg
	case QuadKind:
		seg.P0, seg.P2 = seg.P2, seg.P0
		return seg
	case CubicKind:
		seg.P0, seg.P1, seg.P2, seg.P3 = seg.P3, seg.P2, seg.P1, seg.P0
		return seg
	default:
		return PathSegment{}
	}
}

// / Assumes split at extrema.
func (seg PathSegment) windingInner(pt Point) int {
	start := seg.Eval(0)
	end := seg.Eval(1)
	var sign int
	if end.Y > start.Y {
		if pt.Y < start.Y || pt.Y >= end.Y {
			return 0
		}
		sign = -1
	} else if end.Y < start.Y {
		if pt.Y < end.Y || pt.Y >= start.Y {
			return 0
		}
		sign = 1
	} else {
		return 0
	}
	switch seg.Kind {
	case LineKind:
		if pt.X < min(start.X, end.X) {
			return 0
		}
		if pt.X >= max(start.X, end.X) {
			return sign
		}
		// line equation ax + by = c
		a := end.Y - start.Y
		b := start.X - end.X
		c := a*start.X + b*start.Y
		if (a*pt.X+b*pt.Y-c)*float64(sign) <= 0.0 {
			return sign
		} else {
			return 0
		}
	case QuadKind:
		quad := seg.Quad()
		p1 := quad.P1
		if pt.X < min(start.X, end.X, p1.X) {
			return 0
		}
		if pt.X >= max(start.X, end.X, p1.X) {
			return sign
		}
		a := end.Y - 2.0*p1.Y + start.Y
		b := 2.0 * (p1.Y - start.Y)
		c := start.Y - pt.Y
		solution, n := SolveQuadratic(c, b, a)
		for _, t := range solution[:n] {
			if t >= 0.0 && t <= 1.0 {
				x := quad.Eval(t).X
				if pt.X >= x {
					return sign
				} else {
					return 0
				}
			}
		}
		return 0
	case CubicKind:
		cubic := seg.Cubic()
		p1 := cubic.P1
		p2 := cubic.P2
		if pt.X < min(start.X, end.X, p1.X, p2.X) {
			return 0
		}
		if pt.X >= max(start.X, end.X, p1.X, p2.X) {
			return sign
		}
		a := end.Y - 3.0*p2.Y + 3.0*p1.Y - start.Y
		b := 3.0 * (p2.Y - 2.0*p1.Y + start.Y)
		c := 3.0 * (p1.Y - start.Y)
		d := start.Y - pt.Y
		solution, n := SolveCubic(d, c, b, a)
		for _, t := range solution[:n] {
			if t >= 0.0 && t <= 1.0 {
				x := cubic.Eval(t).X
				if pt.X >= x {
					return sign
				} else {
					return 0
				}
			}
		}
		return 0
	default:
		return 0
	}
}

// / Compute the winding number contribution of a single segment.
// /
// / Cast a ray to the left and count intersections.
func (seg PathSegment) Winding(pt Point) int {
	exs, n := ExtremaRanges(seg)
	var w int
	for _, ex := range exs[:n] {
		w += seg.Subsegment(ex[0], ex[1]).windingInner(pt)
	}
	return w
}

// / An intersection of a [`Line`] and a [`PathSeg`].
// /
// / This can be generated with the [`PathSeg::intersect_line`] method.
type LineIntersection struct {
	/// The 'time' that the intersection occurs, on the line.
	///
	/// This value is in the range 0..1.
	LineT float64
	/// The 'time' that the intersection occurs, on the path segment.
	///
	/// This value is nominally in the range 0..1, although it may slightly exceed
	/// that range at the boundaries of segments.
	SegmentT float64
}

func (li LineIntersection) IsInf() bool {
	return math.IsInf(li.LineT, 0) || math.IsInf(li.SegmentT, 0)
}

func (li LineIntersection) IsNaN() bool {
	return math.IsNaN(li.LineT) || math.IsNaN(li.SegmentT)
}

// / Compute intersections against a line.
// /
// / Returns a vector of the intersections. For each intersection,
// / the `t` value of the segment and line are given.
// /
// / Note: This test is designed to be inclusive of points near the endpoints
// / of the segment. This is so that testing a line against multiple
// / contiguous segments of a path will be guaranteed to catch at least one
// / of them. In such cases, use higher level logic to coalesce the hits
// / (the `t` value may be slightly outside the range of 0..1).
// /
// / # Examples
// /
// / ```
// / # use kurbo::*;
// / let seg = PathSeg::Line(Line::new((0.0, 0.0), (2.0, 0.0)));
// / let line = Line::new((1.0, 2.0), (1.0, -2.0));
// / let intersection = seg.intersect_line(line);
// / assert_eq!(intersection.len(), 1);
// / let intersection = intersection[0];
// / assert_eq!(intersection.segment_t, 0.5);
// / assert_eq!(intersection.line_t, 0.5);
// /
// / let point = seg.eval(intersection.segment_t);
// / assert_eq!(point, Point::new(1.0, 0.0));
// / ```
func (seg PathSegment) IntersectLine(line Line) ([3]LineIntersection, int) {
	switch seg.Kind {
	case LineKind:
		return seg.Line().IntersectLine(line)
	case QuadKind:
		return seg.Quad().IntersectLine(line)
	case CubicKind:
		return seg.Cubic().IntersectLine(line)
	default:
		return [3]LineIntersection{}, 0
	}
}

// MinDistance encodes the minimum distance between two Bézier curves, as returned by
// [PathSegment.MinDist].
type MinDistance struct {
	// The shortest distance between any two points on the two curves.
	Distance float64
	// The position of the nearest point on the first curve, as a parameter.
	T1 float64
	// The position of the nearest point on the second curve, as a parameter.
	T2 float64
}

func (seg PathSegment) vec2s() ([4]Vec2, int) {
	switch seg.Kind {
	case LineKind:
		return [4]Vec2{
			Vec2(seg.P0),
			Vec2(seg.P1),
		}, 2
	case QuadKind:
		return [4]Vec2{
			Vec2(seg.P0),
			Vec2(seg.P1),
			Vec2(seg.P2),
		}, 3
	case CubicKind:
		return [4]Vec2{
			Vec2(seg.P0),
			Vec2(seg.P1),
			Vec2(seg.P2),
			Vec2(seg.P3),
		}, 4
	default:
		return [4]Vec2{}, 0
	}
}

// MinDist returns the minimum distance between two path segments.
func (seg PathSegment) MinDist(other PathSegment, accuracy float64) MinDistance {
	vecs1, n1 := seg.vec2s()
	vecs2, n2 := other.vec2s()
	ret := minDistParam(
		vecs1[:n1],
		vecs2[:n2],
		[2]float64{0.0, 1.0},
		[2]float64{0.0, 1.0},
		accuracy,
		math.Inf(1),
	)
	distance, t1, t2 := ret[0], ret[1], ret[2]
	return MinDistance{
		Distance: math.Sqrt(distance),
		T1:       t1,
		T2:       t2,
	}
}

// / Compute endpoint tangents of a path segment.
// /
// / This version is robust to the path segment not being a regular curve.
func (seg PathSegment) Tangents() (Vec2, Vec2) {
	switch seg.Kind {
	case LineKind:
		return seg.Line().Tangents()
	case QuadKind:
		return seg.Quad().Tangents()
	case CubicKind:
		return seg.Cubic().Tangents()
	default:
		panic(fmt.Sprintf("invalid PathSegment kind %v", seg.Kind))
	}
}

// / A Bézier path.
// /
// / These docs assume basic familiarity with Bézier curves; for an introduction,
// / see Pomax's wonderful [A Primer on Bézier Curves].
// /
// / This path can contain lines, quadratics ([`QuadBez`]) and cubics
// / ([`CubicBez`]), and may contain multiple subpaths.
// /
// / # Elements and Segments
// /
// / A Bézier path can be represented in terms of either 'elements' ([`PathEl`])
// / or 'segments' ([`PathSeg`]). Elements map closely to how Béziers are
// / generally used in PostScript-style drawing APIs; they can be thought of as
// / instructions for drawing the path. Segments more directly describe the
// / path itself, with each segment being an independent line or curve.
// /
// / These different representations are useful in different contexts.
// / For tasks like drawing, elements are a natural fit, but when doing
// / hit-testing or subdividing, we need to have access to the segments.
// /
// / Conceptually, a `BezPath` contains zero or more subpaths. Each subpath
// / *always* begins with a `MoveTo`, then has zero or more `LineTo`, `QuadTo`,
// / and `CurveTo` elements, and optionally ends with a `ClosePath`.
// /
// / ```
// / use kurbo::{BezPath, Rect, Shape, Vec2};
// / let accuracy = 0.1;
// / let rect = Rect::from_origin_size((0., 0.,), (10., 10.));
// / // these are equivalent
// / let path1 = rect.to_path(accuracy);
// / let path2: BezPath = rect.path_elements(accuracy).collect();
// /
// / // extend a path with another path:
// / let mut path = rect.to_path(accuracy);
// / let shifted_rect = rect + Vec2::new(5.0, 10.0);
// / path.extend(shifted_rect.to_path(accuracy));
// / ```
// /
// / You can iterate the elements of a `BezPath` with the [`iter`] method,
// / and the segments with the [`segments`] method:
// /
// / ```
// / use kurbo::{BezPath, Line, PathEl, PathSeg, Point, Rect, Shape};
// / let accuracy = 0.1;
// / let rect = Rect::from_origin_size((0., 0.,), (10., 10.));
// / // these are equivalent
// / let path = rect.to_path(accuracy);
// / let first_el = PathEl::MoveTo(Point::ZERO);
// / let first_seg = PathSeg::Line(Line::new((0., 0.), (10., 0.)));
// / assert_eq!(path.iter().next(), Some(first_el));
// / assert_eq!(path.segments().next(), Some(first_seg));
// / ```
// / In addition, if you have some other type that implements
// / `Iterator<Item=PathEl>`, you can adapt that to an iterator of segments with
// / the [`segments` free function].
// /
// / # Advanced functionality
// /
// / In addition to the basic API, there are several useful pieces of advanced
// / functionality available on `BezPath`:
// /
// / - [`flatten`] does Bézier flattening, converting a curve to a series of
// /   line segments
// / - [`intersect_line`] computes intersections of a path with a line, useful
// /   for things like subdividing
// /
// / [A Primer on Bézier Curves]: https://pomax.github.io/bezierinfo/
// / [`iter`]: BezPath::iter
// / [`segments`]: BezPath::segments
// / [`flatten`]: BezPath::flatten
// / [`intersect_line`]: PathSeg::intersect_line
// / [`segments` free function]: segments
// / [`FromIterator<PathEl>`]: std::iter::FromIterator
// / [`Extend<PathEl>`]: std::iter::Extend
type BezPath []PathElement

var _ Shape = BezPath{}

func (p BezPath) PathElements(tolerance float64) iter.Seq[PathElement] {
	return slices.Values([]PathElement(p))
}

// Transform returns a new path with an affine transformation to the path. See
// [BezPath.ApplyTransform] for a version that modifies the path in-place.
func (p BezPath) Transform(aff Affine) BezPath {
	els := make([]PathElement, len(p))
	for i := range p {
		els[i] = p[i].Transform(aff)
	}
	return els
}

// Pop removes and returns the last element in the path. If the path contains no more
// elements, false is returned.
func (p *BezPath) Pop() (PathElement, bool) {
	if len(*p) == 0 {
		return PathElement{}, false
	}
	el := (*p)[len(*p)-1]
	*p = (*p)[:len(*p)-1]
	return el, true
}

// Push adds an element to the path.
func (p *BezPath) Push(el PathElement) {
	*p = append(*p, el)
}

// Push pushes a "move to" element onto the path.
func (p *BezPath) MoveTo(pt Point) { p.Push(MoveTo(pt)) }

// / Push a "line to" element onto the path.
// /
// / Will panic with a debug assert when the path is empty and there is no
// / "move to" element on the path.
// /
// / If `line_to` is called immediately after `close_path` then the current
// / subpath starts at the initial point of the previous subpath.
func (p *BezPath) LineTo(pt Point) { p.Push(LineTo(pt)) }

// / Push a "quad to" element onto the path.
// /
// / Will panic with a debug assert when the path is empty and there is no
// / "move to" element on the path.
// /
// / If `quad_to` is called immediately after `close_path` then the current
// / subpath starts at the initial point of the previous subpath.
func (p *BezPath) QuadTo(p1, p2 Point) { p.Push(QuadTo(p1, p2)) }

// / Push a "curve to" element onto the path.
// /
// / Will panic with a debug assert when the path is empty and there is no
// / "move to" element on the path.
// /
// / If `curve_to` is called immediately after `close_path` then the current
// / subpath starts at the initial point of the previous subpath.
func (p *BezPath) CubicTo(p1, p2, p3 Point) { p.Push(CubicTo(p1, p2, p3)) }

// / Push a "close path" element onto the path.
// /
// / Will panic with a debug assert when the path is empty and there is no
// / "move to" element on the path.
func (p *BezPath) ClosePath() { p.Push(ClosePath()) }

// Segments returns an iterator over the path's segments.
func (p BezPath) Segments() iter.Seq[PathSegment] { return Segments(slices.Values(p)) }

// Elements returns an iterator over the path's elements.
func (p BezPath) Elements() iter.Seq[PathElement] { return slices.Values(p) }

// Truncate truncates the path, keeping the first n elements.
func (p *BezPath) Truncate(n int) {
	if n >= len(*p) {
		return
	}
	*p = (*p)[:n]
}

// Flatten flattens the path to a sequence of lines. See [Flatten] for details on the
// process.
func (p BezPath) Flatten(tolerance float64) iter.Seq[PathElement] {
	return Flatten(p.PathElements(0), tolerance)
}

// / Signed area.
func (p BezPath) SignedArea() float64 {
	return SegmentsSignedArea(p.Segments())
}

func (p BezPath) Perimeter(accuracy float64) float64 {
	return p.Arclen(accuracy)
}

func (p BezPath) Arclen(accuracy float64) float64 {
	return SegmentsPerimeter(p.Segments(), accuracy)
}

// / Winding number of point.
func (p BezPath) Winding(pt Point) int {
	return SegmentsWinding(p.Segments(), pt)
}

func (p BezPath) BoundingBox() Rect {
	return SegmentsBoundingBox(p.Segments())
}

// Segment returns the segment at the given element index, if any.
//
// If you need to access all segments, [BezPath.Segments] provides a better
// API. This function is intended for random access of specific elements, for clients
// that require this specifically.
//
// This returns the segment that ends at the provided element
// index. In effect this means it is 1-indexed: since no segment ends at
// the first element (which is presumed to be a [MoveTo]) Segment(0) will
// always return false.
func (p BezPath) Segment(idx int) (PathSegment, bool) {
	if idx == 0 || idx >= len(p) {
		return PathSegment{}, false
	}
	var last Point
	switch prev := p[idx-1]; prev.Kind {
	case MoveToKind:
		last = prev.P0
	case LineToKind:
		last = prev.P0
	case QuadToKind:
		last = prev.P1
	case CubicToKind:
		last = prev.P2
	default:
		return PathSegment{}, false
	}

	switch el := p[idx]; el.Kind {
	case LineToKind:
		return Line{last, el.P0}.Seg(), true
	case QuadToKind:
		return QuadBez{last, el.P0, el.P1}.Seg(), true
	case CubicToKind:
		return CubicBez{last, el.P0, el.P1, el.P2}.Seg(), true
	case ClosePathKind:
		for i := idx - 1; i >= 0; i-- {
			el := p[i]
			if el.Kind == MoveToKind && el.P0 != last {
				return Line{last, el.P0}.Seg(), true
			}
		}
		return PathSegment{}, false

	default:
		return PathSegment{}, false
	}
}

// HasSegments reports whether the path contains any segments. A path that consists only
// of MoveTo and ClosePath elements has no segments.
func (p BezPath) HasSegments() bool {
	for i := range p {
		el := p[i]
		if el.Kind != MoveToKind && el.Kind != ClosePathKind {
			return true
		}
	}
	return false
}

// ApplyTransform destructively applies an affine transformation to the path. See
// [BezPath.Transform] for a version that returns a new path instead.
func (p *BezPath) ApplyTransform(aff Affine) {
	for i := range *p {
		(*p)[i] = (*p)[i].Transform(aff)
	}
}

func (p BezPath) IsInf() bool {
	for i := range p {
		if p[i].IsInf() {
			return true
		}
	}
	return false
}

func (p BezPath) IsNaN() bool {
	for i := range p {
		if p[i].IsNaN() {
			return true
		}
	}
	return false
}

// ControlBox returns a rectangle that conservatively encloses the path.
//
// Unlike [BezPath.BoundingBox], this uses control points directly rather than computing
// tight bounds for curve elements.
func (p BezPath) ControlBox() Rect {
	first := true
	var cbox Rect
	addPt := func(pt Point) {
		if first {
			first = false
			cbox = NewRectFromPoints(pt, pt)
		} else {
			cbox = cbox.UnionPoint(pt)
		}
	}
	for i := range p {
		el := p[i]
		switch el.Kind {
		case MoveToKind, LineToKind:
			addPt(el.P0)
		case QuadToKind:
			addPt(el.P0)
			addPt(el.P1)
		case CubicToKind:
			addPt(el.P0)
			addPt(el.P1)
			addPt(el.P2)
		case ClosePathKind:
		}
	}

	return cbox
}

// / Convert the path to an SVG path string representation.
// /
// / The current implementation doesn't take any special care to produce a
// / short string (reducing precision, using relative movement).
func (p BezPath) SVG(opts SVGOptions) string {
	return SVG(p.Elements(), opts)
}

func (p BezPath) WriteSVG(w io.Writer, opts SVGOptions) error {
	return WriteSVG(w, p.Elements(), opts)
}

// / Returns a new path with the winding direction of all subpaths reversed.
func (p BezPath) ReverseSubpaths() BezPath {
	elements := p
	startIdx := 1
	startPt := Point{}
	reversed := BezPath(make([]PathElement, 0, len(elements)))
	// Pending move is used to capture degenerate subpaths that should
	// remain in the reversed output.
	pendingMove := false
	for ix, el := range elements {
		switch el.Kind {
		case MoveToKind:
			pt := el.P0
			if pendingMove {
				reversed.Push(MoveTo(startPt))
			}
			if startIdx < ix {
				reverseSubpath(startPt, elements[startIdx:ix], &reversed)
			}
			pendingMove = true
			startPt = pt
			startIdx = ix + 1
		case ClosePathKind:
			if startIdx <= ix {
				reverseSubpath(startPt, elements[startIdx:ix], &reversed)
			}
			reversed.Push(ClosePath())
			startIdx = ix + 1
			pendingMove = false
		default:
			pendingMove = false
		}
	}
	if startIdx < len(elements) {
		reverseSubpath(startPt, elements[startIdx:], &reversed)
	} else if pendingMove {
		reversed.Push(MoveTo(startPt))
	}
	return reversed
}

// EndPoint returns the end point of the path element, or false if none exists. It exists
// for all kinds except for [ClosePathKind].
func (el PathElement) EndPoint() (Point, bool) {
	switch el.Kind {
	case MoveToKind:
		return el.P0, true
	case LineToKind:
		return el.P0, true
	case QuadToKind:
		return el.P1, true
	case CubicToKind:
		return el.P2, true
	default:
		return Point{}, false
	}
}

// / Helper for reversing a subpath.
// /
// / The `els` parameter must not contain any `MoveTo` or `ClosePath` elements.
func reverseSubpath(startPt Point, els []PathElement, reversed *BezPath) {
	var endPt Point
	if len(els) > 0 {
		endPt, _ = els[len(els)-1].EndPoint()
	} else {
		endPt = startPt
	}
	reversed.Push(MoveTo(endPt))
	for ix := len(els) - 1; ix >= 0; ix-- {
		el := &els[ix]

		var endPt Point
		if ix > 0 {
			endPt, _ = els[ix-1].EndPoint()
		} else {
			endPt = startPt
		}
		switch el.Kind {
		case LineToKind:
			reversed.Push(LineTo(endPt))
		case QuadToKind:
			reversed.Push(QuadTo(el.P0, endPt))
		case CubicToKind:
			reversed.Push(CubicTo(el.P1, el.P0, endPt))
		default:
			panic("reverseSubpath expects MoveTo and ClosePath to be removed")
		}
	}
}

// Flatten flattens a sequence of arbitrary path elements to a sequence of lines that
// approximate the original curve.
//
// The tolerance value controls the maximum distance between the curved input
// segments and their polyline approximations. (In technical terms, this is the
// [Hausdorff distance]). The algorithm attempts to bound this distance
// by tolerance but this is not absolutely guaranteed. The appropriate value
// depends on the use, but for antialiased rendering, a value of 0.25 has been
// determined to give good results. The number of segments tends to scale as the
// inverse square root of tolerance.
//
// This algorithm is based on the blog post [Flattening quadratic Béziers],
// but with some refinements. For one, there is a more careful approximation at cusps. For
// two, the algorithm is extended to work with cubic Béziers as well, by first subdividing
// into quadratics and then computing the subdivision of each quadratic. However, as a
// clever trick, these quadratics are subdivided fractionally, and their endpoints are not
// included.
//
// [Flattening quadratic Béziers]: https://raphlinus.github.io/graphics/curves/2019/12/23/flatten-quadbez.html
// [Hausdorff distance]: https://en.wikipedia.org/wiki/Hausdorff_distance
func Flatten(seq iter.Seq[PathElement], tolerance float64) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		// / Proportion of tolerance budget that goes to cubic to quadratic conversion.
		const toQuadTol = 0.1

		sqrtTol := math.Sqrt(tolerance)
		var lastPt option[Point]
		quadBuf := []struct {
			q      QuadBez
			params flattenParams
		}{}
		for el := range seq {
			switch el.Kind {
			case MoveToKind:
				lastPt.set(el.P0)
				if !yield(el) {
					return
				}
			case LineToKind:
				lastPt.set(el.P0)
				if !yield(el) {
					return
				}
			case QuadToKind:
				p1, p2 := el.P0, el.P1
				if p0 := lastPt.value; lastPt.isSet {
					q := QuadBez{p0, p1, p2}
					params := q.estimateSubdiv(sqrtTol)
					n := max(int(math.Ceil(0.5*params.val/sqrtTol)), 1)
					step := 1.0 / float64(n)
					for i := 1; i < n; i++ {
						u := float64(i) * step
						t := q.determineSubdivT(&params, u)
						p := q.Eval(t)
						if !yield(LineTo(p)) {
							return
						}
					}
					if !yield(LineTo(p2)) {
						return
					}
				}
				lastPt.set(p2)
			case CubicToKind:
				p1, p2, p3 := el.P0, el.P1, el.P2
				if p0 := lastPt.value; lastPt.isSet {
					c := CubicBez{p0, p1, p2, p3}

					// Subdivide into quadratics, and estimate the number of
					// subdivisions required for each, summing to arrive at an
					// estimate for the number of subdivisions for the cubic.
					// Also retain these parameters for later.
					quadBuf = quadBuf[:0]
					sqrtRemainTol := sqrtTol * math.Sqrt(1.0-toQuadTol)
					sum := 0.0
					for quad := range c.Quadratics(tolerance * toQuadTol) {
						q := quad.Segment
						params := q.estimateSubdiv(sqrtRemainTol)
						sum += params.val
						quadBuf = append(quadBuf, struct {
							q      QuadBez
							params flattenParams
						}{q, params})
					}
					n := max(int(math.Ceil(0.5*sum/sqrtRemainTol)), 1)

					// Iterate through the quadratics, outputting the points of
					// subdivisions that fall within that quadratic.
					step := sum / float64(n)
					i := 1
					valSum := 0.0
					for _, thingy := range quadBuf {
						q := thingy.q
						params := thingy.params
						target := float64(i) * step
						recipVal := 1.0 / params.val
						for target < valSum+params.val {
							u := (target - valSum) * recipVal
							t := q.determineSubdivT(&params, u)
							p := q.Eval(t)
							if !yield(LineTo(p)) {
								return
							}
							i += 1
							if i == n+1 {
								break
							}
							target = float64(i) * step
						}
						valSum += params.val
					}
					if !yield(LineTo(p3)) {
						return
					}
				}
				lastPt.set(p3)
			case ClosePathKind:
				lastPt.clear()
				if !yield(el) {
					return
				}
			}
		}
	}
}

func SegmentsPerimeter(seq iter.Seq[PathSegment], accuracy float64) float64 {
	var sum float64
	for s := range seq {
		sum += s.Arclen(accuracy)
	}
	return sum
}

func SegmentsSignedArea(seq iter.Seq[PathSegment]) float64 {
	var sum float64
	for s := range seq {
		sum += s.SignedArea()
	}
	return sum
}

func SegmentsBoundingBox(seq iter.Seq[PathSegment]) Rect {
	var bbox Rect
	first := true
	for s := range seq {
		sbbox := BoundingBox(s)
		if first {
			first = false
			bbox = sbbox
		} else {
			bbox = bbox.Union(sbbox)
		}
	}
	return bbox
}

func SegmentsWinding(seq iter.Seq[PathSegment], pt Point) int {
	var sum int
	for s := range seq {
		sum += s.Winding(pt)
	}
	return sum
}
