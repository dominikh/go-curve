package curve

import (
	"iter"
	"math"
	"slices"
)

type Rect struct {
	X0, Y0 float64
	X1, Y1 float64
}

var _ ClosedShape = Rect{}

// NewRectFromPoints returns a rectangle with the extents of p0 and p1, ensuring that
// width and height are non-negative.
func NewRectFromPoints(p0, p1 Point) Rect {
	return Rect{p0.X, p0.Y, p1.X, p1.Y}.Abs()
}

// NewRectFromOrigin returns a rectangle with the given size, extending to the right and
// down (for positive sizes) from the origin. Width and height are ensured to be
// non-negative.
func NewRectFromOrigin(origin Point, size Size) Rect {
	return NewRectFromPoints(origin, origin.Translate(Vec2(size.AsVec2())))
}

// NewRectFromCenter returns a rectangle with the given size, centered around the center
// point.
func NewRectFromCenter(center Point, size Size) Rect {
	return Rect{
		X0: center.X - size.Width,
		Y0: center.Y - size.Height,
		X1: center.X + size.Width,
		Y1: center.Y + size.Height,
	}
}

// WithOrigin returns a new rectangle with the same size as r and a new origin.
func (r Rect) WithOrigin(origin Point) Rect {
	return NewRectFromOrigin(origin, r.Size())
}

// WithSize returns a new rectangle with the same origin as r and a new size.
func (r Rect) WithSize(size Size) Rect {
	return NewRectFromOrigin(r.Origin(), size)
}

// TODO(dh): Inset

// Abs returns a new rectangle with the same extents as r, but ensuring that width and
// height are non-negative.
func (r Rect) Abs() Rect {
	return Rect{
		X0: min(r.X0, r.X1),
		Y0: min(r.Y0, r.Y1),
		X1: max(r.X0, r.X1),
		Y1: max(r.Y0, r.Y1),
	}
}

func (r Rect) MinX() float64 { return min(r.X0, r.X1) }
func (r Rect) MaxX() float64 { return max(r.X0, r.X1) }
func (r Rect) MinY() float64 { return min(r.Y0, r.Y1) }
func (r Rect) MaxY() float64 { return max(r.Y0, r.Y1) }

// Origin returns the origin of the rectangle.
//
// This is the top left corner in a y-down space and with
// non-negative width and height.
func (r Rect) Origin() Point {
	return Point{
		X: r.X0,
		Y: r.Y0,
	}
}

// Width returns the rectangle's width, defined as X1 − X0. It may be negative.
func (r Rect) Width() float64 {
	return r.X1 - r.X0
}

// Height returns the rectangle's heigth, defined as Y1 − Y0. It may be negative.
func (r Rect) Height() float64 {
	return r.Y1 - r.Y0
}

func (r Rect) Size() Size {
	return Size{
		Width:  r.Width(),
		Height: r.Height(),
	}
}

func (r Rect) Center() Point {
	return Point{
		X: 0.5 * (r.X0 + r.X1),
		Y: 0.5 * (r.Y0 + r.Y1),
	}
}

func (r Rect) Contains(pt Point) bool {
	return pt.X >= r.X0 &&
		pt.X < r.X1 &&
		pt.Y >= r.Y0 &&
		pt.Y < r.Y1
}

// Union returns the smallest rectangle enclosing r and o.
//
// Results are valid only if width and height are non-negative.
func (r Rect) Union(o Rect) Rect {
	return Rect{
		X0: min(r.X0, o.X0),
		Y0: min(r.Y0, o.Y0),
		X1: max(r.X1, o.X1),
		Y1: max(r.Y1, o.Y1),
	}
}

// UnionPoint computes the union with one point.
//
// This method includes the perimeter of zero-area rectangles.
// Thus, a succession of UnionPoint operations on a series of
// points yields their enclosing rectangle.
//
// Results are valid only if width and height are non-negative.
func (r Rect) UnionPoint(pt Point) Rect {
	return Rect{
		X0: min(r.X0, pt.X),
		Y0: min(r.Y0, pt.Y),
		X1: max(r.X1, pt.X),
		Y1: max(r.Y1, pt.Y),
	}
}

// Intersect returns the intersection of two rectangles.
//
// The result is zero-area if either input has negative width or
// height. The result always has non-negative width and height.
func (r Rect) Intersect(o Rect) Rect {
	x0 := max(r.X0, o.X0)
	y0 := max(r.Y0, o.Y0)
	x1 := min(r.X1, o.X1)
	y1 := min(r.Y1, o.Y1)
	return Rect{
		X0: x0,
		Y0: y0,
		X1: max(x0, x1),
		Y1: max(y0, y1),
	}
}

// Inflate expands a rectangle by a constant amount in both directions.
//
// The logic simply applies the amount in each direction. If rectangle
// area or added dimensions are negative, this could give odd results.
func (r Rect) Inflate(width, height float64) Rect {
	return Rect{
		X0: r.X0 - width,
		Y0: r.Y0 - height,
		X1: r.X1 + width,
		Y1: r.Y1 + height,
	}
}

// Round returns a new rectangle with each coordinate value rounded to the nearest
// integer.
func (r Rect) Round() Rect {
	return Rect{
		X0: math.Round(r.X0),
		Y0: math.Round(r.Y0),
		X1: math.Round(r.X1),
		Y1: math.Round(r.Y1),
	}
}

// Ceil returns a new rectangle with each coordinate value rounded to the least integer
// value greater than or equal to it.
func (r Rect) Ceil() Rect {
	return Rect{
		X0: math.Ceil(r.X0),
		Y0: math.Ceil(r.Y0),
		X1: math.Ceil(r.X1),
		Y1: math.Ceil(r.Y1),
	}
}

// Floor returns a new rectangle, with each coordinate value rounded down to the
// nearest integer.
func (r Rect) Floor() Rect {
	return Rect{
		X0: math.Floor(r.X0),
		Y0: math.Floor(r.Y0),
		X1: math.Floor(r.X1),
		Y1: math.Floor(r.Y1),
	}
}

// Expand returns a new rectangle, with each coordinate value rounded away from
// the center of the rectangle to the nearest integer, unless they are already
// an integer. That is to say this function will return the smallest possible
// rectangle with integer coordinates that is a superset of the input
// rectangle.
func (r Rect) Expand() Rect {
	var x0, y0, x1, y1 float64
	if r.X0 < r.X1 {
		x0 = math.Floor(r.X0)
		x1 = math.Ceil(r.X1)
	} else {
		x0 = math.Ceil(r.X0)
		x1 = math.Floor(r.X1)
	}
	if r.Y0 < r.Y1 {
		y0 = math.Floor(r.Y0)
		y1 = math.Ceil(r.Y1)
	} else {
		y0 = math.Ceil(r.Y0)
		y1 = math.Floor(r.Y1)
	}
	return Rect{
		X0: x0,
		Y0: y0,
		X1: x1,
		Y1: y1,
	}
}

// Trunc returns a new rectangle, with each coordinate value rounded towards the
// center of the rectangle to the nearest integer, unless they are already an
// integer. That is to say this function will return the biggest possible
// rectangle with integer coordinates that is a subset of the input rectangle.
func (r Rect) Trunc() Rect {
	var x0, y0, x1, y1 float64
	if r.X0 < r.X1 {
		x0 = math.Ceil(r.X0)
		x1 = math.Floor(r.X1)
	} else {
		x0 = math.Floor(r.X0)
		x1 = math.Ceil(r.X1)
	}
	if r.Y0 < r.Y1 {
		y0 = math.Ceil(r.Y0)
		y1 = math.Floor(r.Y1)
	} else {
		y0 = math.Floor(r.Y0)
		y1 = math.Ceil(r.Y1)
	}
	return Rect{
		X0: x0,
		Y0: y0,
		X1: x1,
		Y1: y1,
	}
}

// ScaleFromOrigin scales the rectangle by the factor f with respect to the
// origin (the point (0, 0)).
func (r Rect) ScaleFromOrigin(f float64) Rect {
	return Rect{
		X0: r.X0 * f,
		Y0: r.Y0 * f,
		X1: r.X1 * f,
		Y1: r.Y1 * f,
	}
}

// AspectRatio returns the aspect ratio of the rectangle.
//
// This is defined as the height divided by the width. It measures the
// "squareness" of the rectangle (a value of 1 is square).
//
// If the width is 0 the output will be "sign(y1 - y0) * infinity".
//
// If The width and height are 0, the result will be NaN.
func (r Rect) AspectRatio() float64 {
	return r.Size().AspectRatio()
}

func (r Rect) IsInf() bool {
	return math.IsInf(r.X0, 0) ||
		math.IsInf(r.X1, 0) ||
		math.IsInf(r.Y0, 0) ||
		math.IsInf(r.Y1, 0)
}

func (r Rect) IsNaN() bool {
	return math.IsNaN(r.X0) ||
		math.IsNaN(r.X1) ||
		math.IsNaN(r.Y0) ||
		math.IsNaN(r.Y1)
}

func (r Rect) Translate(v Vec2) Rect {
	return Rect{
		X0: r.X0 + v.X,
		Y0: r.Y0 + v.Y,
		X1: r.X1 + v.X,
		Y1: r.Y1 + v.Y,
	}
}

func (r Rect) Area() float64 {
	return r.Width() * r.Height()
}

func (r Rect) BoundingBox() Rect {
	return r.Abs()
}

func (r Rect) Perimeter(accuracy float64) float64 {
	return 2 * (math.Abs(r.Width()) + math.Abs(r.Height()))
}

func (r Rect) Winding(pt Point) int {
	// Note: this function is carefully designed so that if the plane is
	// tiled with rectangles, the winding number will be nonzero for exactly
	// one of them.
	xmin := min(r.X0, r.X1)
	xmax := max(r.X0, r.X1)
	ymin := min(r.Y0, r.Y1)
	ymax := max(r.Y0, r.Y1)
	if pt.X >= xmin && pt.X < xmax && pt.Y >= ymin && pt.Y < ymax {
		if r.X1 > r.X0 != (r.Y1 > r.Y0) {
			return -1
		} else {
			return 1
		}
	} else {
		return 0
	}
}

func (r Rect) Path(tolerance float64) BezPath { return slices.Collect(r.PathElements(tolerance)) }

func (r Rect) PathElements(tolerance float64) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		_ = yield(MoveTo(Pt(r.X0, r.Y0))) &&
			yield(LineTo(Pt(r.X1, r.Y0))) &&
			yield(LineTo(Pt(r.X1, r.Y1))) &&
			yield(LineTo(Pt(r.X0, r.Y1))) &&
			yield(ClosePath())
	}
}

// RoundedRect creates a new [RoundedRect] from this rectangle and the provided
// corner radii.
func (r Rect) RoundedRect(radii RoundedRectRadii) RoundedRect {
	r = r.Abs()
	shortestSide := min(r.Width(), r.Height())
	radii = radii.Abs().Clamp(shortestSide / 2)
	return RoundedRect{
		Rect:  r,
		Radii: radii,
	}
}

// ContainedRectWithAspectRatio returns the largest possible rectangle that is
// fully contained in this rectangle, with the given aspect ratio.
//
// The aspect ratio is specified fractionally, as height / width.
//
// The resulting rectangle will be centered if it is smaller than the input
// rectangle.
func (r Rect) ContainedRectWithAspectRatio(aspectRatio float64) Rect {
	width, height := r.Width(), r.Height()
	rAspect := height / width

	// TODO the parameter 1e-9 was chosen quickly and may not be optimal.
	if math.Abs(rAspect-aspectRatio) < 1e-9 {
		return r
	} else if math.Abs(rAspect) < math.Abs(aspectRatio) {
		// shrink x to fit
		newWidth := height / aspectRatio
		gap := (width - newWidth) * 0.5
		x0 := r.X0 + gap
		x1 := r.X1 - gap
		return Rect{x0, r.Y0, x1, r.Y1}
	} else {
		// shrink y to fit
		newHeight := width * aspectRatio
		gap := (height - newHeight) * 0.5
		y0 := r.Y0 + gap
		y1 := r.Y1 - gap
		return Rect{r.X0, y0, r.X1, y1}
	}
}

type RoundedRect struct {
	Rect
	Radii RoundedRectRadii
}

var _ ClosedShape = RoundedRect{}

func NewRoundedRect(x0, y0, x1, y1, radius float64) RoundedRect {
	return RoundedRect{
		Rect{x0, y0, x1, y1},
		RoundedRectRadii{radius, radius, radius, radius},
	}
}

func (r RoundedRect) Area() float64 {
	// A corner is a quarter-circle, i.e.
	// .............#
	// .       ######
	// .    #########
	// .  ###########
	// . ############
	// .#############
	// ##############
	// |-----r------|
	// For each corner, we need to subtract the square that bounds this
	// quarter-circle, and add back in the area of quarter circle.

	corner := func(radius float64) float64 {
		return (math.Pi/4 - 1) * radius * radius
	}

	// Start with the area of the bounding rectangle. For each corner,
	// subtract the area of the corner under the quarter-circle, and add
	// back the area of the quarter-circle.
	return r.Rect.Area() +
		corner(r.Radii.TopLeft) +
		corner(r.Radii.TopRight) +
		corner(r.Radii.BottomRight) +
		corner(r.Radii.BottomLeft)
}

func dropFirst[T any](seq iter.Seq[T]) iter.Seq[T] {
	return func(yield func(T) bool) {
		first := true
		for el := range seq {
			if first {
				first = false
				continue
			}
			if !yield(el) {
				break
			}
		}
	}
}

func (r RoundedRect) PathElements(tolerance float64) iter.Seq[PathElement] {
	buildArcIter := func(i int, center Point, ellipseRadii Vec2) iter.Seq[PathElement] {
		a := Arc{
			Center:     center,
			Radii:      ellipseRadii,
			StartAngle: math.Pi / 2 * float64(i),
			SweepAngle: math.Pi / 2,
			XRotation:  0.0,
		}
		return dropFirst(a.PathElements(tolerance))
	}

	arcs := [...]iter.Seq[PathElement]{
		buildArcIter(
			2,
			Point{
				X: r.Rect.X0 + r.Radii.TopLeft,
				Y: r.Rect.Y0 + r.Radii.TopLeft,
			},
			Vec2{
				X: r.Radii.TopLeft,
				Y: r.Radii.TopLeft,
			},
		),
		buildArcIter(
			3,
			Point{
				X: r.Rect.X1 - r.Radii.TopRight,
				Y: r.Rect.Y0 + r.Radii.TopRight,
			},
			Vec2{
				X: r.Radii.TopRight,
				Y: r.Radii.TopRight,
			},
		),
		buildArcIter(
			0,
			Point{
				X: r.Rect.X1 - r.Radii.BottomRight,
				Y: r.Rect.Y1 - r.Radii.BottomRight,
			},
			Vec2{
				X: r.Radii.BottomRight,
				Y: r.Radii.BottomRight,
			},
		),
		buildArcIter(
			1,
			Point{
				X: r.Rect.X0 + r.Radii.BottomLeft,
				Y: r.Rect.Y1 - r.Radii.BottomLeft,
			},
			Vec2{
				X: r.Radii.BottomLeft,
				Y: r.Radii.BottomLeft,
			},
		),
	}

	rect := []PathElement{
		LineTo(Pt(
			r.Rect.X1-r.Radii.TopRight,
			r.Rect.Y0,
		)),
		LineTo(Pt(
			r.Rect.X1,
			r.Rect.Y1-r.Radii.BottomRight,
		)),
		LineTo(Pt(
			r.Rect.X0+r.Radii.BottomLeft,
			r.Rect.Y1,
		)),
		ClosePath(),
	}

	return func(yield func(PathElement) bool) {
		e := MoveTo(Pt(
			r.Rect.X0,
			r.Rect.Y0+r.Radii.TopLeft,
		))
		if !yield(e) {
			return
		}

		// Generate the arc curve elements.
		// When we've reached the end of the arc, add a line towards next arc.
		for i, arc := range arcs {
			for e := range arc {
				if !yield(e) {
					return
				}
			}
			e := rect[i]
			if !yield(e) {
				return
			}
		}
	}
}

func (r RoundedRect) Perimeter(accuracy float64) float64 {
	corner := func(radius float64) float64 {
		return (-2.0 + math.Pi/2) * radius * radius
	}

	// Start with the full perimeter. For each corner, subtract the
	// border surrounding the rounded corner and add the quarter-circle
	// perimeter.
	return r.Rect.Perimeter(1.0) +
		corner(r.Radii.TopLeft) +
		corner(r.Radii.TopRight) +
		corner(r.Radii.BottomRight) +
		corner(r.Radii.BottomLeft)
}

func (rr RoundedRect) Winding(pt Point) int {
	center := rr.Center()

	// 1. Translate the point relative to the center of the rectangle.
	pt.X -= center.X
	pt.Y -= center.Y

	// 2. Pick a radius value to use based on which quadrant the point is
	//    in.
	var radius float64
	switch {
	case pt.X < 0.0 && pt.Y < 0.0:
		radius = rr.Radii.TopLeft
	case pt.X >= 0.0 && pt.Y < 0.0:
		radius = rr.Radii.TopRight
	case pt.X >= 0.0 && pt.Y >= 0.0:
		radius = rr.Radii.BottomRight
	case pt.X < 0.0 && pt.Y >= 0.0:
		radius = rr.Radii.BottomLeft
	}

	// 3. This is the width and height of a rectangle with one corner at
	//    the center of the rounded rectangle, and another corner at the
	//    center of the relevant corner circle.
	insideHalfWidth := max(rr.Width()/2.0-radius, 0)
	insideHalfHeight := max(rr.Height()/2.0-radius, 0)

	// 4. Three things are happening here.
	//
	//    First, the x- and y-values are being reflected into the positive
	//    (bottom-right quadrant). The radius has already been determined,
	//    so it doesn't matter what quadrant is used.
	//
	//    After reflecting, the points are clamped so that their x- and y-
	//    values can't be lower than the x- and y- values of the center of
	//    the corner circle, and the coordinate system is transformed
	//    again, putting (0, 0) at the center of the corner circle.
	px := max(math.Abs(pt.X)-insideHalfWidth, 0.0)
	py := max(math.Abs(pt.Y)-insideHalfHeight, 0.0)

	// 5. The transforms above clamp all input points such that they will
	//    be inside the rounded rectangle if the corresponding output point
	//    (px, py) is inside a circle centered around the origin with the
	//    given radius.
	inside := px*px+py*py <= radius*radius
	if inside {
		return 1
	} else {
		return 0
	}
}

var _ Shape = RoundedRect{}

func (r RoundedRect) IsInf() bool {
	return r.Rect.IsInf() || r.Radii.IsInf()
}

func (r RoundedRect) IsNaN() bool {
	return r.Rect.IsNaN() || r.Radii.IsNaN()
}

type RoundedRectRadii struct {
	TopLeft     float64
	TopRight    float64
	BottomRight float64
	BottomLeft  float64
}

func (r RoundedRectRadii) Abs() RoundedRectRadii {
	return RoundedRectRadii{
		TopLeft:     math.Abs(r.TopLeft),
		TopRight:    math.Abs(r.TopRight),
		BottomLeft:  math.Abs(r.BottomLeft),
		BottomRight: math.Abs(r.BottomRight),
	}
}

func (r RoundedRectRadii) Clamp(max float64) RoundedRectRadii {
	return RoundedRectRadii{
		TopLeft:     min(r.TopLeft, max),
		TopRight:    min(r.TopRight, max),
		BottomLeft:  min(r.BottomLeft, max),
		BottomRight: min(r.BottomRight, max),
	}
}

func (r RoundedRectRadii) IsInf() bool {
	return math.IsInf(r.TopLeft, 0) ||
		math.IsInf(r.TopRight, 0) ||
		math.IsInf(r.BottomRight, 0) ||
		math.IsInf(r.BottomLeft, 0)
}

func (r RoundedRectRadii) IsNaN() bool {
	return math.IsNaN(r.TopLeft) ||
		math.IsNaN(r.TopRight) ||
		math.IsNaN(r.BottomRight) ||
		math.IsNaN(r.BottomLeft)
}
