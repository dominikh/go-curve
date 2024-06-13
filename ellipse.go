package curve

import (
	"iter"
	"math"
)

type Ellipse struct {
	inner Affine
}

var _ ClosedShape = Ellipse{}

// / Create A new ellipse with a given center, radii, and rotation.
// /
// / The returned ellipse will be the result of taking a circle, stretching
// / it by the `radii` along the x and y axes, then rotating it from the
// / x axis by `rotation` radians, before finally translating the center
// / to `center`.
// /
// / Rotation is clockwise in a y-down coordinate system. For more on
// / rotation, see [`Affine::rotate`].
func NewEllipse(center Point, radii Vec2, xRotation float64) Ellipse {
	rx, ry := radii.Splat()
	return newEllipse(Vec2(center), rx, ry, xRotation)
}

// / Returns the largest ellipse that can be bounded by this [`Rect`].
// /
// / This uses the absolute width and height of the rectangle.
// /
// / This ellipse is always axis-aligned; to apply rotation you can call
// / [`with_rotation`] with the result.
// /
// / [`with_rotation`]: Ellipse::with_rotation
func NewEllipseFromRect(rect Rect) Ellipse {
	center := Vec2(rect.Center())
	width, height := rect.Size().Scale(1.0 / 2.0).Splat()
	return newEllipse(center, width, height, 0.0)
}

// NewEllipseFromAffine creates an ellipse from an affine transformation of the unit
// circle.
func NewEllipseFromAffine(aff Affine) Ellipse {
	return Ellipse{inner: aff}
}

func NewEllipseFromCircle(c Circle) Ellipse {
	return NewEllipse(c.Center, Vec(c.Radius, c.Radius), 0)
}

// / Create a new `Ellipse` centered on the provided point.
func (e Ellipse) WithCenter(center Point) Ellipse {
	return Ellipse{inner: e.inner.WithTranslation(Vec2(center))}
}

// / Create a new `Ellipse` with the provided radii.
func (e Ellipse) WithRadii(radii Vec2) Ellipse {
	_, rotation := e.inner.svd()
	translation := e.inner.Translation()
	return newEllipse(translation, radii.X, radii.Y, rotation)
}

// / Create a new `Ellipse`, with the rotation replaced by `rotation`
// / radians.
// /
// / The rotation is clockwise, for a y-down coordinate system. For more
// / on rotation, See [`Affine::rotate`].
func (e Ellipse) WithRotation(rotation float64) Ellipse {
	scale, _ := e.inner.svd()
	translation := e.inner.Translation()
	return newEllipse(translation, scale.X, scale.Y, rotation)
}

func newEllipse(center Vec2, scaleX, scaleY, xRotation float64) Ellipse {
	// Since the circle is symmetric about the x and y axes, using absolute values for the
	// radii results in the same ellipse. For simplicity we make this change here.
	return Ellipse{
		inner: Translate(center).
			Mul(Rotate(xRotation)).
			Mul(Scale(math.Abs(scaleX), math.Abs(scaleY))),
	}
}

// Contains implements ClosedShape.
func (e Ellipse) Contains(pt Point) bool {
	return e.Winding(pt) != 0
}

func (e Ellipse) IsInf() bool {
	return e.inner.IsInf()
}

func (e Ellipse) IsNaN() bool {
	return e.inner.IsNaN()
}

// Area implements ClosedShape.
func (e Ellipse) Area() float64 {
	x, y := e.Radii().Splat()
	return math.Pi * x * y
}

// BoundingBox implements ClosedShape.
func (e Ellipse) BoundingBox() Rect {
	// Compute a tight bounding box of the ellipse.
	//
	// See https://www.iquilezles.org/www/articles/ellipses/ellipses.htm. We can get the two
	// radius vectors by applying the affine map to the two impulses (1, 0) and (0, 1) which gives
	// (a, b) and (c, d) if the affine map is
	//
	//  a | c | e
	// -----------
	//  b | d | f
	//
	//  We can then use the method in the link with the translation to get the bounding box.

	aff := e.inner.Coefficients()
	a2 := aff[0] * aff[0]
	b2 := aff[1] * aff[1]
	c2 := aff[2] * aff[2]
	d2 := aff[3] * aff[3]
	cx := aff[4]
	cy := aff[5]
	rangeX := math.Sqrt(a2 + c2)
	rangeY := math.Sqrt(b2 + d2)
	return Rect{
		X0: cx - rangeX,
		Y0: cy - rangeY,
		X1: cx + rangeX,
		Y1: cy + rangeY,
	}
}

// PathElements implements ClosedShape.
func (e Ellipse) PathElements(tolerance float64) iter.Seq[PathElement] {
	radii, xRotation := e.inner.svd()
	return Arc{
		Center:     e.Center(),
		Radii:      radii,
		StartAngle: 0.0,
		SweepAngle: 2 * math.Pi,
		XRotation:  xRotation,
	}.PathElements(tolerance)
}

// Perimeter implements ClosedShape.
func (e Ellipse) Perimeter(accuracy float64) float64 {
	// TODO rather than delegate to the bezier path, it is possible to use various series
	// expansions to compute the perimeter to any accuracy. I believe Ramanujan authored the
	// quickest to converge. See
	// https://www.mathematica-journal.com/2009/11/23/on-the-perimeter-of-an-ellipse/
	// and https://en.wikipedia.org/wiki/Ellipse#Circumference
	//
	return SegmentsPerimeter(Segments(e.PathElements(0.1)), accuracy)
}

// Winding implements ClosedShape.
func (e Ellipse) Winding(pt Point) int {
	// Strategy here is to apply the inverse map to the point and see if it is in the unit
	// circle.
	inv := e.inner.Invert()
	if Vec2(pt.Transform(inv)).Hypot2() < 1.0 {
		return 1
	} else {
		return 0
	}
}

// Center returns the center of the ellipse.
func (e Ellipse) Center() Point {
	return Point(e.inner.Translation())
}

// Radii returns the two radii of the ellipse.
//
// The first number is the horizontal radius and the second is the
// vertical radius, before rotation.
func (e Ellipse) Radii() Vec2 {
	radii, _ := e.inner.svd()
	return radii
}

// / The ellipse's rotation, in radians.
// /
// / This allows all possible ellipses to be drawn by always starting with
// / an ellipse with the two radii on the x and y axes.
func (e Ellipse) Rotation() float64 {
	_, rot := e.inner.svd()
	return rot
}

// RadiiRotation returns the radii and the rotation of this ellipse.
//
// This is equivalent to, but more efficiant than, using [Ellipse.Radii] and
// [Ellipse.Rotation].
func (e Ellipse) RadiiRotation() (Vec2, float64) {
	return e.inner.svd()
}

func (e Ellipse) Translate(v Vec2) Ellipse {
	return Ellipse{
		inner: Translate(v).Mul(e.inner),
	}
}

func (e Ellipse) Transform(aff Affine) Ellipse {
	return Ellipse{
		inner: e.inner.Mul(aff),
	}
}
