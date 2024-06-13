package curve

import (
	"iter"
	"math"
)

// / Construct an affine transform from coefficients.
// /
// / If the coefficients are `(a, b, c, d, e, f)`, then the resulting
// / transformation represents this augmented matrix:
// /
// / ```text
// / | a c e |
// / | b d f |
// / | 0 0 1 |
// / ```
// /
// / Note that this convention is transposed from PostScript and
// / Direct2D, but is consistent with the
// / [Wikipedia](https://en.wikipedia.org/wiki/Affine_transformation)
// / formulation of affine transformation as augmented matrix. The
// / idea is that `(A * B) * v == A * (B * v)`, where `*` is the
// / [`Mul`](std::ops::Mul) trait.
type Affine struct {
	// We represent Affine as a struct instead of an array because Go applies fuck-all
	// optimizations to arrays, while structs benefit from SROA.

	N0, N1, N2, N3, N4, N5 float64
}

// / The identity transform.
var Identity = Affine{1, 0, 0, 1, 0, 0}

// / A transform that is flipped on the y-axis. Useful for converting between
// / y-up and y-down spaces.
var FlipY = Affine{1, 0, 0, -1, 0, 0}

// / A transform that is flipped on the x-axis.
var FlipX = Affine{-1, 0, 0, 1, 0, 0}

// / An affine transform representing non-uniform scaling
// / with different scale values for x and y
func Scale(x, y float64) Affine {
	return Affine{x, 0, 0, y, 0, 0}
}

// / An affine transform representing translation.
func Translate(v Vec2) Affine {
	return Affine{1, 0, 0, 1, v.X, v.Y}
}

// / An affine transform representing rotation.
// /
// / The convention for rotation is that a positive angle rotates a
// / positive X direction into positive Y. Thus, in a Y-down coordinate
// / system (as is common for graphics), it is a clockwise rotation, and
// / in Y-up (traditional for math), it is anti-clockwise.
// /
// / The angle, `th`, is expressed in radians.
func Rotate(th float64) Affine {
	sin, cos := math.Sincos(th)
	return Affine{cos, sin, -sin, cos, 0, 0}
}

// / An affine transform representing a rotation of `th` radians about `center`.
// /
// / See [`Affine::rotate()`] for more info.
func RotateAbout(th float64, center Point) Affine {
	c := Vec2(center)
	return Translate(c.Negate()).ThenRotate(th).ThenTranslate(c)
}

// / An affine transformation representing a skew.
// /
// / The `skew_x` and `skew_y` parameters represent skew factors for the
// / horizontal and vertical directions, respectively.
// /
// / This is commonly used to generate a faux oblique transform for
// / font rendering. In this case, you can slant the glyph 20 degrees
// / clockwise in the horizontal direction (assuming a Y-up coordinate
// / system):
// /
// / ```
// / let oblique_transform = kurbo::Affine::skew(20f64.to_radians().tan(), 0.0);
// / ```
func Skew(x, y float64) Affine {
	return Affine{1, y, x, 1, 0, 0}
}

// / Create an affine transform that represents reflection about the line `point + direction * t, t in (-infty, infty)`
// /
// / # Examples
// /
// / ```
// / # use kurbo::{Point, Vec2, Affine};
// / # fn assert_near(p0: Point, p1: Point) {
// / #     assert!((p1 - p0).hypot() < 1e-9, "{p0:?} != {p1:?}");
// / # }
// / let point = Point::new(1., 0.);
// / let vec = Vec2::new(1., 1.);
// / let map = Affine::reflect(point, vec);
// / assert_near(map * Point::new(1., 0.), Point::new(1., 0.));
// / assert_near(map * Point::new(2., 1.), Point::new(2., 1.));
// / assert_near(map * Point::new(2., 2.), Point::new(3., 1.));
// / ```
func Reflect(pt Point, direction Vec2) Affine {
	n := Vec2{
		X: direction.Y,
		Y: -direction.X,
	}.Normalize()

	// Compute Householder reflection matrix
	x2 := n.X * n.X
	xy := n.X * n.Y
	y2 := n.Y * n.Y
	// Here we also add in the post translation, because it doesn't require any further calc.
	aff := Affine{
		1.0 - 2.0*x2,
		-2.0 * xy,
		-2.0 * xy,
		1.0 - 2.0*y2,
		pt.X,
		pt.Y,
	}
	return aff.PreTranslate(Vec2(pt).Negate())
}

// Coefficients returns the the coefficients of the transform.
func (aff Affine) Coefficients() [6]float64 {
	return [6]float64{aff.N0, aff.N1, aff.N2, aff.N3, aff.N4, aff.N5}
}

// NewAffine creates a new affine transformation from an array of coefficients.
// Alternatively, you can initialize the fields of [Affine] manually.
func NewAffine(n [6]float64) Affine {
	return Affine{n[0], n[1], n[2], n[3], n[4], n[5]}
}

func (aff Affine) Mul(o Affine) Affine {
	return Affine{
		aff.N0*o.N0 + aff.N2*o.N1,
		aff.N1*o.N0 + aff.N3*o.N1,
		aff.N0*o.N2 + aff.N2*o.N3,
		aff.N1*o.N2 + aff.N3*o.N3,
		aff.N0*o.N4 + aff.N2*o.N5 + aff.N4,
		aff.N1*o.N4 + aff.N3*o.N5 + aff.N5,
	}
}

// / A rotation by th followed by aff.
// /
// / Equivalent to "aff * Rotate(th)"
func (aff Affine) PreRotate(th float64) Affine {
	return aff.Mul(Rotate(th))
}

// / `self` followed by a rotation of `th`.
// /
// / Equivalent to "Rotate(th) * aff"
func (aff Affine) ThenRotate(th float64) Affine {
	return Rotate(th).Mul(aff)
}

// / A rotation by `th` about `center` followed by `self`.
// /
// / Equivalent to "aff * RotateAbout(th)"
func (aff Affine) PreRotateAbout(th float64, center Point) Affine {
	return RotateAbout(th, center).Mul(aff)
}

// / `self` followed by a rotation of `th` about `center`.
// /
// / Equivalent to "RotateAbout(th, center) * aff"
func (aff Affine) ThenRotateAbout(th float64, center Point) Affine {
	return RotateAbout(th, center).Mul(aff)
}

// / A scale by `(scale_x, scale_y)` followed by `self`.
// /
// / Equivalent to `self * Scale(scale_x, scale_y)`
func (aff Affine) PreScale(x, y float64) Affine {
	return aff.Mul(Scale(x, y))
}

// / `self` followed by a scale of `(scale_x, scale_y)`.
// /
// / Equivalent to `Scale(scale_x, scale_y) * self`
func (aff Affine) ThenScale(x, y float64) Affine {
	return Scale(x, y).Mul(aff)
}

// / A translation of `trans` followed by `self`.
// /
// / Equivalent to `self * Translate(trans)`
func (aff Affine) PreTranslate(v Vec2) Affine {
	return aff.Mul(Translate(v))
}

// / `self` followed by a translation of `trans`.
// /
// / Equivalent to `Translate(trans) * self`
func (aff Affine) ThenTranslate(v Vec2) Affine {
	aff.N4 += v.X
	aff.N5 += v.Y
	return aff
}

// / Creates an affine transformation that takes the unit square to the given rectangle.
// /
// / Useful when you want to draw into the unit square but have your output fill any rectangle.
// / In this case push the `Affine` onto the transform stack.
func MapUnitSquare(rect Rect) Affine {
	return Affine{
		rect.Width(),
		0, 0,
		rect.Height(),
		rect.X0,
		rect.Y0,
	}
}

// / Compute the determinant of this transform.
func (aff Affine) Determinant() float64 {
	return aff.N0*aff.N3 - aff.N1*aff.N2
}

// / Compute the inverse transform.
// /
// / Produces NaN values when the determinant is zero.
func (aff Affine) Invert() Affine {
	invDet := 1 / aff.Determinant()
	return Affine{
		+invDet * aff.N3,
		-invDet * aff.N1,
		-invDet * aff.N2,
		+invDet * aff.N0,
		+invDet * (aff.N2*aff.N5 - aff.N3*aff.N4),
		+invDet * (aff.N1*aff.N4 - aff.N0*aff.N5),
	}
}

func (aff Affine) IsInf() bool {
	return math.IsInf(aff.N0, 0) ||
		math.IsInf(aff.N1, 0) ||
		math.IsInf(aff.N2, 0) ||
		math.IsInf(aff.N3, 0) ||
		math.IsInf(aff.N4, 0) ||
		math.IsInf(aff.N5, 0)
}

func (aff Affine) IsNaN() bool {
	return math.IsNaN(aff.N0) ||
		math.IsNaN(aff.N1) ||
		math.IsNaN(aff.N2) ||
		math.IsNaN(aff.N3) ||
		math.IsNaN(aff.N4) ||
		math.IsNaN(aff.N5)
}

// Compute the singular value decomposition of the linear transformation (ignoring the
// translation).
//
// All non-degenerate linear transformations can be represented as
//
//  1. a rotation about the origin.
//  2. a scaling along the x and y axes
//  3. another rotation about the origin
//
// composed together. Decomposing a 2x2 matrix in this way is called a "singular value
// decomposition" and is written `U Σ V^T`, where U and V^T are orthogonal (rotations) and Σ
// is a diagonal matrix (a scaling).
//
// Since currently this function is used to calculate ellipse radii and rotation from an
// affine map on the unit circle, we don't calculate V^T, since a rotation of the unit (or
// any) circle about its center always results in the same circle. This is the reason that an
// ellipse mapped using an affine map is always an ellipse.
//
// Will return NaNs if the matrix (or equivalently the linear map) is singular.
//
// First part of the return tuple is the scaling, second part is the angle of rotation (in
// radians)
func (aff Affine) svd() (scale Vec2, th float64) {
	a := aff.N0
	a2 := a * a
	b := aff.N1
	b2 := b * b
	c := aff.N2
	c2 := c * c
	d := aff.N3
	d2 := d * d
	ab := a * b
	cd := c * d
	th = math.Atan2(0.5*(2.0*(ab+cd)), a2-b2+c2-d2)
	s1 := a2 + b2 + c2 + d2
	s2 := math.Sqrt(math.Pow(a2-b2+c2-d2, 2) + 4.0*math.Pow(ab+cd, 2))
	return Vec2{
		X: math.Sqrt(0.5 * (s1 + s2)),
		Y: math.Sqrt(0.5 * (s1 - s2)),
	}, th
}

// / Returns the translation part of this affine map (`(self.0[4], self.0[5])`).
func (aff Affine) Translation() Vec2 {
	return Vec2{
		X: aff.N4,
		Y: aff.N5,
	}
}

// / Replaces the translation portion of this affine map
// /
// / The translation can be seen as being applied after the linear part of the map.
func (aff Affine) WithTranslation(v Vec2) Affine {
	aff.N4 = v.X
	aff.N5 = v.Y
	return aff
}

// / Compute the bounding box of a transformed rectangle.
// /
// / Returns the minimal `Rect` that encloses the given `Rect` after affine transformation.
// / If the transform is axis-aligned, then this bounding box is "tight", in other words the
// / returned `Rect` is the transformed rectangle.
// /
// / The returned rectangle always has non-negative width and height.
func (aff Affine) TransformRectBoundingBox(rect Rect) Rect {
	p00 := Pt(rect.X0, rect.Y0).Transform(aff)
	p01 := Pt(rect.X0, rect.Y1).Transform(aff)
	p10 := Pt(rect.X1, rect.Y0).Transform(aff)
	p11 := Pt(rect.X1, rect.Y1).Transform(aff)
	return NewRectFromPoints(p00, p01).Union(NewRectFromPoints(p10, p11))
}

func Transform[T interface{ Transform(Affine) T }](seq iter.Seq[T], aff Affine) iter.Seq[T] {
	return func(yield func(T) bool) {
		for v := range seq {
			if !yield(v.Transform(aff)) {
				break
			}
		}
	}
}
