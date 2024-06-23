package curve

import (
	"fmt"
	"math"
)

type Point struct {
	X float64
	Y float64
}

// Pt returns the point (x, y).
func Pt(x, y float64) Point {
	return Point{X: x, Y: y}
}

func (pt Point) Splat() (float64, float64) {
	return pt.X, pt.Y
}

func (pt Point) String() string {
	return fmt.Sprintf("(%g, %g)", pt.X, pt.Y)
}

func (pt Point) Translate(o Vec2) Point {
	return Point{
		X: pt.X + o.X,
		Y: pt.Y + o.Y,
	}
}

func (pt Point) Transform(aff Affine) Point {
	return Point{
		X: aff.N0*pt.X + aff.N2*pt.Y + aff.N4,
		Y: aff.N1*pt.X + aff.N3*pt.Y + aff.N5,
	}
}

// Sub computes p−o.
// To subtract a vector from p, use Add and negate the vector.
func (pt Point) Sub(o Point) Vec2 {
	return Vec2{
		X: pt.X - o.X,
		Y: pt.Y - o.Y,
	}
}

// Lerp linearly interpolates between two points.
func (pt Point) Lerp(o Point, t float64) Point {
	return Point(Vec2(pt).Lerp(Vec2(o), t))
}

// Midpoint returns the midpoint of two points.
func (pt Point) Midpoint(o Point) Point {
	return Point{
		X: 0.5 * (pt.X + o.X),
		Y: 0.5 * (pt.Y + o.Y),
	}
}

// Distance returns the euclidean distance between two points.
func (pt Point) Distance(o Point) float64 {
	x := pt.X - o.X
	y := pt.Y - o.Y
	return math.Hypot(x, y)
}

// DistanceSquared returns the squared euclidean distance between two points.
func (pt Point) DistanceSquared(o Point) float64 {
	x := pt.X - o.X
	y := pt.Y - o.Y
	return x*x + y*y
}

// Round returns a new point with x and y rounded to the nearest integers.
func (pt Point) Round() Point {
	return Point{
		X: math.Round(pt.X),
		Y: math.Round(pt.Y),
	}
}

// Ceil returns a new point with x and y rounded up to the nearest integers.
func (pt Point) Ceil() Point {
	return Point{
		X: math.Ceil(pt.X),
		Y: math.Ceil(pt.Y),
	}
}

// Floor returns a new point with x and y rounded down to the nearest integers.
func (pt Point) Floor() Point {
	return Point{
		X: math.Floor(pt.X),
		Y: math.Floor(pt.Y),
	}
}

// Expand returns a new point with x and y rounded away from zero to the nearest integers.
func (pt Point) Expand() Point {
	return Point{
		X: expand(pt.X),
		Y: expand(pt.Y),
	}
}

// Trunc returns a new point with x and y rounded towards zero to the nearest integers.
func (pt Point) Trunc() Point {
	return Point{
		X: math.Trunc(pt.X),
		Y: math.Trunc(pt.Y),
	}
}

// IsInf reports whether at least one of x and y is infinite.
func (pt Point) IsInf() bool {
	return math.IsInf(pt.X, 0) || math.IsInf(pt.Y, 0)
}

// IsNaN reports whether at least one of x and y is NaN.
func (pt Point) IsNaN() bool {
	return math.IsNaN(pt.X) || math.IsNaN(pt.Y)
}
