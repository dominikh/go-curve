package curve

import (
	"iter"
	"math"
	"slices"
)

type Circle struct {
	Center Point
	Radius float64
}

var _ ClosedShape = Circle{}

// Contains implements ClosedShape.
func (c Circle) Contains(pt Point) bool {
	return c.Winding(pt) != 0
}

func (c Circle) Path(tolerance float64) BezPath { return slices.Collect(c.PathElements(tolerance)) }

func (c Circle) PathElements(tolerance float64) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		scaledError := math.Abs(c.Radius) / tolerance
		var n int
		var armLength float64
		if scaledError < 1.0/1.9608e-4 {
			// Solution from http://spencermortensen.com/articles/bezier-circle/
			n = 4
			armLength = 0.551915024494
		} else {
			// This is empirically determined to fall within error tolerance.
			n = int(math.Ceil(math.Pow(1.1163*scaledError, 1.0/6.0)))
			// Note: this isn't minimum error, but it is simple and we can easily
			// estimate the error.
			armLength = (4.0 / 3.0) * math.Tan(math.Pi/2/(float64(n)))
		}

		x, y := c.Center.Splat()
		r := c.Radius
		if !yield(MoveTo(Pt(x+r, y))) {
			return
		}
		deltaTh := 2.0 * math.Pi / float64(n)
		for ix := 1; ix <= n; ix++ {
			a := armLength
			th1 := deltaTh * float64(ix)
			th0 := th1 - deltaTh
			s0, c0 := math.Sincos(th0)
			var s1, c1 float64
			if ix == n {
				s1 = 0.0
				c1 = 1.0
			} else {
				s1, c1 = math.Sincos(th1)
			}
			if !yield(CubicTo(
				Pt(x+r*(c0-a*s0), y+r*(s0+a*c0)),
				Pt(x+r*(c1+a*s1), y+r*(s1-a*c1)),
				Pt(x+r*c1, y+r*s1),
			)) {
				return
			}
		}
		if !yield(ClosePath()) {
			return
		}
	}
}

// Segment returns a circle segment by cutting out parts of this circle.
func (c Circle) Segment(innerRadius float64, startAngle, sweepAngle float64) CircleSegment {
	return CircleSegment{
		Center:      c.Center,
		OuterRadius: c.Radius,
		InnerRadius: innerRadius,
		StartAngle:  startAngle,
		SweepAngle:  sweepAngle,
	}
}

func (c Circle) IsInf() bool {
	return c.Center.IsInf() || math.IsInf(c.Radius, 0)
}

func (c Circle) IsNaN() bool {
	return c.Center.IsNaN() || math.IsNaN(c.Radius)
}

func (c Circle) Translate(v Vec2) Circle {
	return Circle{
		Center: c.Center.Translate(v),
		Radius: c.Radius,
	}
}

func (c Circle) Area() float64 {
	return math.Pi * c.Radius * c.Radius
}

func (c Circle) BoundingBox() Rect {
	r := math.Abs(c.Radius)
	x := c.Center.X
	y := c.Center.Y
	return Rect{
		X0: x - r,
		Y0: y - r,
		X1: x + r,
		Y1: y + r,
	}
}

func (c Circle) Perimeter(accuracy float64) float64 {
	return math.Abs(2 * math.Pi * c.Radius)
}

func (c Circle) Winding(pt Point) int {
	if pt.Sub(c.Center).Hypot2() < c.Radius*c.Radius {
		return 1
	} else {
		return 0
	}
}

func (c Circle) Transform(aff Affine) Ellipse {
	return NewEllipseFromCircle(c).Transform(aff)
}

// CircleSegment represents a segment of a circle.
//
// If InnerRadius > 0, then the shape will be a doughnut segment.
type CircleSegment struct {
	Center      Point
	OuterRadius float64
	InnerRadius float64
	StartAngle  float64
	SweepAngle  float64
}

var _ ClosedShape = CircleSegment{}

// Contains implements ClosedShape.
func (cs CircleSegment) Contains(pt Point) bool {
	return cs.Winding(pt) != 0
}

func (cs CircleSegment) Path(tolerance float64) BezPath {
	return slices.Collect(cs.PathElements(tolerance))
}

// PathElements implements Shape.
func (cs CircleSegment) PathElements(tolerance float64) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		if !yield(MoveTo(pointOnCircle(cs.Center, cs.InnerRadius, cs.StartAngle))) {
			return
		}

		// First radius
		if !yield(LineTo(pointOnCircle(cs.Center, cs.OuterRadius, cs.StartAngle))) {
			return
		}

		// Outer arc
		a := Arc{
			Center:     cs.Center,
			Radii:      Vec2{cs.OuterRadius, cs.OuterRadius},
			StartAngle: cs.StartAngle,
			SweepAngle: cs.SweepAngle,
			XRotation:  0.0,
		}
		for el := range dropFirst(a.PathElements(tolerance)) {
			if !yield(el) {
				return
			}
		}

		// Second radius
		if !yield(LineTo(pointOnCircle(cs.Center, cs.InnerRadius, cs.StartAngle+cs.SweepAngle))) {
			return
		}

		// Inner arc
		a = Arc{
			Center:     cs.Center,
			Radii:      Vec2{cs.InnerRadius, cs.InnerRadius},
			StartAngle: cs.StartAngle + cs.SweepAngle,
			SweepAngle: -cs.SweepAngle,
			XRotation:  0.0,
		}
		for el := range dropFirst(a.PathElements(tolerance)) {
			if !yield(el) {
				return
			}
		}
	}
}

func pointOnCircle(center Point, radius float64, angle float64) Point {
	sin, cos := math.Sincos(angle)
	return center.Translate(
		Vec2{
			X: cos * radius,
			Y: sin * radius,
		})
}

func (cs CircleSegment) IsInf() bool {
	return cs.Center.IsInf() ||
		math.IsInf(cs.OuterRadius, 0) ||
		math.IsInf(cs.InnerRadius, 0) ||
		math.IsInf(cs.StartAngle, 0) ||
		math.IsInf(cs.SweepAngle, 0)
}

func (cs CircleSegment) IsNaN() bool {
	return cs.Center.IsNaN() ||
		math.IsNaN(cs.OuterRadius) ||
		math.IsNaN(cs.InnerRadius) ||
		math.IsNaN(cs.StartAngle) ||
		math.IsNaN(cs.SweepAngle)
}

func (cs CircleSegment) Translate(v Vec2) CircleSegment {
	cs.Center = cs.Center.Translate(v)
	return cs
}

func (cs CircleSegment) Area() float64 {
	return 0.5 * math.Abs(cs.OuterRadius*cs.OuterRadius-cs.InnerRadius*cs.InnerRadius) * cs.SweepAngle
}

func (cs CircleSegment) BoundingBox() Rect {
	// todo this is currently not tight
	r := max(cs.InnerRadius, cs.OuterRadius)
	x := cs.Center.X
	y := cs.Center.Y
	return Rect{
		X0: x - r,
		Y0: y - r,
		X1: x + r,
		Y1: y + r,
	}
}

func (cs CircleSegment) Perimeter(accuracy float64) float64 {
	return 2.0*
		math.Abs(cs.OuterRadius-cs.InnerRadius) +
		cs.SweepAngle*(cs.InnerRadius+cs.OuterRadius)
}

func (cs CircleSegment) Winding(pt Point) int {
	angle := pt.Sub(cs.Center).Angle()
	if angle < cs.StartAngle || angle > cs.StartAngle+cs.SweepAngle {
		return 0
	}
	dist2 := pt.Sub(cs.Center).Hypot2()
	if dist2 < cs.OuterRadius*cs.OuterRadius && dist2 > cs.InnerRadius*cs.InnerRadius ||

		dist2 < cs.InnerRadius*cs.InnerRadius && dist2 > cs.OuterRadius*cs.OuterRadius {
		return 1
	} else {
		return 0
	}
}
