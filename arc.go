package curve

import (
	"iter"
	"math"
)

type Arc struct {
	Center     Point
	Radii      Vec2
	StartAngle float64
	SweepAngle float64
	XRotation  float64
}

var _ ClosedShape = Arc{}

// Contains implements ClosedShape.
func (a Arc) Contains(pt Point) bool {
	return a.Winding(pt) != 0
}

func (a Arc) PathElements(tolerance float64) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		p0 := sampleEllipse(a.Radii, a.XRotation, a.StartAngle)
		if !yield(MoveTo(a.Center.Translate(p0))) {
			return
		}

		scaledError := max(a.Radii.X, a.Radii.Y) / tolerance
		// Number of subdivisions per ellipse based on error tolerance.
		// Note: this may slightly underestimate the error for quadrants.
		nError := max(math.Pow(1.1163*scaledError, 1.0/6.0), 3.999_999)
		n := math.Ceil(nError * math.Abs(a.SweepAngle) * (1.0 / (2.0 * math.Pi)))
		angleStep := a.SweepAngle / n
		armLen := math.Copysign((4.0/3.0)*math.Tan(math.Abs(0.25*angleStep)), a.SweepAngle)
		angle0 := a.StartAngle
		p0 = sampleEllipse(a.Radii, a.XRotation, angle0)

		for range int(n) {
			angle1 := angle0 + angleStep
			p1 := p0.Add(sampleEllipse(a.Radii, a.XRotation, angle0+math.Pi/2).Mul(armLen))
			p3 := sampleEllipse(a.Radii, a.XRotation, angle1)
			p2 := p3.Sub(sampleEllipse(a.Radii, a.XRotation, angle1+math.Pi/2).Mul(armLen))

			angle0 = angle1
			p0 = p3

			if !yield(CubicTo(
				a.Center.Translate(p1),
				a.Center.Translate(p2),
				a.Center.Translate(p3),
			)) {
				break
			}
		}
	}
}

// / Take the ellipse radii, how the radii are rotated, and the sweep angle, and return a
// / point on the ellipse.
func sampleEllipse(radii Vec2, xRotation float64, angle float64) Vec2 {
	sin, cos := math.Sincos(angle)
	u := radii.X * cos
	v := radii.Y * sin
	return rotatePt(Vec2{u, v}, xRotation)
}

// / Rotate `pt` about the origin by `angle` radians.
func rotatePt(pt Vec2, angle float64) Vec2 {
	sin, cos := math.Sincos(angle)
	return Vec2{
		X: pt.X*cos - pt.Y*sin,
		Y: pt.X*sin + pt.Y*cos,
	}
}

func (a Arc) Area() float64 {
	return math.Pi * a.Radii.X * a.Radii.Y
}

func (a Arc) BoundingBox() Rect {
	panic("not implemented")
}

// / The perimeter of the ellipse.
// /
// / Note: Finding the perimeter of an ellipse is [fairly involved][wikipedia],
// / so for now we just approximate by using the bezier curve representation.
// /
// / [wikipedia]: https://en.wikipedia.org/wiki/Ellipse#Circumference
func (a Arc) Perimeter(accuracy float64) float64 {
	panic("not implemented")
}

func (a Arc) Winding(pt Point) int {
	panic("not implemented")
}

func (a Arc) Translate(v Vec2) Arc {
	a.Center = a.Center.Translate(v)
	return a
}
