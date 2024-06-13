package curve

import (
	"math"
	"testing"
)

func assertNear(t *testing.T, p0 Point, p1 Point, epsilon float64) {
	t.Helper()
	if d := p1.Sub(p0).Hypot(); d > epsilon {
		t.Fatalf("got %s, expected %s", p0, p1)
	}
}

func TestAffineBasic(t *testing.T) {
	const epsilon = 1e-9
	p := Pt(3, 4)

	assertNear(t, p.Transform(Identity), p, epsilon)
	assertNear(t, p.Transform(Scale(2, 2)), Pt(6, 8), epsilon)
	assertNear(t, p.Transform(Rotate(0)), p, epsilon)
	assertNear(t, p.Transform(Rotate(math.Pi/2)), Pt(-4, 3), epsilon)
	assertNear(t, p.Transform(Translate(Vec(5, 6))), Pt(8, 10), epsilon)
	assertNear(t, p.Transform(Skew(0, 0)), p, epsilon)
	assertNear(t, p.Transform(Skew(2, 4)), Pt(11, 16), epsilon)
}

func TestAffineMul(t *testing.T) {
	const epsilon = 1e-9
	a1 := Affine{1, 2, 3, 4, 5, 6}
	a2 := Affine{0.1, 1.2, 2.3, 3.4, 4.5, 5.6}

	px := Pt(1, 0)
	py := Pt(0, 1)
	pxy := Pt(1, 1)

	assertNear(t, px.Transform(a2).Transform(a1), px.Transform(a1.Mul(a2)), epsilon)
	assertNear(t, py.Transform(a2).Transform(a1), py.Transform(a1.Mul(a2)), epsilon)
	assertNear(t, pxy.Transform(a2).Transform(a1), pxy.Transform(a1.Mul(a2)), epsilon)
}

func TestAffineInvert(t *testing.T) {
	const epsilon = 1e-9
	a := Affine{0.1, 1.2, 2.3, 3.4, 4.5, 5.6}
	aInv := a.Invert()

	px := Pt(1, 0)
	py := Pt(0, 1)
	pxy := Pt(1, 1)

	assertNear(t, px.Transform(aInv).Transform(a), px, epsilon)
	assertNear(t, py.Transform(aInv).Transform(a), py, epsilon)
	assertNear(t, pxy.Transform(aInv).Transform(a), pxy, epsilon)
	assertNear(t, px.Transform(a).Transform(aInv), px, epsilon)
	assertNear(t, py.Transform(a).Transform(aInv), py, epsilon)
	assertNear(t, pxy.Transform(a).Transform(aInv), pxy, epsilon)
}

func TestReflection(t *testing.T) {
	affineAssertNear := func(a0, a1 Affine) {
		a0a := a0.Coefficients()
		a1a := a1.Coefficients()
		for i := range 6 {
			if d := math.Abs(a0a[i] - a1a[i]); d > 1e-9 {
				t.Fatalf("%g > %g", d, 1e-9)
			}
		}
	}

	affineAssertNear(Reflect(Point{}, Vec(1, 0)), Affine{1, 0, 0, -1, 0, 0})
	affineAssertNear(Reflect(Point{}, Vec(0, 1)), Affine{-1, 0, 0, 1, 0, 0})
	affineAssertNear(Reflect(Point{}, Vec(1, 1)), Affine{0, 1, 1, 0, 0, 0})

	const epsilon = 1e-9
	{
		// No translation
		point := Pt(0, 0)
		vec := Vec(1, 1)
		aff := Reflect(point, vec)
		assertNear(t, Pt(0, 0).Transform(aff), Pt(0, 0), epsilon)
		assertNear(t, Pt(1, 1).Transform(aff), Pt(1, 1), epsilon)
		assertNear(t, Pt(1, 2).Transform(aff), Pt(2, 1), epsilon)
	}

	{
		// With translation
		point := Pt(1, 0)
		vec := Vec(1, 1)
		aff := Reflect(point, vec)
		assertNear(t, Pt(1, 0).Transform(aff), Pt(1, 0), epsilon)
		assertNear(t, Pt(2, 1).Transform(aff), Pt(2, 1), epsilon)
		assertNear(t, Pt(2, 2).Transform(aff), Pt(3, 1), epsilon)
	}
}
