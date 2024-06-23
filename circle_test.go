package curve

import (
	"math"
	"slices"
	"testing"
)

func TestCircleAreaSign(t *testing.T) {
	approxEqual := func(x, y float64) bool {
		return math.Abs(x-y) < 1e-7
	}

	center := Pt(5, 5)
	c := Circle{center, 5}
	if a := c.Area(); !approxEqual(a, 25*math.Pi) {
		t.Errorf("got area %v, expected %v", a, 25.0*math.Pi)
	}

	if w := c.Winding(center); w != 1 {
		t.Errorf("got winding number %d, expected 1", w)
	}

	p := BezPath(slices.Collect(c.PathElements(1e-9)))
	if ca, pa := c.Area(), p.SignedArea(); !approxEqual(ca, pa) {
		t.Errorf("got areas %v and %v, expected them to be equal", ca, pa)
	}
	if cw, pw := c.Winding(center), p.Winding(center); cw != pw {
		t.Errorf("got winding numbers %d and %d, expected them to be equal", cw, pw)
	}

	cNegRadius := Circle{center, -5}
	if a := cNegRadius.Area(); !approxEqual(a, 25.0*math.Pi) {
		t.Errorf("got area %v, expected %v", a, 25.0*math.Pi)
	}

	if w := cNegRadius.Winding(center); w != 1 {
		t.Errorf("got winding number %d, expected 1", w)
	}

	pNegRadius := BezPath(slices.Collect(cNegRadius.PathElements(1e-9)))
	if ca, pa := cNegRadius.Area(), pNegRadius.SignedArea(); !approxEqual(ca, pa) {
		t.Errorf("got areas %v and %v, expected them to be equal", ca, pa)
	}
	if cw, pw := cNegRadius.Winding(center), pNegRadius.Winding(center); cw != pw {
		t.Errorf("got winding numbers %d and %d, expected them to be equal", cw, pw)
	}
}
