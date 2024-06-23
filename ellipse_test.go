package curve

import (
	"math"
	"slices"
	"testing"
)

func TestEllipseAreaSign(t *testing.T) {
	approxEqual := func(x, y float64) bool {
		return math.Abs(x-y) < 1e-7
	}

	center := Pt(5.0, 5.0)
	e := NewEllipse(center, Vec(5.0, 5.0), 1.0)
	if a := e.Area(); !approxEqual(a, 25.0*math.Pi) {
		t.Errorf("got area %v, expected %v", a, 25.0*math.Pi)
	}
	e = NewEllipse(center, Vec(5.0, 10.0), 1.0)
	if a := e.Area(); !approxEqual(a, 50.0*math.Pi) {
		t.Errorf("got area %v, expected %v", a, 50.0*math.Pi)
	}

	if w := e.Winding(center); w != 1 {
		t.Errorf("got winding number %d, expected 1", w)
	}

	p := BezPath(slices.Collect(e.PathElements(1e-9)))
	if ea, pa := e.Area(), p.SignedArea(); !approxEqual(ea, pa) {
		t.Errorf("got areas %v and %v, expected them to be the same", ea, pa)
	}

	if ew, pw := e.Winding(center), p.Winding(center); ew != pw {
		t.Errorf("got winding numbers %d and %d, expected them to be the same", ew, pw)
	}

	eNegRadius := NewEllipse(center, Vec(-5.0, 10.0), 1.0)
	if a := eNegRadius.Area(); !approxEqual(a, 50.0*math.Pi) {
		t.Errorf("got area %v, expected %v", a, 50.0*math.Pi)
	}

	if w := eNegRadius.Winding(center); w != 1 {
		t.Errorf("got winding number %d, expected 1", w)
	}

	pNegRadius := BezPath(slices.Collect(eNegRadius.PathElements(1e-9)))
	if ea, pa := eNegRadius.Area(), pNegRadius.SignedArea(); !approxEqual(ea, pa) {
		t.Errorf("got areas %v and %v, expected them to be the same", ea, pa)
	}

	if ew, pw := eNegRadius.Winding(center), pNegRadius.Winding(center); ew != pw {
		t.Errorf("got winding numbers %d and %d, expected them to be the same", ew, pw)
	}
}
