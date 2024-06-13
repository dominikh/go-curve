package curve

import (
	"math"
	"testing"

	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestLineArclen(t *testing.T) {
	l := Line{Pt(0.0, 0.0), Pt(1.0, 1.0)}
	want := math.Sqrt(2.0)
	epsilon := 1e-9
	if d := l.Arclen(epsilon) - want; d > epsilon {
		t.Errorf("%g > %g", d, epsilon)
	}

	ts := l.SolveForArclen(want/3.0, epsilon)
	if d := math.Abs(ts - 1.0/3.0); d > epsilon {
		t.Errorf("%g > %g", d, epsilon)
	}
}

func TestLineIsInf(t *testing.T) {
	if (Line{Pt(0.0, 0.0), Pt(1.0, 1.0)}).IsInf() {
		t.Error("line is infinite but shouldn't be")
	}

	if !(Line{Pt(0.0, 0.0), Pt(math.Inf(1), 1.0)}).IsInf() {
		t.Errorf("line is finite but shouldn't be")
	}

	if !(Line{Pt(0.0, 0.0), Pt(0.0, math.Inf(1))}).IsInf() {
		t.Errorf("line is finite but shouldn't be")
	}
}

func TestIntersectLine(t *testing.T) {
	hLine := Line{Pt(0.0, 0.0), Pt(100.0, 0.0)}
	vLine := Line{Pt(10.0, -10.0), Pt(10.0, 10.0)}
	xs, n := hLine.IntersectLine(vLine)
	want := []LineIntersection{{0.5, 0.1}}
	diff(t, xs[:n], want, cmpopts.EquateApprox(0, 1e-7))

	vLine = Line{Pt(-10.0, -10.0), Pt(-10.0, 10.0)}
	if xs, n := hLine.IntersectLine(vLine); n != 0 {
		t.Errorf("expected no intersections, got %v", xs[:n])
	}

	vLine = Line{Pt(10.0, 10.0), Pt(10.0, 20.0)}
	if xs, n := hLine.IntersectLine(vLine); n != 0 {
		t.Errorf("expected no intersections, got %v", xs[:n])
	}
}
