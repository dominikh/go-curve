package curve

import (
	"testing"
)

func TestPointArithmetic(t *testing.T) {
	diff(t, Pt(0, 0).Translate(Vec(-10, 0)), Pt(-10, 0))
}

func TestPointDistance(t *testing.T) {
	p1 := Pt(0, 10)
	p2 := Pt(0, 5)
	if d := p1.Distance(p2); d != 5 {
		t.Errorf("got distance %v, want 5", d)
	}

	p3 := Pt(-11, 1)
	p4 := Pt(-7, -2)
	if d := p3.Distance(p4); d != 5 {
		t.Errorf("got distance %v, want 5", d)
	}
}
