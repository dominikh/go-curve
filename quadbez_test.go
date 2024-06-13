package curve

import (
	"math"
	"testing"

	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestQuadBezArclen(t *testing.T) {
	q := QuadBez{
		Pt(0.0, 0.0),
		Pt(0.0, 0.5),
		Pt(1.0, 1.0),
	}
	want := 0.5*math.Sqrt(5.0) + 0.25*math.Log(2.0+math.Sqrt(5.0))
	for i := range 12 {
		accuracy := math.Pow(0.1, float64(i))
		est := q.Arclen(accuracy)
		error := math.Abs(est - want)
		if error > accuracy {
			t.Errorf("got error %g for desired accuracy of %g", error, accuracy)
		}
	}
}

func TestQuadBezArclenPathological(t *testing.T) {
	q := QuadBez{
		Pt(-1.0, 0.0),
		Pt(1.03, 0.0),
		Pt(1.0, 0.0),
	}
	const want = 2.0008737864167325 // A rough empirical calculation
	const accuracy = 1e-11
	est := q.Arclen(accuracy)
	error := math.Abs(est - want)
	if error > accuracy {
		t.Errorf("got error %g for desired accuracy of %g", error, accuracy)
	}
}

func TestQuadBezSubsegment(t *testing.T) {
	q := QuadBez{
		Pt(3.1, 4.1),
		Pt(5.9, 2.6),
		Pt(5.3, 5.8),
	}
	t0 := 0.1
	t1 := 0.8
	qs := q.Subsegment(t0, t1)
	epsilon := 1e-12
	n := 10
	for i := range n + 1 {
		tt := float64(i) / float64(n)
		ts := t0 + tt*(t1-t0)
		assertNear(t, q.Eval(ts), qs.Eval(tt), epsilon)
	}
}

func TestQuadBezDifferentiate(t *testing.T) {
	q := QuadBez{
		Pt(0.0, 0.0),
		Pt(0.0, 0.5),
		Pt(1.0, 1.0),
	}
	deriv := q.Differentiate()
	const n = 10
	for i := range n + 1 {
		ts := float64(i) / float64(n)
		const delta = 1e-6
		p := q.Eval(ts)
		p1 := q.Eval(ts + delta)
		dApprox := p1.Sub(p).Mul(1.0 / delta)
		d := Vec2(deriv.Eval(ts))
		if error := d.Sub(dApprox).Hypot(); error > delta*2 {
			t.Errorf("got difference of %g, want at most %g", error, delta*2)
		}
	}
}

func TestQuadBezRaise(t *testing.T) {
	q := QuadBez{
		Pt(3.1, 4.1),
		Pt(5.9, 2.6),
		Pt(5.3, 5.8),
	}
	c := q.Raise()
	qd := q.Differentiate()
	cd := c.Differentiate()
	const epsilon = 1e-12
	const n = 10

	for i := range n + 1 {
		ts := float64(i) / float64(n)
		assertNear(t, q.Eval(ts), c.Eval(ts), epsilon)
		assertNear(t, qd.Eval(ts), cd.Eval(ts), epsilon)
	}
}

func TestQuadbezSignedArea(t *testing.T) {
	approxEqual := func(x, y float64) bool {
		return math.Abs(x-y) < 1e-12
	}

	// y = 1 - x^2
	q := QuadBez{Pt(1.0, 0.0), Pt(0.5, 1.0), Pt(0.0, 1.0)}
	if v := q.SignedArea(); !approxEqual(v, 2.0/3.0) {
		t.Errorf("got %v, want %v", v, 2.0/3.0)
	}
	if v := q.Transform(Rotate(0.5)).SignedArea(); !approxEqual(v, 2.0/3.0) {
		t.Errorf("got %v, want %v", v, 2.0/3.0)
	}
	if v := q.Transform(Translate(Vec(0.0, 1.0))).SignedArea(); !approxEqual(v, 3.5/3.0) {
		t.Errorf("got %v, want %v", v, 3.5/3.0)
	}
	if v := q.Transform(Translate(Vec(1.0, 0.0))).SignedArea(); !approxEqual(v, 3.5/3.0) {
		t.Errorf("got %v, want %v", v, 3.5/3.0)
	}
}

func TestQuadbezNearest(t *testing.T) {
	verify := func(q QuadBez, pt Point, want float64) {
		t.Helper()
		_, got := q.Nearest(pt, 1e-3)
		if math.Abs(got-want) > 1e-6 {
			t.Errorf("got %v, want %v", got, want)
		}
	}
	q := QuadBez{Pt(-1.0, 1.0), Pt(0.0, -1.0), Pt(1.0, 1.0)}
	verify(q, Pt(0.0, 0.0), 0.5)
	verify(q, Pt(0.0, 0.1), 0.5)
	verify(q, Pt(0.0, -0.1), 0.5)
	verify(q, Pt(0.5, 0.25), 0.75)
	verify(q, Pt(1.0, 1.0), 1.0)
	verify(q, Pt(1.1, 1.1), 1.0)
	verify(q, Pt(-1.1, 1.1), 0.0)
	a := Rotate(0.5)
	verify(q.Transform(a), Pt(0.5, 0.25).Transform(a), 0.75)
}

func TestQuadbezExtrema(t *testing.T) {
	approx := cmpopts.EquateApprox(0, 1e-6)

	// y = x^2
	q := QuadBez{Pt(-1.0, 1.0), Pt(0.0, -1.0), Pt(1.0, 1.0)}
	extrema, n := q.Extrema()
	want := []float64{0.5}
	diff(t, extrema[:n], want, approx)

	q = QuadBez{Pt(0.0, 0.5), Pt(1.0, 1.0), Pt(0.5, 0.0)}
	extrema, n = q.Extrema()
	want = []float64{1.0 / 3.0, 2.0 / 3.0}
	diff(t, extrema[:n], want, approx)

	// Reverse direction
	q = QuadBez{Pt(0.5, 0.0), Pt(1.0, 1.0), Pt(0.0, 0.5)}
	extrema, n = q.Extrema()
	want = []float64{1.0 / 3.0, 2.0 / 3.0}
	diff(t, extrema[:n], want, approx)
}

func TestIntersectQuad(t *testing.T) {
	q := QuadBez{Pt(0.0, -10.0), Pt(10.0, 20.0), Pt(20.0, -10.0)}
	vLine := Line{Pt(10.0, -10.0), Pt(10.0, 10.0)}
	xs, n := q.IntersectLine(vLine)
	want := []LineIntersection{{0.75, 0.5}}
	diff(t, xs[:n], want, cmpopts.EquateApprox(0, 1e-6))

	hLine := Line{Pt(0.0, 0.0), Pt(100.0, 0.0)}
	if _, n := q.IntersectLine(hLine); n != 2 {
		t.Errorf("got %d intersections, want 2", n)
	}
}

func TestQuadBezNearestLowOrder(t *testing.T) {
	// This test exposes a degenerate case in the solver used internally
	// by the "nearest" calculation - the cubic term is zero.
	verify := func(q QuadBez, pt Point, want float64) {
		t.Helper()
		_, got := q.Nearest(pt, 1e-3)
		if math.Abs(got-want) > 1e-6 {
			t.Errorf("got %v, want %v", got, want)
		}
	}

	q := QuadBez{Pt(-1.0, 0.0), Pt(0.0, 0.0), Pt(1.0, 0.0)}

	verify(q, Pt(0.0, 0.0), 0.5)
	verify(q, Pt(0.0, 1.0), 0.5)
}
