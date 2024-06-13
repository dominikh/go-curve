package curve

import (
	"fmt"
	"math"
	"testing"

	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestCubicBezDeriv(t *testing.T) {
	// y = x^2
	c := CubicBez{
		Pt(0.0, 0.0),
		Pt(1.0/3.0, 0.0),
		Pt(2.0/3.0, 1.0/3.0),
		Pt(1.0, 1.0),
	}
	deriv := c.Differentiate()

	const n = 10
	const delta = 1e-6
	for i := range n + 1 {
		ts := float64(i) / float64(n)
		p := c.Eval(ts)
		p1 := c.Eval(ts + delta)
		dApprox := p1.Sub(p).Mul(1.0 / delta)
		d := Vec2(deriv.Eval(ts))
		if l := d.Sub(dApprox).Hypot(); l >= delta*2 {
			t.Errorf("got difference of %g, want at most %g", l, delta*2)
		}
	}
}

func TestCubicBezToQuadratics(t *testing.T) {
	// y = x^3
	c := CubicBez{
		Pt(0.0, 0.0),
		Pt(1.0/3.0, 0.0),
		Pt(2.0/3.0, 0.0),
		Pt(1.0, 1.0),
	}
	for i := range 10 {
		accuracy := math.Pow(0.1, float64(i))
		worst := 0.0
		for seq := range c.Quadratics(accuracy) {
			t0, t1, q := seq.Start, seq.End, seq.Segment
			const epsilon = 1e-12
			if delta := q.Start().Sub(c.Eval(t0)).Hypot(); delta > epsilon {
				t.Fatalf("%g > %g", delta, epsilon)
			}
			if delta := q.End().Sub(c.Eval(t1)).Hypot(); delta > epsilon {
				t.Fatalf("%g > %g", delta, epsilon)
			}
			const n = 4
			for j := range n + 1 {
				ts := float64(j) / float64(n)
				p := q.Eval(ts)
				error := math.Abs(p.Y - math.Pow(p.X, 3))
				if worst > error {
					error = worst
				}
				if error > accuracy {
					t.Fatalf("got error %g for desired accuracy of %g", error, accuracy)
				}
			}
		}
	}
}

func TestCubicBezToQuadraticsDegenerate(t *testing.T) {
	// ensure ToQuadratics returns something given colinear points
	c := CubicBez{
		Pt(0.0, 9.0),
		Pt(6.0, 6.0),
		Pt(12.0, 3.0),
		Pt(18.0, 0.0),
	}
	var n int
	for range c.Quadratics(1e-6) {
		n++
	}
	if n != 1 {
		t.Errorf("got %d quadratics, expected 1", n)
	}
}

func BenchmarkCubicBezToQuadratics(b *testing.B) {
	shape := CubicBez{
		P0: Pt(20, 40),
		P1: Pt(40, 80),
		P2: Pt(-40, 40),
		P3: Pt(42, 62),
	}
	for i := range 11 {
		acc := 1.0 / math.Pow(10, float64(2*i))
		b.Run(fmt.Sprintf("1e-%d", 2*i), func(b *testing.B) {
			for range b.N {
				for range shape.Quadratics(acc) {
				}
			}
		})
	}
}

func TestIntersectCubic(t *testing.T) {
	c := CubicBez{Pt(0.0, -10.0), Pt(10.0, 20.0), Pt(20.0, -20.0), Pt(30.0, 10.0)}
	vLine := Line{Pt(10.0, -10.0), Pt(10.0, 10.0)}
	xs, n := c.IntersectLine(vLine)
	want := []LineIntersection{{16.0 / 27.0, 1.0 / 3.0}}
	diff(t, want, xs[:n], cmpopts.EquateApprox(0, 1e-8))

	hLine := Line{Pt(0.0, 0.0), Pt(100.0, 0.0)}
	if _, n := c.IntersectLine(hLine); n != 3 {
		t.Errorf("got %d intersections, want 3", n)
	}
}

func TestCubicBezExtrema(t *testing.T) {
	// y = x^2
	q := CubicBez{Pt(0.0, 0.0), Pt(0.0, 1.0), Pt(1.0, 1.0), Pt(1.0, 0.0)}
	extrema, n := q.Extrema()
	if n != 1 {
		t.Fatalf("got %d extrema, expected 1", n)
	}
	if want := 0.5; math.Abs(extrema[0]-want) > 1e-6 {
		t.Errorf("got extrema %v, want %v", extrema[0], want)
	}

	q = CubicBez{Pt(0.4, 0.5), Pt(0.0, 1.0), Pt(1.0, 0.0), Pt(0.5, 0.4)}
	extrema, n = q.Extrema()
	if n != 4 {
		t.Fatalf("got %d extrema, expected 4", n)
	}
}

func TestCubicNearest(t *testing.T) {
	verify := func(c CubicBez, pt Point, want float64) {
		t.Helper()
		_, got := c.Nearest(pt, 1e-6)
		if math.Abs(got-want) > 1e-6 {
			t.Errorf("got %v, want %v", got, want)
		}
	}

	// y = x^3
	c := CubicBez{Pt(0.0, 0.0), Pt(1.0/3.0, 0.0), Pt(2.0/3.0, 0.0), Pt(1.0, 1.0)}
	verify(c, Pt(0.1, 0.001), 0.1)
	verify(c, Pt(0.2, 0.008), 0.2)
	verify(c, Pt(0.3, 0.027), 0.3)
	verify(c, Pt(0.4, 0.064), 0.4)
	verify(c, Pt(0.5, 0.125), 0.5)
	verify(c, Pt(0.6, 0.216), 0.6)
	verify(c, Pt(0.7, 0.343), 0.7)
	verify(c, Pt(0.8, 0.512), 0.8)
	verify(c, Pt(0.9, 0.729), 0.9)
	verify(c, Pt(1.0, 1.0), 1.0)
	verify(c, Pt(1.1, 1.1), 1.0)
	verify(c, Pt(-0.1, 0.0), 0.0)
	a := Rotate(0.5)
	verify(c.Transform(a), Pt(0.1, 0.001).Transform(a), 0.1)
}

func TestCubicBezInflections(t *testing.T) {
	approx := cmpopts.EquateApprox(0, 1e-6)

	c := CubicBez{
		Pt(0.0, 0.0),
		Pt(0.8, 1.0),
		Pt(0.2, 1.0),
		Pt(1.0, 0.0),
	}
	inflections, n := c.Inflections()
	want := []float64{
		0.311018,
		0.688982,
	}
	diff(t, want, inflections[:n], approx)

	c = CubicBez{Pt(0.0, 0.0), Pt(1.0, 1.0), Pt(2.0, -1.0), Pt(3.0, 0.0)}
	inflections, n = c.Inflections()
	want = []float64{0.5}
	diff(t, want, inflections[:n])

	c = CubicBez{Pt(0.0, 0.0), Pt(1.0, 1.0), Pt(2.0, 1.0), Pt(3.0, 0.0)}
	inflections, n = c.Inflections()
	diff(t, []float64{}, inflections[:n])
}

var pointComparer = cmp.Comparer(func(p1, p2 Point) bool {
	return p1.Distance(p2) <= 1e-12
})

func TestCubicBezApproxSpline(t *testing.T) {
	c1 := CubicBez{
		Pt(550.0, 258.0),
		Pt(1044.0, 482.0),
		Pt(2029.0, 1841.0),
		Pt(1934.0, 1554.0),
	}

	{
		quad, _ := c1.tryApproxQuadratic(344.0)
		expected := QuadBez{
			Pt(550.0, 258.0),
			Pt(1673.665720592873, 767.5164401068898),
			Pt(1934.0, 1554.0),
		}
		diff(t, expected, quad)
	}

	{
		_, ok := c1.tryApproxQuadratic(343.0)
		if ok {
			t.Errorf("expected approximation to fail")
		}
	}

	{
		spline, _ := c1.approxQuadSplineN(2, 343.0)
		expected := QuadBSpline{
			Pt(550.0, 258.0),
			Pt(920.5, 426.0),
			Pt(2005.25, 1769.25),
			Pt(1934.0, 1554.0),
		}
		diff(t, expected, spline, pointComparer)
	}

	{
		spline, _ := c1.ApproxQuadSpline(5.0)
		expected := QuadBSpline{
			Pt(550.0, 258.0),
			Pt(673.5, 314.0),
			Pt(984.8777777777776, 584.2666666666667),
			Pt(1312.6305555555557, 927.825),
			Pt(1613.1194444444443, 1267.425),
			Pt(1842.7055555555555, 1525.8166666666666),
			Pt(1957.75, 1625.75),
			Pt(1934.0, 1554.0),
		}
		diff(t, expected, spline, pointComparer)
	}
}

func TestCubicBezCubicsToQuadraticSplines(t *testing.T) {
	curves := []CubicBez{
		CubicBez{
			Pt(550.0, 258.0),
			Pt(1044.0, 482.0),
			Pt(2029.0, 1841.0),
			Pt(1934.0, 1554.0),
		},
		CubicBez{
			Pt(859.0, 384.0),
			Pt(1998.0, 116.0),
			Pt(1596.0, 1772.0),
			Pt(8.0, 1824.0),
		},
		CubicBez{
			Pt(1090.0, 937.0),
			Pt(418.0, 1300.0),
			Pt(125.0, 91.0),
			Pt(104.0, 37.0),
		},
	}
	converted, ok := CubicsToQuadraticSplines(curves, 5.0)
	if !ok {
		t.Fatal("could not convert cubics to b-splines")
	}
	if len(converted) != 3 {
		t.Fatalf("got %d splines, want 3", len(converted))
	}
	for i, s := range converted {
		if len(s) != 8 {
			t.Fatalf("got %d points in spline %d, want 8", len(s), i)
		}
	}
	diff(t, Pt(673.5, 314.0), converted[0][1], pointComparer)
	diff(t, Pt(88639.0/90.0, 52584.0/90.0), converted[0][2], pointComparer)
}

func TestCubicBezApproxSplineMatchesPython(t *testing.T) {
	// Ensure rounding behavior for division matches fonttools
	// cu2qu.
	// See https://github.com/linebender/kurbo/issues/272
	cubic := CubicBez{
		Pt(408.0, 321.0),
		Pt(408.0, 452.0),
		Pt(342.0, 560.0),
		Pt(260.0, 560.0),
	}
	spline, ok := cubic.ApproxQuadSpline(1.0)
	if !ok {
		t.Fatal("could not convert cubic to b-spline")
	}
	diff(
		t,
		QuadBSpline{
			Pt(408.0, 321.0),
			// Previous behavior produced 386.49999999999994 for the
			// y coordinate leading to inconsistent rounding.
			Pt(408.0, 386.5),
			Pt(368.16666666666663, 495.0833333333333),
			Pt(301.0, 560.0),
			Pt(260.0, 560.0),
		},
		spline,
		/* not using pointComparer, we want an exact match */
	)
}

func TestCubicToQuadraticMatchesPython(t *testing.T) {
	// from https://github.com/googlefonts/fontmake-rs/issues/217
	cubic := CubicBez{
		Pt(796.0, 319.0),
		Pt(727.0, 314.0),
		Pt(242.0, 303.0),
		Pt(106.0, 303.0),
	}

	// FontTools can approximate this curve successfully in 7 splits, we can too
	if _, ok := cubic.approxQuadSplineN(7, 1.0); !ok {
		t.Error("could not approximate curve in 7 splits")
	}

	// FontTools can solve this with accuracy 0.001, we can too
	if _, ok := CubicsToQuadraticSplines([]CubicBez{cubic}, 0.001); !ok {
		t.Error("could not approximate curve with 0.001 accuracy")
	}
}

func TestCubicsToQuadraticSplinesMatchesPython(t *testing.T) {
	// https://github.com/linebender/kurbo/pull/273
	light := CubicBez{
		Pt(378.0, 608.0),
		Pt(378.0, 524.0),
		Pt(355.0, 455.0),
		Pt(266.0, 455.0),
	}
	regular := CubicBez{
		Pt(367.0, 607.0),
		Pt(367.0, 511.0),
		Pt(338.0, 472.0),
		Pt(243.0, 472.0)}
	bold := CubicBez{
		Pt(372.425, 593.05),
		Pt(372.425, 524.95),
		Pt(355.05, 485.95),
		Pt(274.0, 485.95),
	}
	qsplines, _ := CubicsToQuadraticSplines([]CubicBez{light, regular, bold}, 1.0)
	diff(
		t,
		[]QuadBSpline{
			QuadBSpline{
				Pt(378.0, 608.0),
				Pt(378.0, 566.0),
				Pt(359.0833333333333, 496.5),
				Pt(310.5, 455.0),
				Pt(266.0, 455.0),
			},
			QuadBSpline{
				Pt(367.0, 607.0),
				Pt(367.0, 559.0),
				// Previous behavior produced 496.5 for the y coordinate
				Pt(344.5833333333333, 499.49999999999994),
				Pt(290.5, 472.0),
				Pt(243.0, 472.0),
			},
			QuadBSpline{
				Pt(372.425, 593.05),
				Pt(372.425, 559.0),
				Pt(356.98333333333335, 511.125),
				Pt(314.525, 485.95),
				Pt(274.0, 485.95),
			},
		},
		qsplines,
	)
}

func TestCubicBezArclen(t *testing.T) {
	c := CubicBez{
		Pt(0.0, 0.0),
		Pt(1.0/3.0, 0.0),
		Pt(2.0/3.0, 1.0/3.0),
		Pt(1.0, 1.0),
	}
	trueArclen := 0.5*math.Sqrt(5.0) + 0.25*math.Log(2.0+math.Sqrt(5.0))
	for i := range 12 {
		accuracy := math.Pow(0.1, float64(i))
		diff(t, trueArclen, c.Arclen(accuracy), cmpopts.EquateApprox(0, accuracy))
	}
}

func TestCubicBezInvArclen(t *testing.T) {
	// y = x^2 / 100
	c := CubicBez{
		Pt(0.0, 0.0),
		Pt(100.0/3.0, 0.0),
		Pt(200.0/3.0, 100.0/3.0),
		Pt(100.0, 100.0),
	}
	trueArclen := 100.0 * (0.5*math.Sqrt(5.0) + 0.25*math.Log(2.0+math.Sqrt(5.0)))
	for i := range 12 {
		accuracy := math.Pow(0.1, float64(i))
		n := 10
		for j := range n + 1 {
			arc := float64(j) * (1.0 / float64(n) * trueArclen)
			tt := SolveForArclen(c, arc, accuracy*0.5)
			actualArc := c.Subsegment(0.0, tt).Arclen(accuracy * 0.5)
			diff(t, arc, actualArc, cmpopts.EquateApprox(0, accuracy))
		}
	}
	// corner case: user passes accuracy larger than total arc length
	accuracy := trueArclen * 1.1
	arc := trueArclen * 0.5
	tt := SolveForArclen(c, arc, accuracy)
	actualArc := c.Subsegment(0.0, tt).Arclen(accuracy)
	diff(t, arc, actualArc, cmpopts.EquateApprox(0, 2*accuracy))
}

func TestCubicBezInvArclenAccuracy(t *testing.T) {
	c := CubicBez{
		Pt(0.2, 0.73),
		Pt(0.35, 1.08),
		Pt(0.85, 1.08),
		Pt(1.0, 0.73),
	}
	trueT := SolveForArclen(c, 0.5, 1e-12)
	for i := 1; i < 12; i++ {
		accuracy := math.Pow(0.1, float64(i))
		approxT := SolveForArclen(c, 0.5, accuracy)
		diff(t, trueT, approxT, cmpopts.EquateApprox(0, accuracy))
	}
}

func TestCubicBezSignedAreaLinear(t *testing.T) {
	// y = 1 - x
	c := CubicBez{
		Pt(1.0, 0.0),
		Pt(2.0/3.0, 1.0/3.0),
		Pt(1.0/3.0, 2.0/3.0),
		Pt(0.0, 1.0),
	}
	const epsilon = 1e-12

	diff(t, 0.5, c.SignedArea())
	diff(t, 0.5, c.Transform(Rotate(0.5)).SignedArea(), cmpopts.EquateApprox(0, epsilon))
	diff(t, 1.0, c.Transform(Translate(Vec(0.0, 1.0))).SignedArea(), cmpopts.EquateApprox(0, epsilon))
	diff(t, 1.0, c.Transform(Translate(Vec(1.0, 0.0))).SignedArea(), cmpopts.EquateApprox(0, epsilon))
}

func TestCubicBezSignedArea(t *testing.T) {
	// y = 1 - x^3
	c := CubicBez{
		Pt(1.0, 0.0),
		Pt(2.0/3.0, 1.0),
		Pt(1.0/3.0, 1.0),
		Pt(0.0, 1.0),
	}
	const epsilon = 1e-12
	diff(t, 0.75, c.SignedArea(), cmpopts.EquateApprox(0, epsilon))
	diff(t, 0.75, c.Transform(Rotate(0.5)).SignedArea(), cmpopts.EquateApprox(0, epsilon))
	diff(t, 1.25, c.Transform(Translate(Vec(0.0, 1.0))).SignedArea(), cmpopts.EquateApprox(0, epsilon))
	diff(t, 1.25, c.Transform(Translate(Vec(1.0, 0.0))).SignedArea(), cmpopts.EquateApprox(0, epsilon))
}
