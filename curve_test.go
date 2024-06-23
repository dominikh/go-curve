package curve

import (
	"math"
	"slices"
	"sort"
	"testing"
)

func checkRoots(t *testing.T, roots, expected []float64) {
	if len(roots) != len(expected) {
		t.Fatalf("got %d roots, expected %d", len(roots), len(expected))
	}
	const epsilon = 1e-12
	sort.Float64s(roots)
	sort.Float64s(expected)
	for i := range roots {
		if math.Abs(roots[i]-expected[i]) > epsilon {
			t.Errorf("root %d is %v but we expected %v", i, roots[i], expected[i])
		}
	}
}

func TestSolveCubic(t *testing.T) {
	slice := func(roots [3]float64, n int) []float64 {
		return roots[:n]
	}
	checkRoots(t, slice(SolveCubic(-5, 0, 0, 1)), []float64{math.Cbrt(5)})
	checkRoots(t, slice(SolveCubic(-5.0, -1.0, 0.0, 1.0)), []float64{1.90416085913492})
	checkRoots(t, slice(SolveCubic(0.0, -1.0, 0.0, 1.0)), []float64{-1.0, 0.0, 1.0})
	checkRoots(t, slice(SolveCubic(-2.0, -3.0, 0.0, 1.0)), []float64{-1.0, 2.0})
	checkRoots(t, slice(SolveCubic(2.0, -3.0, 0.0, 1.0)), []float64{-2.0, 1.0})
	checkRoots(t, slice(SolveCubic(2.0-1e-12, 5.0, 4.0, 1.0)),
		[]float64{
			-1.9999999999989995,
			-1.0000010000848456,
			-0.9999989999161546,
		},
	)
	checkRoots(t, slice(SolveCubic(2.0+1e-12, 5.0, 4.0, 1.0)), []float64{-2.0})
}

func TestSolveQuadratic(t *testing.T) {
	slice := func(roots [2]float64, n int) []float64 {
		return roots[:n]
	}
	checkRoots(t, slice(SolveQuadratic(-5.0, 0.0, 1.0)), []float64{-math.Sqrt(5), math.Sqrt(5)})
	checkRoots(t, slice(SolveQuadratic(5.0, 0.0, 1.0)), []float64{})
	checkRoots(t, slice(SolveQuadratic(5.0, 1.0, 0.0)), []float64{-5.0})
	checkRoots(t, slice(SolveQuadratic(1.0, 2.0, 1.0)), []float64{-1.0})
}

func TestSolveQuartic(t *testing.T) {
	// These test cases are taken from Orellana and De Michele paper (Table 1).
	testWithRoots := func(coeffs [4]float64, roots []float64, relErr float64) {
		t.Helper()

		// Note: in paper, coefficients are in decreasing order.
		actual, n := SolveQuartic(coeffs[3], coeffs[2], coeffs[1], coeffs[0], 1.0)
		sort.Float64s(actual[:n])
		if n != len(roots) {
			t.Fatalf("got %d roots, expected %d", n, len(roots))
		}
		for i := range actual[:n] {
			if math.Abs(actual[i]-roots[i]) > relErr*math.Abs(roots[i]) {
				t.Errorf("root %d is %v but we expected %v", i, actual[i], roots[i])
			}
		}
	}

	testVietaRoots := func(x1, x2, x3, x4 float64, roots []float64, relErr float64) {
		t.Helper()
		a := -(x1 + x2 + x3 + x4)
		b := x1*(x2+x3) + x2*(x3+x4) + x4*(x1+x3)
		c := -x1*x2*(x3+x4) - x3*x4*(x1+x2)
		d := x1 * x2 * x3 * x4
		testWithRoots([4]float64{a, b, c, d}, roots, relErr)
	}

	testVieta := func(x1, x2, x3, x4, relErr float64) {
		t.Helper()
		testVietaRoots(x1, x2, x3, x4, []float64{x1, x2, x3, x4}, relErr)
	}

	// case 1
	testVieta(1.0, 1e3, 1e6, 1e9, 1e-16)
	// case 2
	testVieta(2.0, 2.001, 2.002, 2.003, 1e-6)
	// case 3
	testVieta(1e47, 1e49, 1e50, 1e53, 2e-16)
	// case 4
	testVieta(-1.0, 1.0, 2.0, 1e14, 1e-16)
	// case 5
	testVieta(-2e7, -1.0, 1.0, 1e7, 1e-16)
	// case 6
	testWithRoots(
		[4]float64{-9000002.0, -9999981999998.0, 19999982e6, -2e13},
		[]float64{-1e6, 1e7},
		1e-16,
	)
	// case 7
	testWithRoots(
		[4]float64{2000011.0, 1010022000028.0, 11110056e6, 2828e10},
		[]float64{-7.0, -4.0},
		1e-16,
	)
	// case 8
	testWithRoots(
		[4]float64{-100002011.0, 201101022001.0, -102200111000011.0, 11000011e8},
		[]float64{11.0, 1e8},
		1e-16,
	)
	// cases 9-13 have no real roots
	// case 14
	testVietaRoots(1000.0, 1000.0, 1000.0, 1000.0, []float64{1000.0, 1000.0}, 1e-16)
	// case 15
	testVietaRoots(1e-15, 1000.0, 1000.0, 1000.0, []float64{1e-15, 1000.0, 1000.0}, 1e-15)
	// case 16 no real roots
	// case 17
	testVieta(10000.0, 10001.0, 10010.0, 10100.0, 1e-6)
	// case 19
	testVietaRoots(1.0, 1e30, 1e30, 1e44, []float64{1.0, 1e30, 1e44}, 1e-16)
	// case 20
	// FAILS, error too big
	testVieta(1.0, 1e7, 1e7, 1e14, 1e-7)
	// case 21 doesn't pick up double root
	// case 22
	testVieta(1.0, 10.0, 1e152, 1e154, 3e-16)
	// case 23
	testWithRoots(
		[4]float64{1.0, 1.0, 3.0 / 8.0, 1e-3},
		[]float64{-0.497314148060048, -0.00268585193995149},
		2e-15,
	)
	// case 24
	const s = 1e30
	testWithRoots(
		[4]float64{-(1.0 + 1.0/s), 1.0/s - s*s, s*s + s, -s},
		[]float64{-s, 1e-30, 1.0, s},
		2e-16,
	)
}

func TestSolveITP(t *testing.T) {
	f := func(x float64) float64 { return x*x*x - x - 2.0 }
	x := SolveITP(f, 1.0, 2.0, 1e-12, 0, 0.2, f(1.0), f(2.0))
	if n := math.Abs(f(x)); n > 6e-12 {
		t.Errorf("%v > 6e-12", n)
	}
}

func TestSolveForArclen(t *testing.T) {
	c := CubicBez{
		Pt(0.0, 0.0),
		Pt(100.0/3.0, 0.0),
		Pt(200.0/3.0, 100.0/3.0),
		Pt(100.0, 100.0),
	}
	const target = 100.0
	SolveITP(
		func(t float64) float64 { return c.Subsegment(0.0, t).Arclen(1e-9) - target },
		0.0,
		1.0,
		1e-6,
		1,
		0.2,
		-target,
		c.Arclen(1e-9)-target,
	)
}

func TestSVGSingle(t *testing.T) {
	segments := []PathSegment{
		CubicBez{
			Pt(10.0, 10.0),
			Pt(20.0, 20.0),
			Pt(30.0, 30.0),
			Pt(40.0, 40.0),
		}.Seg(),
	}
	var path BezPath = slices.Collect(Elements(slices.Values(segments)))
	want := "M10,10 C20,20 30,30 40,40"
	got := path.SVG(SVGOptions{})
	diff(t, got, want)
}

func TestSVGTwoNoMove(t *testing.T) {
	segments := []PathSegment{
		CubicBez{
			Pt(10.0, 10.0),
			Pt(20.0, 20.0),
			Pt(30.0, 30.0),
			Pt(40.0, 40.0),
		}.Seg(),
		CubicBez{
			Pt(40.0, 40.0),
			Pt(30.0, 30.0),
			Pt(20.0, 20.0),
			Pt(10.0, 10.0),
		}.Seg(),
	}
	var path BezPath = slices.Collect(Elements(slices.Values(segments)))
	want := "M10,10 C20,20 30,30 40,40 C30,30 20,20 10,10"
	got := path.SVG(SVGOptions{})
	diff(t, got, want)
}

func TestSVGTwoMove(t *testing.T) {
	segments := []PathSegment{
		CubicBez{
			Pt(10.0, 10.0),
			Pt(20.0, 20.0),
			Pt(30.0, 30.0),
			Pt(40.0, 40.0),
		}.Seg(),
		CubicBez{
			Pt(50.0, 50.0),
			Pt(30.0, 30.0),
			Pt(20.0, 20.0),
			Pt(10.0, 10.0),
		}.Seg(),
	}
	var path BezPath = slices.Collect(Elements(slices.Values(segments)))
	want := "M10,10 C20,20 30,30 40,40 M50,50 C30,30 20,20 10,10"
	got := path.SVG(SVGOptions{})
	diff(t, got, want)
}
