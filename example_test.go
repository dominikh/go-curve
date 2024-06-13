package curve_test

import (
	"fmt"

	"honnef.co/go/curve"
)

type Quartic struct{}

var _ curve.FittableCurve = Quartic{}

// eval evaluates f(x) = (x⁴ + x³ - 13x² - x) / 10
func eval(x float64) float64 {
	return (x*x*x*x + x*x*x - 13*x*x - x) / 10
}

// evalDeriv evaluates the derivative of [eval], f'(x) = (4x³ + 3x² - 26x - 1) / 10
func evalDeriv(x float64) float64 {
	return (4*x*x*x + 3*x*x - 26*x - 1) / 10
}

// BreakCusp implements curve.ParamCurveFit.
func (q Quartic) BreakCusp(start float64, end float64) (float64, bool) {
	return 0, false
}

// SamplePtDeriv implements curve.ParamCurveFit.
func (q Quartic) SamplePtDeriv(t float64) (curve.Point, curve.Vec2) {
	// t ∈ [0, 1] but we want to plot our function from x=-4.5 to x=3.5
	x := (1-t)*-4.5 + t*3.5
	// We negate y and dy because our coordinate system is y-down
	y := -eval(x)
	dx := 1.0
	dy := -evalDeriv(x)
	return curve.Pt(x, y), curve.Vec(dx, dy)
}

// SamplePtTangent implements curve.ParamCurveFit.
func (q Quartic) SamplePtTangent(t float64, sign float64) curve.CurveFitSample {
	p, tangent := q.SamplePtDeriv(t)
	return curve.CurveFitSample{
		Point:   p,
		Tangent: tangent,
	}
}

func ExampleFitToBezPathOpt() {
	// We approximate the curve to an accuracy of 5e-3 because the final image will be
	// scaled up significantly, probably by a factor of 100. In a real use case we'd
	// probably use 1e-3, but 5e-3 keeps the example output quite a bit shorter.
	p := curve.FitToBezPathOpt(Quartic{}, 0.005)

	// We'll draw the curve we computed as an SVG document.
	fmt.Println(`<svg viewBox="-4.5 -8 8 15" xmlns="http://www.w3.org/2000/svg">`)
	// Draw a background spanning the x range we've plotted.
	fmt.Println(`<rect x="-4.5" y="-8" width="8" height="23" fill="#CCC" />`)
	// Draw an x axis at y=0.
	fmt.Println(`<path d="M-4.5,0 L3.5,0" fill="none" stroke="black" stroke-width="0.1" />`)

	colors := []string{"red", "green", "blue", "purple", "pink", "magenta", "aqua"}
	i := 0
	for seg := range curve.Segments(p.Elements()) {
		// All segments in our path cubic beziers
		c := colors[i]
		i = (i + 1) % len(colors)
		// We display every segment as a separate path with a separate color so the
		// piecewise approximation of the original curve is evident.
		fmt.Printf(`<path d="M%f,%f C%f,%f %f,%f %f,%f" stroke="%s" fill="none" stroke-width="0.1" />`,
			seg.P0.X, seg.P0.Y,
			seg.P1.X, seg.P1.Y,
			seg.P2.X, seg.P2.Y,
			seg.P3.X, seg.P3.Y,
			c,
		)
		fmt.Println()
	}
	fmt.Println("</svg>")

	// Output:
	// <svg viewBox="-4.5 -8 8 15" xmlns="http://www.w3.org/2000/svg">
	// <rect x="-4.5" y="-8" width="8" height="23" fill="#CCC" />
	// <path d="M-4.5,0 L3.5,0" fill="none" stroke="black" stroke-width="0.1" />
	// <path d="M-4.500000,-6.018750 C-4.284341,-1.969751 -4.049187,0.666184 -4.049187,0.666184" stroke="red" fill="none" stroke-width="0.1" />
	// <path d="M-4.049187,0.666184 C-4.049187,0.666184 -3.829442,3.129388 -3.585367,4.436979" stroke="green" fill="none" stroke-width="0.1" />
	// <path d="M-3.585367,4.436979 C-3.585367,4.436979 -3.411264,5.369710 -3.218875,5.747386" stroke="blue" fill="none" stroke-width="0.1" />
	// <path d="M-3.218875,5.747386 C-3.218875,5.747386 -3.093661,5.993191 -2.957352,6.011324" stroke="purple" fill="none" stroke-width="0.1" />
	// <path d="M-2.957352,6.011324 C-2.828198,6.028505 -2.684500,5.841225 -2.684500,5.841225" stroke="pink" fill="none" stroke-width="0.1" />
	// <path d="M-2.684500,5.841225 C-2.478873,5.573236 -2.210352,4.823255 -2.210352,4.823255" stroke="magenta" fill="none" stroke-width="0.1" />
	// <path d="M-2.210352,4.823255 C-2.095232,4.501726 -1.599559,2.920835 -1.599559,2.920835" stroke="aqua" fill="none" stroke-width="0.1" />
	// <path d="M-1.599559,2.920835 C-1.202513,1.654502 -0.934874,1.048021 -0.934874,1.048021" stroke="red" fill="none" stroke-width="0.1" />
	// <path d="M-0.934874,1.048021 C-0.682836,0.476892 -0.448107,0.221195 -0.448107,0.221195" stroke="green" fill="none" stroke-width="0.1" />
	// <path d="M-0.448107,0.221195 C-0.340043,0.103478 -0.232610,0.048044 -0.232610,0.048044" stroke="blue" fill="none" stroke-width="0.1" />
	// <path d="M-0.232610,0.048044 C-0.059986,-0.041027 0.116203,0.028999 0.116203,0.028999" stroke="purple" fill="none" stroke-width="0.1" />
	// <path d="M0.116203,0.028999 C0.325542,0.112201 0.554976,0.429316 0.554976,0.429316" stroke="pink" fill="none" stroke-width="0.1" />
	// <path d="M0.554976,0.429316 C0.554976,0.429316 0.703106,0.634056 0.877215,0.961362" stroke="magenta" fill="none" stroke-width="0.1" />
	// <path d="M0.877215,0.961362 C0.877215,0.961362 0.953841,1.105412 1.367244,1.961852" stroke="aqua" fill="none" stroke-width="0.1" />
	// <path d="M1.367244,1.961852 C1.367244,1.961852 1.697153,2.645317 1.895755,2.888713" stroke="red" fill="none" stroke-width="0.1" />
	// <path d="M1.895755,2.888713 C1.895755,2.888713 2.070322,3.102654 2.218996,3.105883" stroke="green" fill="none" stroke-width="0.1" />
	// <path d="M2.218996,3.105883 C2.361299,3.108974 2.489654,2.921687 2.489654,2.921687" stroke="blue" fill="none" stroke-width="0.1" />
	// <path d="M2.489654,2.921687 C2.678812,2.645679 2.848129,1.939653 2.848129,1.939653" stroke="purple" fill="none" stroke-width="0.1" />
	// <path d="M2.848129,1.939653 C3.085887,0.948240 3.298233,-0.950093 3.298233,-0.950093" stroke="pink" fill="none" stroke-width="0.1" />
	// <path d="M3.298233,-0.950093 C3.298233,-0.950093 3.401349,-1.871932 3.500000,-3.018750" stroke="magenta" fill="none" stroke-width="0.1" />
	// </svg>
}
