package curve

import (
	"fmt"
	"io"
	"iter"
	"math"
	"strconv"
	"strings"
)

// MaxExtrema is the maximum number of extrema that can be reported by
// [Extremer].
//
// This is 4 to support cubic Béziers. If other curves are used, they should be
// subdivided to limit the number of extrema.
const MaxExtrema = 4

// DefaultAccuracy is a default value for methods that take an accuracy
// argument. It is suitable for general-purpose use, such as 2D graphics.
const DefaultAccuracy = 1e-6

// Extremer describes parametrized curves that report their extrema.
type Extremer interface {
	// Extrema computes the extrema of the curve.
	//
	// Only extrema within the interior of the curve count.
	// At most four extrema can be reported, which is sufficient for
	// cubic Béziers.
	//
	// The extrema should be reported in increasing parameter order.
	Extrema() ([MaxExtrema]float64, int)
}

// ExtremaRanges returns parameter ranges, each of which is monotonic within the
// range.
func ExtremaRanges(e Extremer) ([MaxExtrema + 1][2]float64, int) {
	var ret [5][2]float64
	var retN int
	var t0 float64

	ex, n := e.Extrema()
	for _, t := range ex[:n] {
		ret[retN] = [2]float64{t0, t}
		retN++
		t0 = t
	}
	ret[retN] = [2]float64{t0, 1}
	retN++
	return ret, retN
}

// BoundingBox returns the smallest (axis-aligned) rectangle that encloses the
// curve in the range [0, 1].
func BoundingBox(c interface {
	Extremer
	ParametricCurve
}) Rect {
	bbox := NewRectFromPoints(c.Eval(0), c.Eval(1))
	ex, n := c.Extrema()
	for _, t := range ex[:n] {
		bbox = bbox.UnionPoint(c.Eval(t))
	}
	return bbox
}

type ClosedShape interface {
	Shape
	// Area returns the signed area of the closed shape.
	//
	// The convention for positive area is that y increases when x is positive.
	// Thus, it is clockwise when down is increasing y (the usual convention for
	// graphics), and anticlockwise when up is increasing y (the usual
	// convention for math).
	Area() float64

	// Winding returns the [winding number] of a point.
	//
	// The sign of the winding number is consistent with that of
	// [ClosedShape.Area], meaning it is +1 when the point is inside a positive
	// area shape and -1 when it is inside a negative area shape. Of course,
	// greater magnitude values are also possible when the shape is more
	// complex.
	//
	// [winding number]: https://en.wikipedia.org/wiki/Winding_number
	Winding(pt Point) int

	Contains(pt Point) bool
}

type Shape interface {
	// Perimeter returns the length of a shape's perimeter.
	Perimeter(accuracy float64) float64

	// BoundingBox returns the smallest rectangle that encloses the shape.
	BoundingBox() Rect

	// PathElements returns an iterator over path elements that express the
	// shape as a series of "move to", "line to", "quadratic Bézier to", "cubic
	// Bézier to", and "close path" commands.
	//
	// The tolerance parameter controls the accuracy of conversion of geometric
	// primitives to Bézier curves, as some curves such as circles cannot be
	// represented exactly but only approximated. For drawing as in UI elements,
	// a value of 0.1 is appropriate, as it is unlikely to be visible to the
	// eye. For scientific applications, a smaller value might be appropriate.
	// Note that in general the number of cubic Bézier segments scales as
	// 'tolerance ** (-1/6)'.
	PathElements(tolerance float64) iter.Seq[PathElement]

	Path(tolerance float64) BezPath
}

// ParametricCurve describes a curve parametrized by a scalar.
//
// If the result is interpreted as a point, this represents a curve. But the
// result can be interpreted as a vector as well.
type ParametricCurve interface {
	// Eval evaluates the curve at parameter t. Generally, t is in the range [0, 1].
	Eval(t float64) Point
	// Get a subsegment of the curve for the given parameter range.
	SubsegmentCurve(start, end float64) ParametricCurve
	// Subdivide into (roughly) halves.
	SubdivideCurve() (ParametricCurve, ParametricCurve)
	Start() Point
	End() Point
}

// Arclener describes a parametrized curve that can have its arc length
// measured.
type Arclener interface {
	// Arclen returns the length of the curve.
	//
	// The result is accurate to the given accuracy (subject to roundoff errors
	// for ridiculously low values). Compute time may vary with accuracy, if the
	// curve needs to be subdivided.
	Arclen(accuracy float64) float64
}

// SignedAreaer describes a parametrized curve that can have the signed area
// under it measured.
//
// For a closed path, the signed area of the path is the sum of signed areas of
// the segments. This is a variant of the "shoelace formula."
//
// This can be computed exactly for Béziers thanks to Green's theorem, and also
// for simple curves such as circular arcs. For more exotic curves, it's
// probably best to subdivide to cubics. We leave that to the caller, which is
// why we don't give an accuracy parameter here.
type SignedAreaer interface {
	SignedArea() float64
	// XXX find a name that isn't stupid
}

// expand rounds f away from zero.
func expand(f float64) float64 {
	return math.Copysign(math.Ceil(math.Abs(f)), f)
}

// Elements converts a sequence of path segments to a sequence of path elements.
func Elements(seq iter.Seq[PathSegment]) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		var currentPos option[Point]
		for seg := range seq {
			start := seg.Start()
			if !currentPos.isSet || currentPos.value != start {
				if !yield(MoveTo(start)) {
					return
				}
			}
			if !yield(seg.PathElement()) {
				return
			}
			currentPos.set(seg.End())
		}
	}
}

// Segments converts a sequence of path elements to a sequence of path segments.
func Segments(seq iter.Seq[PathElement]) iter.Seq[PathSegment] {
	return func(yield func(PathSegment) bool) {
		first := true
		var start, last Point
		for el := range seq {
			if first {
				first = false
				switch el.Kind {
				case MoveToKind:
					start = el.P0
				case LineToKind:
					start = el.P0
				case QuadToKind:
					start = el.P1
				case CubicToKind:
					start = el.P2
				case ClosePathKind:
					panic("first path element mustn't be ClosePath")
				}
				last = start
			}

			switch el.Kind {
			case MoveToKind:
				start = el.P0
				last = el.P0
			case LineToKind:
				p := last
				last = el.P0

				if !yield(Line{p, el.P0}.Seg()) {
					return
				}
			case QuadToKind:
				p := last
				last = el.P1
				if !yield(QuadBez{p, el.P0, el.P1}.Seg()) {
					return
				}
			case CubicToKind:
				p := last
				last = el.P2
				if !yield(CubicBez{p, el.P0, el.P1, el.P2}.Seg()) {
					return
				}
			case ClosePathKind:
				if last != start {
					p := last
					last = start
					if !yield(Line{p, start}.Seg()) {
						return
					}
				}
			default:
				panic(fmt.Sprintf("unhandled case %v", el.Kind))
			}
		}
	}
}

// SVGOptions specifies optional settings for [SVG] and [WriteSVG].
type SVGOptions struct {
	// The maximum precision with which to format coordinates. A value of 0
	// chooses the highest precision necessary to unambiguously represent any
	// given coordinate.
	MaxPrecision int
}

// SVG converts a sequence of path elements to a string of SVG path commands.
//
// See [WriteSVG] for a version that writes to an [io.Writer] instead of
// returning a string.
//
// The current implementation doesn't take any special care to produce a
// short string (reducing precision, using relative movement).
func SVG(seq iter.Seq[PathElement], opts SVGOptions) string {
	sb := &strings.Builder{}
	WriteSVG(sb, seq, opts)
	return sb.String()
}

// WriteSVG converts a sequence of path elements to a string of SVG path
// commands and writes it to w.
//
// See [SVG] for a version that returns a string instead.
//
// The current implementation doesn't take any special care to produce a
// short string (reducing precision, using relative movement).
func WriteSVG(w io.Writer, seq iter.Seq[PathElement], opts SVGOptions) error {
	space := []byte(" ")
	z := []byte("Z")
	var err error
	write := func(s []byte) {
		if err != nil {
			return
		}
		_, err = w.Write(s)
	}
	writef := func(s string, v ...any) {
		if err != nil {
			return
		}
		_, err = fmt.Fprintf(w, s, v...)
	}
	format := func(n float64) string {
		maxPrec := opts.MaxPrecision
		if maxPrec <= 0 {
			return strconv.FormatFloat(n, 'f', -1, 64)
		} else {
			s := strconv.FormatFloat(n, 'f', maxPrec, 64)
			return strings.TrimRight(s, "0")
		}
	}
	first := true
	for el := range seq {
		if err != nil {
			return err
		}
		if !first {
			write(space)
		}
		first = false
		switch el.Kind {
		case MoveToKind:
			writef("M%s,%s", format(el.P0.X), format(el.P0.Y))
		case LineToKind:
			writef("L%s,%s", format(el.P0.X), format(el.P0.Y))
		case QuadToKind:
			writef("Q%s,%s %s,%s",
				format(el.P0.X), format(el.P0.Y),
				format(el.P1.X), format(el.P1.Y))
		case CubicToKind:
			writef("C%s,%s %s,%s %s,%s",
				format(el.P0.X), format(el.P0.Y),
				format(el.P1.X), format(el.P1.Y),
				format(el.P2.X), format(el.P2.Y))
		case ClosePathKind:
			write(z)
		default:
			panic("unreachable")
		}
	}
	return err
}

// SolveQuadratic finds real roots of a quadratic equation.
//
// Returns values of x for which c0 + c1 x + c2 x² = 0.0
//
// This function tries to be quite numerically robust. If the equation is nearly
// linear, it will return the root ignoring the quadratic term; the other root
// might be out of representable range. In the degenerate case where all
// coefficients are zero, so that all values of x satisfy the equation, a single
// 0.0 is returned.
func SolveQuadratic(c0, c1, c2 float64) ([2]float64, int) {
	sc0 := c0 / c2
	sc1 := c1 / c2
	if math.IsInf(sc0, 0) || math.IsInf(sc1, 0) {
		// c2 is zero or very small, treat as linear eqn
		root := -c0 / c1
		if !math.IsInf(root, 0) {
			return [2]float64{root}, 1
		} else if c0 == 0.0 && c1 == 0.0 {
			// Degenerate case
			return [2]float64{0}, 1
		} else {
			return [2]float64{}, 0
		}
	}
	arg := sc1*sc1 - 4.0*sc0
	var root1 float64
	if math.IsInf(arg, 0) {
		// Likely, calculation of sc1 * sc1 overflowed. Find one root
		// using sc1 x + x² = 0, other root as sc0 / root1.
		root1 = -sc1
	} else {
		if arg < 0.0 {
			return [2]float64{}, 0
		} else if arg == 0.0 {
			return [2]float64{-0.5 * sc1}, 1
		}
		// See https://math.stackexchange.com/questions/866331
		root1 = -0.5 * (sc1 + math.Copysign(math.Sqrt(arg), sc1))
	}
	root2 := sc0 / root1
	if !math.IsInf(root2, 0) {
		// Sort just to be friendly and make results deterministic.
		if root2 > root1 {
			return [2]float64{root1, root2}, 2
		} else {
			return [2]float64{root2, root1}, 2
		}
	} else {
		return [2]float64{root1}, 1
	}
}

// SolveCubic finds real roots of cubic equations.
//
// The implementation is not (yet) fully robust, but it does handle the case
// where c3 is zero (in that case, solving the quadratic equation).
//
// See: https://momentsingraphics.de/CubicRoots.html
//
// That implementation is in turn based on Jim Blinn's "How to Solve a Cubic
// Equation", which is masterful.
//
// Returns values of x for which c0 + c1 x + c2 x² + c3 x³ = 0.0
//
// The second return value states how many roots were found.
func SolveCubic(c0, c1, c2, c3 float64) ([3]float64, int) {
	c3Recip := 1.0 / c3
	scaledC2 := c2 * (1.0 / 3.0 * c3Recip)
	scaledC1 := c1 * (1.0 / 3.0 * c3Recip)
	scaledC0 := c0 * c3Recip
	if math.IsInf(scaledC0, 0) || math.IsInf(scaledC1, 0) || math.IsInf(scaledC2, 0) {
		// cubic coefficient is zero or nearly so.
		roots, n := SolveQuadratic(c0, c1, c2)
		return [3]float64{roots[0], roots[1]}, n
	}
	c0, c1, c2 = scaledC0, scaledC1, scaledC2
	// (d0, d1, d2) is called "Delta" in article
	d0 := math.FMA(-c2, c2, c1)
	d1 := math.FMA(-c1, c2, c0)
	d2 := c2*c0 - c1*c1
	// d is called "Discriminant"
	d := 4.0*d0*d2 - d1*d1
	// de is called "Depressed.x", Depressed.y = d0
	de := math.FMA(-2.0*c2, d0, d1)
	// TODO: handle the cases where these intermediate results overflow.
	if d < 0.0 {
		sq := math.Sqrt(-0.25 * d)
		r := -0.5 * de
		t1 := math.Cbrt(r+sq) + math.Cbrt(r-sq)
		return [3]float64{t1 - c2}, 1
	} else if d == 0.0 {
		t1 := math.Copysign(math.Sqrt(-d0), de)
		return [3]float64{t1 - c2, -2.0*t1 - c2}, 2
	} else {
		th := math.Atan2(math.Sqrt(d), -de) * (1.0 / 3.0)
		// (thCos, thSin) is called "CubicRoot"
		thSin, thCos := math.Sincos(th)
		// (r0, r1, r2) is called "Root"
		r0 := thCos
		ss3 := thSin * math.Sqrt(3.0)
		r1 := 0.5 * (-thCos + ss3)
		r2 := 0.5 * (-thCos - ss3)
		t := 2.0 * math.Sqrt(-d0)

		return [3]float64{
			math.FMA(t, r0, -c2),
			math.FMA(t, r1, -c2),
			math.FMA(t, r2, -c2),
		}, 3
	}
}

// Dominant root of depressed cubic x^3 + gx + h = 0.0
//
// Section 2.2 of Orellana and De Michele.
func depressedCubicDominant(g float64, h float64) float64 {
	// Note: some of the techniques in here might be useful to improve the
	// cubic solver, and vice versa.
	q := (-1.0 / 3.0) * g
	r := 0.5 * h
	var phi0 float64
	var k option[float64]
	if math.Abs(q) < 1e102 && math.Abs(r) < 1e154 {
		k.clear()
	} else if math.Abs(q) < math.Abs(r) {
		k.set(1.0 - q*((q/r)*(q/r)))
	} else {
		v := ((r/q)*(r/q))/q - 1.0
		if math.Signbit(q) {
			v = -v
		}
		k.set(v)
	}
	if k.isSet && r == 0.0 {
		if g > 0.0 {
			phi0 = 0.0
		} else {
			phi0 = math.Sqrt(-g)
		}
	} else if k.isSet && k.value < 0.0 || !k.isSet && r*r < q*q*q {
		var t float64
		if k.isSet {
			t = r / q / math.Sqrt(q)
		} else {
			t = r / math.Sqrt(q*q*q)
		}
		phi0 = -2.0 * math.Sqrt(q) * math.Copysign(math.Cos(math.Acos(math.Abs(t))*(1.0/3.0)), t)
	} else {
		var a float64
		if k.isSet {
			if math.Abs(q) < math.Abs(r) {
				a = -r * (1.0 + math.Sqrt(k.value))
			} else {
				a = -r - math.Copysign(math.Sqrt(math.Abs(q))*q*math.Sqrt(k.value), r)
			}
		} else {
			a = -r - math.Copysign(math.Sqrt(r*r-q*q*q), r)
		}
		a = math.Cbrt(a)
		var b float64
		if a == 0.0 {
			b = 0.0
		} else {
			b = q / a
		}
		phi0 = a + b
	}
	// Refine with Newton-Raphson iteration
	x := phi0
	f := (x*x+g)*x + h
	const epsM = 2.22045e-16
	if math.Abs(f) < epsM*max(x*x*x, g*x, h) {
		return x
	}
	for range 8 {
		deltaF := 3.0*x*x + g
		if deltaF == 0.0 {
			break
		}
		newX := x - f/deltaF
		newF := (newX*newX+g)*newX + h
		if newF == 0.0 {
			return newX
		}
		if math.Abs(newF) >= math.Abs(f) {
			break
		}
		x = newX
		f = newF
	}
	return x
}

// SolveQuartic finds real roots of quartic equations.
//
// This is a fairly literal implementation of the method described in:
// Algorithm 1010: Boosting Efficiency in Solving Quartic Equations with
// No Compromise in Accuracy, Orellana and De Michele, ACM
// Transactions on Mathematical Software, Vol. 46, No. 2, May 2020.
func SolveQuartic(c0, c1, c2, c3, c4 float64) ([4]float64, int) {
	if c4 == 0.0 {
		ret, n := SolveCubic(c0, c1, c2, c3)
		return [4]float64{ret[0], ret[1], ret[2], 0}, n
	}
	if c0 == 0.0 {
		// Note: appends 0 root at end, doesn't sort. We might want to do that.
		res, n := SolveCubic(c1, c2, c3, c4)
		return [4]float64{res[0], res[1], res[2], 0}, n
	}
	a := c3 / c4
	b := c2 / c4
	c := c1 / c4
	d := c0 / c4
	if result, n, ok := solveQuarticInner(a, b, c, d, false); ok {
		return result, n
	}
	// Do polynomial rescaling
	const kq = 7.16e76
	for _, rescale := range []bool{false, true} {
		if result, n, ok := solveQuarticInner(
			a/kq,
			b/(kq*kq),
			c/(kq*kq*kq),
			d/(kq*kq*kq*kq),
			rescale,
		); ok {
			for i := range result[:n] {
				result[i] = result[i] * kq
			}
			return result, n
		}
	}
	// Overflow happened, just return no roots.
	return [4]float64{}, 0
}

func solveQuarticInner(a float64, b float64, c float64, d float64, rescale bool) ([4]float64, int, bool) {
	vs, ok := factorQuarticInner(a, b, c, d, rescale)
	if !ok {
		return [4]float64{}, 0, false
	}
	var out [4]float64
	var outN int
	for _, v := range vs {
		roots, n := SolveQuadratic(v[1], v[0], 1.0)
		for _, root := range roots[:n] {
			out[outN] = root
			outN++
		}
	}
	return out, outN, true
}

// factorQuarticInner factors a quartic into two quadratics.
//
// Attempt to factor a quartic equation into two quadratic equations. Returns
// false either if there is overflow (in which case rescaling might succeed) or
// the factorization would result in complex coefficients.
//
// Discussion question: distinguish the two cases in return value?
func factorQuarticInner(
	a float64,
	b float64,
	c float64,
	d float64,
	rescale bool,
) ([2][2]float64, bool) {
	calcEpsQ := func(a1, b1, a2, b2 float64) float64 {
		epsA := relativeEpsilon(a1+a2, a)
		epsB := relativeEpsilon(b1+a1*a2+b2, b)
		epsC := relativeEpsilon(b1*a2+a1*b2, c)
		return epsA + epsB + epsC
	}
	calcEpsT := func(a1, b1, a2, b2 float64) float64 {
		return calcEpsQ(a1, b1, a2, b2) + relativeEpsilon(b1*b2, d)
	}
	disc := 9.0*a*a - 24.0*b
	var s float64
	if disc >= 0.0 {
		s = -2.0 * b / (3.0*a + math.Copysign(math.Sqrt(disc), a))
	} else {
		s = -0.25 * a
	}
	aPrime := a + 4.0*s
	bPrime := b + 3.0*s*(a+2.0*s)
	cPrime := c + s*(2.0*b+s*(3.0*a+4.0*s))
	dPrime := d + s*(c+s*(b+s*(a+s)))
	gPrime := 0.0
	hPrime := 0.0
	const kc = 3.49e102
	if rescale {
		aPrimeS := aPrime / kc
		bPrimeS := bPrime / kc
		cPrimeS := cPrime / kc
		dPrimeS := dPrime / kc
		gPrime = aPrimeS*cPrimeS - (4.0/kc)*dPrimeS - (1.0/3.0)*bPrimeS*bPrimeS
		hPrime = (aPrimeS*cPrimeS+(8.0/kc)*dPrimeS-(2.0/9.0)*bPrimeS*bPrimeS)*
			(1.0/3.0)*
			bPrimeS -
			cPrimeS*(cPrimeS/kc) -
			aPrimeS*aPrimeS*dPrimeS
	} else {
		gPrime = aPrime*cPrime - 4.0*dPrime - (1.0/3.0)*bPrime*bPrime
		hPrime = (aPrime*cPrime+8.0*dPrime-(2.0/9.0)*bPrime*bPrime)*(1.0/3.0)*bPrime -
			cPrime*cPrime -
			aPrime*aPrime*dPrime
	}
	if math.IsInf(gPrime, 0) || math.IsInf(hPrime, 0) {
		return [2][2]float64{}, false
	}
	phi := depressedCubicDominant(gPrime, hPrime)
	if rescale {
		phi *= kc
	}
	l1 := a * 0.5
	l3 := (1.0/6.0)*b + 0.5*phi
	delt2 := c - a*l3
	d2Cand1 := (2.0/3.0)*b - phi - l1*l1
	l2Cand1 := 0.5 * delt2 / d2Cand1
	l2Cand2 := 2.0 * (d - l3*l3) / delt2
	d2Cand2 := 0.5 * delt2 / l2Cand2
	d2Cand3 := d2Cand1
	l2Cand3 := l2Cand2
	d2Best := 0.0
	l2Best := 0.0
	epsLBest := 0.0

	things := [][2]float64{
		{d2Cand1, l2Cand1},
		{d2Cand2, l2Cand2},
		{d2Cand3, l2Cand3},
	}
	for i, thing := range things {
		d2, l2 := thing[0], thing[1]
		eps0 := relativeEpsilon(d2+l1*l1+2.0*l3, b)
		eps1 := relativeEpsilon(2.0*(d2*l2+l1*l3), c)
		eps2 := relativeEpsilon(d2*l2*l2+l3*l3, d)
		epsL := eps0 + eps1 + eps2
		if i == 0 || epsL < epsLBest {
			d2Best = d2
			l2Best = l2
			epsLBest = epsL
		}
	}
	d2 := d2Best
	l2 := l2Best
	alpha1 := 0.0
	beta1 := 0.0
	alpha2 := 0.0
	beta2 := 0.0
	if d2 < 0.0 {
		sq := math.Sqrt(-d2)
		alpha1 = l1 + sq
		beta1 = l3 + sq*l2
		alpha2 = l1 - sq
		beta2 = l3 - sq*l2
		if math.Abs(beta2) < math.Abs(beta1) {
			beta2 = d / beta1
		} else if math.Abs(beta2) > math.Abs(beta1) {
			beta1 = d / beta2
		}
		var cands [][2]float64
		if math.Abs(alpha1) != math.Abs(alpha2) {
			if math.Abs(alpha1) < math.Abs(alpha2) {
				a1Cand1 := (c - beta1*alpha2) / beta2
				a1Cand2 := (b - beta2 - beta1) / alpha2
				a1Cand3 := a - alpha2
				// Note: cand 3 is first because it is infallible, simplifying logic
				cands = [][2]float64{
					{a1Cand3, alpha2},
					{a1Cand1, alpha2},
					{a1Cand2, alpha2},
				}
			} else {
				a2Cand1 := (c - alpha1*beta2) / beta1
				a2Cand2 := (b - beta2 - beta1) / alpha1
				a2Cand3 := a - alpha1
				cands = [][2]float64{
					{alpha1, a2Cand3},
					{alpha1, a2Cand1},
					{alpha1, a2Cand2},
				}
			}
			epsQBest := 0.0
			for i, c := range cands {
				a1, a2 := c[0], c[1]
				if !math.IsInf(a1, 0) && !math.IsInf(a2, 0) {
					epsQ := calcEpsQ(a1, beta1, a2, beta2)
					if i == 0 || epsQ < epsQBest {
						alpha1 = a1
						alpha2 = a2
						epsQBest = epsQ
					}
				}
			}
		}
	} else if d2 == 0.0 {
		d3 := d - l3*l3
		alpha1 = l1
		beta1 = l3 + math.Sqrt(-d3)
		alpha2 = l1
		beta2 = l3 - math.Sqrt(-d3)
		if math.Abs(beta1) > math.Abs(beta2) {
			beta2 = d / beta1
		} else if math.Abs(beta2) > math.Abs(beta1) {
			beta1 = d / beta2
		}
		// TODO: handle case d2 is very small?
	} else {
		// This case means no real roots; in the most general case we might want
		// to factor into quadratic equations with complex coefficients.
		return [2][2]float64{}, false
	}
	// Newton-Raphson iteration on alpha/beta coeff's.
	epsT := calcEpsT(alpha1, beta1, alpha2, beta2)
	for range 8 {
		if epsT == 0.0 {
			break
		}
		f0 := beta1*beta2 - d
		f1 := beta1*alpha2 + alpha1*beta2 - c
		f2 := beta1 + alpha1*alpha2 + beta2 - b
		f3 := alpha1 + alpha2 - a
		c1 := alpha1 - alpha2
		detJ := beta1*beta1 - beta1*(alpha2*c1+2.0*beta2) +
			beta2*(alpha1*c1+beta2)
		if detJ == 0.0 {
			break
		}
		inv := 1.0 / detJ
		c2 := beta2 - beta1
		c3 := beta1*alpha2 - alpha1*beta2
		dz0 := c1*f0 + c2*f1 + c3*f2 - (beta1*c2+alpha1*c3)*f3
		dz1 := (alpha1*c1+c2)*f0 -
			beta1*c1*f1 -
			beta1*c2*f2 -
			beta1*c3*f3
		dz2 := -c1*f0 - c2*f1 - c3*f2 + (alpha2*c3+beta2*c2)*f3
		dz3 := -(alpha2*c1+c2)*f0 +
			beta2*c1*f1 +
			beta2*c2*f2 +
			beta2*c3*f3
		a1 := alpha1 - inv*dz0
		b1 := beta1 - inv*dz1
		a2 := alpha2 - inv*dz2
		b2 := beta2 - inv*dz3
		newEpsT := calcEpsT(a1, b1, a2, b2)
		// We break if the new eps is equal, paper keeps going
		if newEpsT < epsT {
			alpha1 = a1
			beta1 = b1
			alpha2 = a2
			beta2 = b2
			epsT = newEpsT
		} else {
			break
		}
	}
	return [2][2]float64{{alpha1, beta1}, {alpha2, beta2}}, true
}

// SolveITP solves an arbitrary function for a zero-crossing.
//
// This uses the [ITP method], as described in the paper [An Enhancement of the
// Bisection Method Average Performance Preserving Minmax Optimality].
//
// The values of ya and yb are given as arguments rather than computed from f,
// as the values may already be known, or they may be less expensive to compute
// as special cases.
//
// It is assumed that ya < 0.0 and yb > 0.0, otherwise unexpected results may
// occur.
//
// The value of epsilon must be larger than 2**-63 * (b - a), otherwise integer
// overflow may occur. The a and b parameters represent the lower and upper
// bounds of the bracket searched for a solution.
//
// The ITP method has tuning parameters. This implementation hardwires k2 to 2,
// both because it avoids an expensive floating point exponentiation and because
// this value has been tested to work well with curve fitting problems.
//
// The n0 parameter controls the relative impact of the bisection and secant
// components. When it is 0, the number of iterations is guaranteed to be no
// more than the number required by bisection (thus, this method is strictly
// superior to bisection). However, when the function is smooth, a value of 1
// gives the secant method more of a chance to engage, so the average number of
// iterations is likely lower, though there can be one more iteration than
// bisection in the worst case.
//
// The k1 parameter is harder to characterize, and interested users are referred
// to the paper, as well as encouraged to do empirical testing. To match the
// paper, a value of 0.2 / (b - a) is suggested, and this is confirmed to give
// good results.
//
// When the function is monotonic, the returned result is guaranteed to be
// within epsilon of the zero crossing. For more detailed analysis, again see
// the paper.
//
// [ITP method]: https://en.wikipedia.org/wiki/ITP_Method
// [An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality]: https://dl.acm.org/doi/10.1145/3423597
func SolveITP(
	f func(float64) float64,
	a float64,
	b float64,
	epsilon float64,
	n0 int,
	k1 float64,
	ya float64,
	yb float64,
) float64 {
	n1_2 := int(max(math.Ceil(math.Log2((b-a)/epsilon))-1.0, 0.0))
	nmax := n0 + n1_2
	scaledEpsilon := epsilon * float64(uint64(1)<<nmax)
	for b-a > 2.0*epsilon {
		x1_2 := 0.5 * (a + b)
		r := scaledEpsilon - 0.5*(b-a)
		xf := (yb*a - ya*b) / (yb - ya)
		sigma := x1_2 - xf
		// This has k2 = 2 hardwired for efficiency.
		delta := k1 * ((b - a) * (b - a))
		var xt float64
		if delta <= math.Abs(x1_2-xf) {
			xt = xf + math.Copysign(delta, sigma)
		} else {
			xt = x1_2
		}
		var xitp float64
		if math.Abs(xt-x1_2) <= r {
			xitp = xt
		} else {
			xitp = x1_2 - math.Copysign(r, sigma)
		}
		yitp := f(xitp)
		if yitp > 0.0 {
			b = xitp
			yb = yitp
		} else if yitp < 0.0 {
			a = xitp
			ya = yitp
		} else {
			return xitp
		}
		scaledEpsilon *= 0.5
	}
	return 0.5 * (a + b)
}

type result[T, E any] struct {
	isOK bool
	ok   T
	err  E
}

// A variant ITP solver that allows fallible functions.
//
// Another difference: it returns the bracket that contains the root, which may
// be important if the function has a discontinuity.
func solveITPFallible[E any](
	f func(float64) result[float64, E],
	a float64,
	b float64,
	epsilon float64,
	n0 uint,
	k1 float64,
	ya float64,
	yb float64,
) result[[2]float64, E] {
	n1_2 := uint(max(math.Ceil(math.Log2((b-a)/epsilon))-1.0, 0.0))
	nmax := n0 + n1_2
	scaledEpsilon := epsilon * float64(uint64(1)<<nmax)
	for b-a > 2.0*epsilon {
		x1_2 := 0.5 * (a + b)
		r := scaledEpsilon - 0.5*(b-a)
		xf := (yb*a - ya*b) / (yb - ya)
		sigma := x1_2 - xf
		// This has k2 = 2 hardwired for efficiency.
		delta := k1 * ((b - a) * (b - a))
		var xt float64
		if delta <= math.Abs(x1_2-xf) {
			xt = xf + math.Copysign(delta, sigma)
		} else {
			xt = x1_2
		}
		var xitp float64
		if math.Abs(xt-x1_2) <= r {
			xitp = xt
		} else {
			xitp = x1_2 - math.Copysign(r, sigma)
		}
		yitp := f(xitp)
		if !yitp.isOK {
			return result[[2]float64, E]{err: yitp.err}
		}
		if yitp.ok > 0.0 {
			b = xitp
			yb = yitp.ok
		} else if yitp.ok < 0.0 {
			a = xitp
			ya = yitp.ok
		} else {
			return result[[2]float64, E]{isOK: true, ok: [2]float64{xitp, xitp}}
		}
		scaledEpsilon *= 0.5
	}
	return result[[2]float64, E]{isOK: true, ok: [2]float64{a, b}}
}

type option[T any] struct {
	isSet bool
	value T
}

func (opt *option[T]) set(v T) {
	opt.isSet = true
	opt.value = v
}

func (opt *option[T]) clear() {
	opt.isSet = false
	opt.value = *new(T)
}

func (opt *option[T]) unwrap() T {
	if !opt.isSet {
		panic("option isn't set")
	}
	return opt.value
}

// FactorQuarticInner factors a quartic into two quadratics.
//
// Returns false either if there is overflow (in which case rescaling might succeed) or
// the factorization would result in complex coefficients.
//
// Discussion question: distinguish the two cases in return value?
func FactorQuarticInner(
	a float64,
	b float64,
	c float64,
	d float64,
	rescale bool,
) ([2][2]float64, bool) {
	calcEpsQ := func(a1, b1, a2, b2 float64) float64 {
		epsA := relativeEpsilon(a1+a2, a)
		epsB := relativeEpsilon(b1+a1*a2+b2, b)
		epsC := relativeEpsilon(b1*a2+a1*b2, c)
		return epsA + epsB + epsC
	}
	calcEpsT := func(a1, b1, a2, b2 float64) float64 {
		return calcEpsQ(a1, b1, a2, b2) + relativeEpsilon(b1*b2, d)
	}
	disc := 9.0*a*a - 24.0*b
	var s float64
	if disc >= 0.0 {
		s = -2.0 * b / (3.0*a + math.Copysign(math.Sqrt(disc), a))
	} else {
		s = -0.25 * a
	}
	aPrime := a + 4.0*s
	bPrime := b + 3.0*s*(a+2.0*s)
	cPrime := c + s*(2.0*b+s*(3.0*a+4.0*s))
	dPrime := d + s*(c+s*(b+s*(a+s)))
	gPrime := 0.0
	hPrime := 0.0
	const kc = 3.49e102
	if rescale {
		aPrimeS := aPrime / kc
		bPrimeS := bPrime / kc
		cPrimeS := cPrime / kc
		dPrimeS := dPrime / kc
		gPrime = aPrimeS*cPrimeS - (4.0/kc)*dPrimeS - (1.0/3.0)*(bPrimeS*bPrimeS)
		hPrime = (aPrimeS*cPrimeS+(8.0/kc)*dPrimeS-(2.0/9.0)*(bPrimeS*bPrimeS))*
			(1.0/3.0)*
			bPrimeS -
			cPrimeS*(cPrimeS/kc) -
			aPrimeS*aPrimeS*dPrimeS
	} else {
		gPrime = aPrime*cPrime - 4.0*dPrime - (1.0/3.0)*(bPrime*bPrime)
		hPrime =
			(aPrime*cPrime+8.0*dPrime-(2.0/9.0)*(bPrime*bPrime))*(1.0/3.0)*bPrime -
				(cPrime * cPrime) -
				(aPrime*aPrime)*dPrime
	}
	if math.IsInf(gPrime, 0) || math.IsInf(hPrime, 0) {
		return [2][2]float64{}, false
	}
	phi := depressedCubicDominant(gPrime, hPrime)
	if rescale {
		phi *= kc
	}
	l1 := a * 0.5
	l3 := (1.0/6.0)*b + 0.5*phi
	delt2 := c - a*l3
	d2Cand1 := (2.0/3.0)*b - phi - l1*l1
	l2Cand1 := 0.5 * delt2 / d2Cand1
	l2Cand2 := 2.0 * (d - l3*l3) / delt2
	d2Cand2 := 0.5 * delt2 / l2Cand2
	d2Cand3 := d2Cand1
	l2Cand3 := l2Cand2
	d2Best := 0.0
	l2Best := 0.0
	epsLBest := 0.0
	for i, cand := range [][2]float64{{d2Cand1, l2Cand1}, {d2Cand2, l2Cand2}, {d2Cand3, l2Cand3}} {
		d2, l2 := cand[0], cand[1]
		eps0 := relativeEpsilon(d2+l1*l1+2.0*l3, b)
		eps1 := relativeEpsilon(2.0*(d2*l2+l1*l3), c)
		eps2 := relativeEpsilon(d2*l2*l2+l3*l3, d)
		epsL := eps0 + eps1 + eps2
		if i == 0 || epsL < epsLBest {
			d2Best = d2
			l2Best = l2
			epsLBest = epsL
		}
	}
	d2 := d2Best
	l2 := l2Best
	alpha1 := 0.0
	beta1 := 0.0
	alpha2 := 0.0
	beta2 := 0.0
	if d2 < 0.0 {
		sq := math.Sqrt(-d2)
		alpha1 = l1 + sq
		beta1 = l3 + sq*l2
		alpha2 = l1 - sq
		beta2 = l3 - sq*l2
		if math.Abs(beta2) < math.Abs(beta1) {
			beta2 = d / beta1
		} else if math.Abs(beta2) > math.Abs(beta1) {
			beta1 = d / beta2
		}
		var cands [][2]float64
		if math.Abs(alpha1) != math.Abs(alpha2) {
			if math.Abs(alpha1) < math.Abs(alpha2) {
				a1Cand1 := (c - beta1*alpha2) / beta2
				a1Cand2 := (b - beta2 - beta1) / alpha2
				a1Cand3 := a - alpha2
				// Note: cand 3 is first because it is infallible, simplifying logic
				cands = [][2]float64{{a1Cand3, alpha2}, {a1Cand1, alpha2}, {a1Cand2, alpha2}}
			} else {
				a2Cand1 := (c - alpha1*beta2) / beta1
				a2Cand2 := (b - beta2 - beta1) / alpha1
				a2Cand3 := a - alpha1
				cands = [][2]float64{{alpha1, a2Cand3}, {alpha1, a2Cand1}, {alpha1, a2Cand2}}
			}
			epsQBest := 0.0
			for i, cand := range cands {
				a1, a2 := cand[0], cand[1]
				if !math.IsInf(a1, 0) && !math.IsInf(a2, 0) {
					epsQ := calcEpsQ(a1, beta1, a2, beta2)
					if i == 0 || epsQ < epsQBest {
						alpha1 = a1
						alpha2 = a2
						epsQBest = epsQ
					}
				}
			}
		}
	} else if d2 == 0.0 {
		d3 := d - l3*l3
		alpha1 = l1
		beta1 = l3 + math.Sqrt(-d3)
		alpha2 = l1
		beta2 = l3 - math.Sqrt(-d3)
		if math.Abs(beta1) > math.Abs(beta2) {
			beta2 = d / beta1
		} else if math.Abs(beta2) > math.Abs(beta1) {
			beta1 = d / beta2
		}
		// TODO: handle case d2 is very small?
	} else {
		// This case means no real roots; in the most general case we might want
		// to factor into quadratic equations with complex coefficients.
		return [2][2]float64{}, false
	}
	// Newton-Raphson iteration on alpha/beta coeff's.
	epsT := calcEpsT(alpha1, beta1, alpha2, beta2)
	for range 8 {
		if epsT == 0.0 {
			break
		}
		f0 := beta1*beta2 - d
		f1 := beta1*alpha2 + alpha1*beta2 - c
		f2 := beta1 + alpha1*alpha2 + beta2 - b
		f3 := alpha1 + alpha2 - a
		c1 := alpha1 - alpha2
		detJ := beta1*beta1 - beta1*(alpha2*c1+2.0*beta2) +
			beta2*(alpha1*c1+beta2)
		if detJ == 0.0 {
			break
		}
		inv := 1.0 / detJ
		c2 := beta2 - beta1
		c3 := beta1*alpha2 - alpha1*beta2
		dz0 := c1*f0 + c2*f1 + c3*f2 - (beta1*c2+alpha1*c3)*f3
		dz1 := (alpha1*c1+c2)*f0 -
			beta1*c1*f1 -
			beta1*c2*f2 -
			beta1*c3*f3
		dz2 := -c1*f0 - c2*f1 - c3*f2 + (alpha2*c3+beta2*c2)*f3
		dz3 := -(alpha2*c1+c2)*f0 +
			beta2*c1*f1 +
			beta2*c2*f2 +
			beta2*c3*f3
		a1 := alpha1 - inv*dz0
		b1 := beta1 - inv*dz1
		a2 := alpha2 - inv*dz2
		b2 := beta2 - inv*dz3
		newEpsT := calcEpsT(a1, b1, a2, b2)
		// We break if the new eps is equal, paper keeps going
		if newEpsT < epsT {
			alpha1 = a1
			beta1 = b1
			alpha2 = a2
			beta2 = b2
			epsT = newEpsT
		} else {
			break
		}
	}
	return [2][2]float64{{alpha1, beta1}, {alpha2, beta2}}, true
}

// Compute epsilon relative to coefficient.
//
// A helper function from the Orellana and De Michele paper.
func relativeEpsilon(raw float64, a float64) float64 {
	if a == 0.0 {
		return math.Abs(raw)
	} else {
		return math.Abs((raw - a) / a)
	}
}

// ArclenSolver can be implemented by types that have a better way of computing
// the solution than the one used by [SolveForArclen].
type ArclenSolver interface {
	SolveForArclen(arclen float64, accuracy float64) float64
}

// SolveForArclen solves for the parameter that has the given arc length from
// the start of the curve.
//
// This implementation uses the [IPT method], as provided by [SolveITP]. This is
// as robust as bisection but typically converges faster. In addition, the
// method takes care to compute arc lengths of increasingly smaller segments of
// the curve, as that is likely faster than repeatedly computing the arc length
// of the segment starting at t=0.
//
// Types can optionally implement [ArclenSolver], in which case this function
// will defer to it.
//
// [ITP method]: https://en.wikipedia.org/wiki/ITP_Method
func SolveForArclen(curve interface {
	ParametricCurve
	Arclener
}, arclen float64, accuracy float64) float64 {
	if curve, ok := curve.(ArclenSolver); ok {
		return curve.SolveForArclen(arclen, accuracy)
	}

	if arclen <= 0.0 {
		return 0.0
	}
	totalArclen := curve.Arclen(accuracy)
	if arclen >= totalArclen {
		return 1.0
	}
	tLast := 0.0
	arclenLast := 0.0
	epsilon := accuracy / totalArclen
	n := 1.0 - min(math.Ceil(math.Log2(epsilon)), 0.0)
	innerAccuracy := accuracy / n
	f := func(t float64) float64 {
		var rangeStart, rangeEnd, dir float64
		if t > tLast {
			rangeStart = tLast
			rangeEnd = t
			dir = 1.0
		} else {
			rangeStart = t
			rangeEnd = tLast
			dir = -1.0
		}
		arc := curve.SubsegmentCurve(rangeStart, rangeEnd).(Arclener).Arclen(innerAccuracy)
		arclenLast += arc * dir
		tLast = t
		return arclenLast - arclen
	}
	return SolveITP(f, 0.0, 1.0, epsilon, 1, 0.2, -arclen, totalArclen-arclen)
}

// Tables of Legendre-Gauss quadrature coefficients, adapted from:
// <https://pomax.github.io/bezierinfo/legendre-gauss.html>

var gaussLegendreCoeffs8 = [...][2]float64{
	{0.3626837833783620, -0.1834346424956498},
	{0.3626837833783620, 0.1834346424956498},
	{0.3137066458778873, -0.5255324099163290},
	{0.3137066458778873, 0.5255324099163290},
	{0.2223810344533745, -0.7966664774136267},
	{0.2223810344533745, 0.7966664774136267},
	{0.1012285362903763, -0.9602898564975363},
	{0.1012285362903763, 0.9602898564975363},
}

var gaussLegendreCoeffs8Half = [...][2]float64{
	{0.3626837833783620, 0.1834346424956498},
	{0.3137066458778873, 0.5255324099163290},
	{0.2223810344533745, 0.7966664774136267},
	{0.1012285362903763, 0.9602898564975363},
}

var gaussLegendreCoeffs16 = [...][2]float64{
	{0.1894506104550685, -0.0950125098376374},
	{0.1894506104550685, 0.0950125098376374},
	{0.1826034150449236, -0.2816035507792589},
	{0.1826034150449236, 0.2816035507792589},
	{0.1691565193950025, -0.4580167776572274},
	{0.1691565193950025, 0.4580167776572274},
	{0.1495959888165767, -0.6178762444026438},
	{0.1495959888165767, 0.6178762444026438},
	{0.1246289712555339, -0.7554044083550030},
	{0.1246289712555339, 0.7554044083550030},
	{0.0951585116824928, -0.8656312023878318},
	{0.0951585116824928, 0.8656312023878318},
	{0.0622535239386479, -0.9445750230732326},
	{0.0622535239386479, 0.9445750230732326},
	{0.0271524594117541, -0.9894009349916499},
	{0.0271524594117541, 0.9894009349916499},
}

var gaussLegendreCoeffs16Half = [...][2]float64{
	{0.1894506104550685, 0.0950125098376374},
	{0.1826034150449236, 0.2816035507792589},
	{0.1691565193950025, 0.4580167776572274},
	{0.1495959888165767, 0.6178762444026438},
	{0.1246289712555339, 0.7554044083550030},
	{0.0951585116824928, 0.8656312023878318},
	{0.0622535239386479, 0.9445750230732326},
	{0.0271524594117541, 0.9894009349916499},
}

var gaussLegendreCoeffs24Half = [...][2]float64{
	{0.1279381953467522, 0.0640568928626056},
	{0.1258374563468283, 0.1911188674736163},
	{0.1216704729278034, 0.3150426796961634},
	{0.1155056680537256, 0.4337935076260451},
	{0.1074442701159656, 0.5454214713888396},
	{0.0976186521041139, 0.6480936519369755},
	{0.0861901615319533, 0.7401241915785544},
	{0.0733464814110803, 0.8200019859739029},
	{0.0592985849154368, 0.8864155270044011},
	{0.0442774388174198, 0.9382745520027328},
	{0.0285313886289337, 0.9747285559713095},
	{0.0123412297999872, 0.9951872199970213},
}
