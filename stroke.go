package curve

import (
	"iter"
	"math"
	"slices"
)

// Join defines the connection between two segments of a stroke.
type Join int

const (
	// A straight line connecting the segments.
	BevelJoin Join = iota
	// The segments are extended to their natural intersection point.
	MiterJoin
	// An arc between the segments.
	RoundJoin
)

// Cap defines the shape to be drawn at the ends of a stroke.
type Cap int

const (
	// Flat cap.
	ButtCap Cap = iota
	// Square cap with dimensions equal to half the stroke width.
	SquareCap
	// Rounded cap with radius equal to half the stroke width.
	RoundCap
)

// Stroke describes the visual style of a stroke.
type Stroke struct {
	// Width of the stroke.
	Width float64
	// Style for connecting segments of the stroke.
	Join Join
	// Limit for miter joins.
	MiterLimit float64
	// Style for capping the beginning of an open subpath.
	StartCap Cap
	// Style for capping the end of an open subpath.
	EndCap Cap
	// Lengths of dashes in alternating on/off order.
	DashPattern []float64
	// Offset of the first dash.
	DashOffset float64
}

// StrokeOpts describes options for path stroking.
type StrokeOpts struct {
	OptLevel OptLevel
}

// OptLevel defines the optimization level for computing strokes.
type OptLevel int

const (
	// Adaptively subdivide segments in half.
	Subdivide OptLevel = iota
	// Compute optimized subdivision points to minimize error.
	Optimized
)

var DefaultStroke = Stroke{
	Width:      1.0,
	Join:       RoundJoin,
	MiterLimit: 4.0,
	StartCap:   RoundCap,
	EndCap:     RoundCap,
}

func (s Stroke) WithWidth(width float64) Stroke      { s.Width = width; return s }
func (s Stroke) WithJoin(join Join) Stroke           { s.Join = join; return s }
func (s Stroke) WithMiterLimit(limit float64) Stroke { s.MiterLimit = limit; return s }
func (s Stroke) WithStartCap(cap Cap) Stroke         { s.StartCap = cap; return s }
func (s Stroke) WithEndCap(cap Cap) Stroke           { s.EndCap = cap; return s }
func (s Stroke) WithCaps(cap Cap) Stroke             { s.StartCap, s.EndCap = cap, cap; return s }
func (s Stroke) WithDashes(offset float64, pattern []float64) Stroke {
	s.DashOffset, s.DashPattern = offset, pattern
	return s
}

type strokeCtx struct {
	yield func(PathElement) bool
	dead  bool
	// forward paths are yielded directly instead of being accumulated like
	// backwardPath is. emittedForward is set to true if any elements have been
	// yielded, and reset to false after a path is completed.
	emittedForward bool
	backwardPath   BezPath
	startPt        Point
	startNorm      Vec2
	startTan       Vec2
	lastPt         Point
	lastTan        Vec2
	// Precomputation of the join threshold to optimize per-join logic.
	// If hypot < (hypot + dot) * joinThresh omit join altogether.
	joinThresh float64
}

// StrokePath expands a stroke into a fill.
//
// The tolerance parameter controls the accuracy of the result. In general,
// the number of subdivisions in the output scales to the -1/6 power of the
// parameter, for example making it 1/64 as big generates twice as many
// segments. The appropriate value depends on the application; if the result
// of the stroke will be scaled up, a smaller value is needed.
//
// This method attempts a fairly high degree of correctness, but ultimately
// is based on computing parallel curves and adding joins and caps, rather than
// computing the rigorously correct parallel sweep (which requires evolutes in
// the general case). See [Nehab 2020] for more discussion.
//
// [Nehab 2020]: https://dl.acm.org/doi/10.1145/3386569.3392392
func StrokePath(
	path iter.Seq[PathElement],
	style Stroke,
	opts StrokeOpts,
	tolerance float64,
) iter.Seq[PathElement] {
	if len(style.DashPattern) == 0 {
		return strokeUndashed(path, style, tolerance, opts)
	} else {
		dashed := Dash(path, style.DashOffset, style.DashPattern)
		return strokeUndashed(dashed, style, tolerance, opts)
	}
}

// Version of stroke expansion for styles with no dashes.
func strokeUndashed(
	path iter.Seq[PathElement],
	style Stroke,
	tolerance float64,
	opts StrokeOpts,
) iter.Seq[PathElement] {
	return func(yield func(v PathElement) bool) {
		ctx := strokeCtx{
			yield:      yield,
			joinThresh: 2.0 * tolerance / style.Width,
		}
		for el := range path {
			p0 := ctx.lastPt
			switch el.Kind {
			case MoveToKind:
				p := el.P0
				ctx.finish(style)
				ctx.startPt = p
				ctx.lastPt = p
			case LineToKind:
				p1 := el.P0
				if p1 != p0 {
					tangent := p1.Sub(p0)
					ctx.doJoin(style, tangent)
					ctx.lastTan = tangent
					ctx.doLine(style, tangent, p1)
				}
			case QuadToKind:
				p1, p2 := el.P0, el.P1
				if p1 != p0 || p2 != p0 {
					q := QuadBez{p0, p1, p2}
					tan0, tan1 := q.Tangents()
					ctx.doJoin(style, tan0)
					ctx.doCubic(style, q.Raise(), tolerance, opts)
					ctx.lastTan = tan1
				}
			case CubicToKind:
				p1, p2, p3 := el.P0, el.P1, el.P2
				if p1 != p0 || p2 != p0 || p3 != p0 {
					c := CubicBez{p0, p1, p2, p3}
					tan0, tan1 := c.Tangents()
					ctx.doJoin(style, tan0)
					ctx.doCubic(style, c, tolerance, opts)
					ctx.lastTan = tan1
				}
			case ClosePathKind:
				if p0 != ctx.startPt {
					tangent := ctx.startPt.Sub(p0)
					ctx.doJoin(style, tangent)
					ctx.lastTan = tangent
					ctx.doLine(style, tangent, ctx.startPt)
				}
				ctx.finishClosed(style)
			}
		}
		ctx.finish(style)
	}
}

func (ctx *strokeCtx) doYield(el PathElement) {
	if ctx.dead {
		return
	}
	ctx.dead = !ctx.yield(el)
	ctx.emittedForward = true
}

// Append backward path to output.
func (ctx *strokeCtx) finish(style Stroke) {
	// TODO: scale
	const tolerance = 1e-3
	if !ctx.emittedForward {
		return
	}
	returnPt, ok := ctx.backwardPath[len(ctx.backwardPath)-1].EndPoint()
	if !ok {
		panic("unreachable")
	}
	d := ctx.lastPt.Sub(returnPt)
	switch style.EndCap {
	case ButtCap:
		ctx.LineTo(returnPt)
	case RoundCap:
		roundCap(ctx, tolerance, ctx.lastPt, d)
	case SquareCap:
		squareCap(ctx, false, ctx.lastPt, d)
	}
	extendReversed(ctx, ctx.backwardPath)
	switch style.StartCap {
	case ButtCap:
		ctx.ClosePath()
	case RoundCap:
		roundCap(ctx, tolerance, ctx.startPt, ctx.startNorm)
	case SquareCap:
		squareCap(ctx, true, ctx.startPt, ctx.startNorm)
	}

	ctx.emittedForward = false
	ctx.backwardPath.Truncate(0)
}

// Finish a closed path
func (ctx *strokeCtx) finishClosed(style Stroke) {
	if !ctx.emittedForward {
		return
	}
	ctx.doJoin(style, ctx.startTan)
	ctx.ClosePath()
	lastPt, ok := ctx.backwardPath[len(ctx.backwardPath)-1].EndPoint()
	if !ok {
		panic("unreachable")
	}
	ctx.MoveTo(lastPt)
	extendReversed(ctx, ctx.backwardPath)
	ctx.ClosePath()
	ctx.emittedForward = false
	ctx.backwardPath.Truncate(0)
}

func (ctx *strokeCtx) doJoin(style Stroke, tan0 Vec2) {
	// TODO: scale
	const tolerance = 1e-3
	scale := 0.5 * style.Width / tan0.Hypot()
	norm := Vec(-tan0.Y, tan0.X).Mul(scale)
	p0 := ctx.lastPt
	if !ctx.emittedForward {
		ctx.MoveTo(p0.Translate(norm.Negate()))
		ctx.backwardPath.MoveTo(p0.Translate(norm))
		ctx.startTan = tan0
		ctx.startNorm = norm
	} else {
		ab := ctx.lastTan
		cd := tan0
		cross := ab.Cross(cd)
		dot := ab.Dot(cd)
		hypot := math.Hypot(cross, dot)
		// possible TODO: a minor speedup could be squaring both sides
		if dot <= 0.0 || math.Abs(cross) >= hypot*ctx.joinThresh {
			switch style.Join {
			case BevelJoin:
				ctx.LineTo(p0.Translate(norm.Negate()))
				ctx.backwardPath.LineTo(p0.Translate(norm))
			case MiterJoin:
				if 2.0*hypot < (hypot+dot)*style.MiterLimit*style.MiterLimit {
					// TODO: maybe better to store lastNorm or derive from path?
					lastScale := 0.5 * style.Width / ab.Hypot()
					lastNorm := Vec(-ab.Y, ab.X).Mul(lastScale)
					if cross > 0.0 {
						fpLast := p0.Translate(lastNorm.Negate())
						fpThis := p0.Translate(norm.Negate())
						h := ab.Cross(fpThis.Sub(fpLast)) / cross
						miterPt := fpThis.Translate(cd.Mul(h).Negate())
						ctx.LineTo(miterPt)
					} else if cross < 0.0 {
						fpLast := p0.Translate(lastNorm)
						fpThis := p0.Translate(norm)
						h := ab.Cross(fpThis.Sub(fpLast)) / cross
						miterPt := fpThis.Translate(cd.Mul(h).Negate())
						ctx.backwardPath.LineTo(miterPt)
					}
				}
				ctx.LineTo(p0.Translate(norm.Negate()))
				ctx.backwardPath.LineTo(p0.Translate(norm))
			case RoundJoin:
				angle := math.Atan2(cross, dot)
				if angle > 0.0 {
					ctx.backwardPath.LineTo(p0.Translate(norm))
					roundJoin(ctx, tolerance, p0, norm, angle)
				} else {
					ctx.LineTo(p0.Translate(norm.Negate()))
					roundJoinRev(&ctx.backwardPath, tolerance, p0, norm.Negate(), -angle)
				}
			}
		}
	}
}

func (ctx *strokeCtx) doLine(style Stroke, tangent Vec2, p1 Point) {
	scale := 0.5 * style.Width / tangent.Hypot()
	norm := Vec(-tangent.Y, tangent.X).Mul(scale)
	ctx.LineTo(p1.Translate(norm.Negate()))
	ctx.backwardPath.LineTo(p1.Translate(norm))
	ctx.lastPt = p1
}

func (ctx *strokeCtx) doCubic(style Stroke, c CubicBez, tolerance float64, opts StrokeOpts) {
	// First, detect degenerate linear case

	// Ordinarily, this is the direction of the chord, but if the chord is very
	// short, we take the longer control arm.
	chord := c.P3.Sub(c.P0)
	chordRef := chord
	chordRefHypot2 := chordRef.Hypot2()
	d01 := c.P1.Sub(c.P0)
	if d01.Hypot2() > chordRefHypot2 {
		chordRef = d01
		chordRefHypot2 = chordRef.Hypot2()
	}
	d23 := c.P3.Sub(c.P2)
	if d23.Hypot2() > chordRefHypot2 {
		chordRef = d23
		chordRefHypot2 = chordRef.Hypot2()
	}
	// Project Bézier onto chord
	p0 := Vec2(c.P0).Dot(chordRef)
	p1 := Vec2(c.P1).Dot(chordRef)
	p2 := Vec2(c.P2).Dot(chordRef)
	p3 := Vec2(c.P3).Dot(chordRef)
	const endpointD = 0.01
	if p3 <= p0 ||
		p1 > p2 ||
		p1 < p0+endpointD*(p3-p0) ||
		p2 > p3-endpointD*(p3-p0) {
		// potentially a cusp inside
		x01 := d01.Cross(chordRef)
		x23 := d23.Cross(chordRef)
		x03 := chord.Cross(chordRef)
		thresh := tolerance * tolerance * chordRefHypot2
		if x01*x01 < thresh && x23*x23 < thresh && x03*x03 < thresh {
			// control points are nearly co-linear
			midpoint := c.P0.Midpoint(c.P3)
			// Mapping back from projection of reference chord
			refVec := chordRef.Mul(1.0 / chordRefHypot2)
			refPt := midpoint.Translate(refVec.Mul(0.5 * (p0 + p3)).Negate())
			ctx.doLinear(style, c, [4]float64{p0, p1, p2, p3}, refPt, refVec)
			return
		}
	}

	// A tuning parameter for regularization. A value too large may distort the curve,
	// while a value too small may fail to generate smooth curves. This is a somewhat
	// arbitrary value, and should be revisited.
	const dimTune = 0.25
	dimension := tolerance * dimTune
	{
		co := NewCubicOffset(c, -0.5*style.Width, dimension)
		forward := fitWithOpts(&co, tolerance, opts)
		for el := range dropFirst(forward) {
			if ctx.dead {
				break
			}
			ctx.doYield(el)
		}
	}
	{
		co := NewCubicOffset(c, 0.5*style.Width, dimension)
		backward := fitWithOpts(&co, tolerance, opts)
		ctx.backwardPath = slices.AppendSeq(ctx.backwardPath, dropFirst(backward))
	}
	ctx.lastPt = c.P3
}

// Do a cubic which is actually linear.
//
// The p argument is the control points projected to the reference chord.
// The ref arguments are the inverse map of a projection back to the client
// coordinate space.
func (ctx *strokeCtx) doLinear(
	style Stroke,
	c CubicBez,
	p [4]float64,
	refPt Point,
	refVec Vec2,
) {
	// Always do round join, to model cusp as limit of finite curvature (see Nehab).
	style = DefaultStroke.WithWidth(style.Width).WithJoin(RoundJoin)
	// Tangents of endpoints (for connecting to joins)
	tan0, tan1 := c.Tangents()
	ctx.lastTan = tan0
	// find cusps
	c0 := p[1] - p[0]
	c1 := 2.0*p[2] - 4.0*p[1] + 2.0*p[0]
	c2 := p[3] - 3.0*p[2] + 3.0*p[1] - p[0]
	roots, n := SolveQuadratic(c0, c1, c2)
	// discard cusps right at endpoints
	const epsilon = 1e-6
	for _, t := range roots[:n] {
		if t > epsilon && t < 1.0-epsilon {
			mt := 1.0 - t
			z := mt*(mt*mt*p[0]+3.0*t*(mt*p[1]+t*p[2])) + t*t*t*p[3]
			p := refPt.Translate(refVec.Mul(z))
			tan := p.Sub(ctx.lastPt)
			ctx.doJoin(style, tan)
			ctx.doLine(style, tan, p)
			ctx.lastTan = tan
		}
	}
	tan := c.P3.Sub(ctx.lastPt)
	ctx.doJoin(style, tan)
	ctx.doLine(style, tan, c.P3)
	ctx.lastTan = tan
	ctx.doJoin(style, tan1)
}

func (ctx *strokeCtx) CubicTo(p1, p2, p3 Point) {
	ctx.doYield(CubicTo(p1, p2, p3))
}

func (ctx *strokeCtx) LineTo(p1 Point) {
	ctx.doYield(LineTo(p1))
}

func (ctx *strokeCtx) ClosePath() {
	ctx.doYield(ClosePath())
}

func (ctx *strokeCtx) QuadTo(p0, p1 Point) {
	ctx.doYield(QuadTo(p0, p1))
}

func (ctx *strokeCtx) MoveTo(pt Point) {
	ctx.doYield(MoveTo(pt))
}

func roundCap(out *strokeCtx, tolerance float64, center Point, norm Vec2) {
	roundJoin(out, tolerance, center, norm, math.Pi)
}

func roundJoin(out *strokeCtx, tolerance float64, center Point, norm Vec2, angle float64) {
	a := Affine{norm.X, norm.Y, -norm.Y, norm.X, center.X, center.Y}
	arc := Arc{Point{}, Vec(1.0, 1.0), math.Pi - angle, angle, 0.0}
	for el := range dropFirst(arc.PathElements(tolerance)) {
		out.CubicTo(
			el.P0.Transform(a),
			el.P1.Transform(a),
			el.P2.Transform(a),
		)
	}
}

func roundJoinRev(out *BezPath, tolerance float64, center Point, norm Vec2, angle float64) {
	a := Affine{norm.X, norm.Y, norm.Y, -norm.X, center.X, center.Y}
	arc := Arc{Point{}, Vec(1.0, 1.0), math.Pi - angle, angle, 0.0}
	for el := range arc.PathElements(tolerance) {
		out.CubicTo(
			el.P0.Transform(a),
			el.P1.Transform(a),
			el.P2.Transform(a),
		)
	}
}

func squareCap(out *strokeCtx, close bool, center Point, norm Vec2) {
	a := Affine{norm.X, norm.Y, -norm.Y, norm.X, center.X, center.Y}
	out.LineTo(Pt(1.0, 1.0).Transform(a))
	out.LineTo(Pt(-1.0, 1.0).Transform(a))
	if close {
		out.ClosePath()
	} else {
		out.LineTo(Pt(-1.0, 0.0).Transform(a))
	}
}

func extendReversed(out *strokeCtx, elements []PathElement) {
	for i := len(elements) - 1; i >= 1; i-- {
		end, ok := elements[i-1].EndPoint()
		if !ok {
			panic("unreachable")
		}
		el := elements[i]
		switch el.Kind {
		case MoveToKind:
			panic("unexpected MoveTo")
		case LineToKind:
			out.LineTo(end)
		case QuadToKind:
			out.QuadTo(el.P0, end)
		case CubicToKind:
			out.CubicTo(el.P1, el.P0, end)
		default:
			panic("unreachable")
		}
	}
}

func fitWithOpts(co *CubicOffset, tolerance float64, opts StrokeOpts) iter.Seq[PathElement] {
	switch opts.OptLevel {
	case Subdivide:
		return FitToBezPath(co, tolerance)
	case Optimized:
		return func(yield func(PathElement) bool) {
			for _, el := range FitToBezPathOpt(co, tolerance) {
				if !yield(el) {
					break
				}
			}
		}
	default:
		return nil
	}
}

const dashAccuracy = 1e-6

type dashState int

const (
	dashStateNeedInput dashState = iota
	dashStateToStash
	dashStateWorking
	dashStateFromStash
)

// Dash returns a dashing iterator. It consumes a sequence of path elements and produces a
// new sequence of path elements representing a dashed version of the original sequence.
//
// Accuracy is currently hard-coded to 1e-6. This is better than generally expected, and
// care is taken to get cusps correct, among other things.
func Dash(
	inner iter.Seq[PathElement],
	dashOffset float64,
	dashes []float64,
) iter.Seq[PathElement] {
	return func(yield func(PathElement) bool) {
		dashIdx := 0
		isActive := true
		dashRemaining := dashes[dashIdx] - dashOffset
		// Find place in dashes array for initial offset.
		for dashRemaining < 0.0 {
			dashIdx = (dashIdx + 1) % len(dashes)
			dashRemaining += dashes[dashIdx]
			isActive = !isActive
		}
		next, stop := iter.Pull(inner)
		diter := &dashIterator{
			innerNext:         next,
			innerStop:         stop,
			inputDone:         false,
			closepathPending:  false,
			dashes:            dashes,
			dashIdx:           dashIdx,
			initDashIdx:       dashIdx,
			initDashRemaining: dashRemaining,
			initIsActive:      isActive,
			isActive:          isActive,
			state:             dashStateNeedInput,
			currentSeg: PathSegment{
				Kind: LineKind,
			},
			t:             0.0,
			dashRemaining: dashRemaining,
			segRemaining:  0.0,
			stashIdx:      0,
		}
		defer diter.innerStop()
		for {
			switch diter.state {
			case dashStateNeedInput:
				if diter.inputDone {
					return
				}
				diter.getInput()
				if diter.inputDone {
					return
				}
				diter.state = dashStateToStash
			case dashStateToStash:
				if el, ok := diter.step(); ok {
					diter.stash = append(diter.stash, el)
				}
			case dashStateWorking:
				if el, ok := diter.step(); ok {
					if !yield(el) {
						return
					}
				}
			case dashStateFromStash:
				if diter.stashIdx < len(diter.stash) {
					el := diter.stash[diter.stashIdx]
					diter.stashIdx++
					if !yield(el) {
						return
					}
				} else {
					diter.stash = diter.stash[:0]
					diter.stashIdx = 0
					if diter.inputDone {
						return
					}
					if diter.closepathPending {
						diter.closepathPending = false
						diter.state = dashStateNeedInput
					} else {
						diter.state = dashStateToStash
					}
				}
			}
		}
	}
}

type dashIterator struct {
	innerNext func() (PathElement, bool)
	innerStop func()

	inputDone         bool
	closepathPending  bool
	dashes            []float64
	dashIdx           int
	initDashIdx       int
	initDashRemaining float64
	initIsActive      bool
	isActive          bool
	state             dashState
	currentSeg        PathSegment
	t                 float64
	dashRemaining     float64
	segRemaining      float64
	startPt           Point
	lastPt            Point
	stash             []PathElement
	stashIdx          int
}

func (di *dashIterator) getInput() {
	for {
		if di.closepathPending {
			di.handleClosepath()
			break
		}
		nextEl, ok := di.innerNext()
		if !ok {
			di.inputDone = true
			di.state = dashStateFromStash
			return
		}
		p0 := di.lastPt
		switch nextEl.Kind {
		case MoveToKind:
			if len(di.stash) != 0 {
				di.state = dashStateFromStash
			}
			di.startPt = nextEl.P0
			di.lastPt = nextEl.P0
			di.resetPhase()
			continue
		case LineToKind:
			l := Line{p0, nextEl.P0}
			di.segRemaining = l.Arclen(dashAccuracy)
			di.currentSeg = l.Seg()
			di.lastPt = nextEl.P0
		case QuadToKind:
			q := QuadBez{p0, nextEl.P0, nextEl.P1}
			di.segRemaining = q.Arclen(dashAccuracy)
			di.currentSeg = q.Seg()
			di.lastPt = nextEl.P1
		case CubicToKind:
			c := CubicBez{p0, nextEl.P0, nextEl.P1, nextEl.P2}
			di.segRemaining = c.Arclen(dashAccuracy)
			di.currentSeg = c.Seg()
			di.lastPt = nextEl.P2
		case ClosePathKind:
			di.closepathPending = true
			if p0 != di.startPt {
				l := Line{p0, di.startPt}
				di.segRemaining = l.Arclen(dashAccuracy)
				di.currentSeg = l.Seg()
				di.lastPt = di.startPt
			} else {
				di.handleClosepath()
			}
		}
		break
	}
	di.t = 0.0
}

// Move arc length forward to next event.
func (di *dashIterator) step() (PathElement, bool) {
	var result PathElement
	var hasResult bool
	if di.state == dashStateToStash && len(di.stash) == 0 {
		if di.isActive {
			result = MoveTo(di.currentSeg.Eval(0))
			hasResult = true
		} else {
			di.state = dashStateWorking
		}
	} else if di.dashRemaining < di.segRemaining {
		// next transition is a dash transition
		seg := di.currentSeg.Subsegment(di.t, 1.0)
		t1 := SolveForArclen(seg, di.dashRemaining, dashAccuracy)
		if di.isActive {
			subseg := seg.Subsegment(0.0, t1)
			result = subseg.PathElement()
			hasResult = true
			di.state = dashStateWorking
		} else {
			p := seg.Eval(t1)
			result = MoveTo(p)
			hasResult = true
		}
		di.isActive = !di.isActive
		di.t += t1 * (1.0 - di.t)
		di.segRemaining -= di.dashRemaining
		di.dashIdx++
		if di.dashIdx == len(di.dashes) {
			di.dashIdx = 0
		}
		di.dashRemaining = di.dashes[di.dashIdx]
	} else {
		if di.isActive {
			seg := di.currentSeg.Subsegment(di.t, 1.0)
			result = seg.PathElement()
			hasResult = true
		}
		di.dashRemaining -= di.segRemaining
		di.getInput()
	}
	return result, hasResult
}

func (di *dashIterator) handleClosepath() {
	if di.state == dashStateToStash {
		// Have looped back without breaking a dash, just play it back
		di.stash = append(di.stash, ClosePath())
	} else if di.isActive {
		// connect with path in stash, skip MoveTo.
		di.stashIdx = 1
	}
	di.state = dashStateFromStash
	di.resetPhase()
}

func (di *dashIterator) resetPhase() {
	di.dashIdx = di.initDashIdx
	di.dashRemaining = di.initDashRemaining
	di.isActive = di.initIsActive
}
