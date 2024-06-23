package curve

import (
	"iter"
	"math"
)

// SplitArclen subdivides each subpath into segments of arc length l. Each group
// of segments is delimited by a zero value path segment.
func SplitArclen(inner iter.Seq[PathSegment], l float64) iter.Seq[PathSegment] {
	return splitArclenMaxGroups(inner, l, -1)
}

// SplitN subdivides each subpath into n segments of identical arc length. Each
// group of segments is delimited by a zero value path segment.
//
// inner must not be a single-use iterator.
func SplitN(inner iter.Seq[PathSegment], n int) iter.Seq[PathSegment] {
	if n <= 1 {
		return inner
	}

	const splitAccuracy = 1e-6
	totalLength := 0.0
	for seg := range inner {
		totalLength += seg.Arclen(splitAccuracy)
	}

	splitLength := totalLength / float64(n)
	if math.IsInf(splitLength, 0) || math.IsNaN(splitLength) {
		return inner
	}

	return splitArclenMaxGroups(inner, splitLength, n)
}

func splitArclenMaxGroups(inner iter.Seq[PathSegment], l float64, n int) iter.Seq[PathSegment] {
	// TODO(dh): should we support splitting the remainder across the start and
	// end of the curve? Currently any remainder is at the end.
	if l == 0 {
		return inner
	}

	const splitAccuracy = 1e-6
	splitLength := l
	if math.IsInf(splitLength, 0) || math.IsNaN(splitLength) {
		return inner
	}

	// We cannot use SplitArcLength(inner, splitLength) because due to rounding
	// errors we might end up with n+1 groups.
	return func(yield func(PathSegment) bool) {
		remainingLength := splitLength
		remainingSegs := n
		var prevEnd option[Point]
		var out []PathSegment
		for seg := range inner {
			if prevEnd.isSet && seg.P0 != prevEnd.value {
				// New subpath
				remainingLength = splitLength
				remainingSegs = n
				if !yield(PathSegment{}) {
					return
				}
			}
			switch seg.Kind {
			case LineKind:
				prevEnd.set(seg.P1)
			case QuadKind:
				prevEnd.set(seg.P2)
			case CubicKind:
				prevEnd.set(seg.P3)
			}

			for {
				if a := seg.Arclen(splitAccuracy); a < remainingLength {
					remainingLength -= a
					if !yield(seg) {
						return
					}
					break
				} else {
					var t float64
					if remainingSegs == 1 {
						t = 1
					} else {
						t = SolveForArclen(seg, remainingLength, splitAccuracy)
					}
					out = append(out, seg.Subsegment(0, t))
					if !yield(seg.Subsegment(0, t)) {
						return
					}
					if remainingSegs > 0 {
						remainingSegs--
					}
					if !yield(PathSegment{}) {
						return
					}
					remainingLength = splitLength
					if t >= 1 {
						break
					} else {
						seg = seg.Subsegment(t, 1)
					}
				}
			}
		}
	}
}
