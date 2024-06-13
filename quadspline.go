package curve

import "iter"

// QuadBSpline is a quadratic B-spline. It is encoded as [P₁, C₁, C₂, C₃, C₄, ..., Pₙ],
// where Pᵢ are on-curve points and Cᵢ are off-curve control points. Only the first and
// last on-curve points are explicit. All other on-curve points are implicit and defined
// as Pᵢ = (Cᵢ₋₁ + Cᵢ) / 2. For example, P₂ lies halfway between C₁ and C₂. This format
// matches the one used by glyf tables in TrueType fonts.
type QuadBSpline []Point

// Quads returns an iterator over the implied sequence of quadratic Bézier segments. The
// returned segments are guaranteed to be G1 continuous.
func (q QuadBSpline) Quads() iter.Seq[QuadBez] {
	return func(yield func(QuadBez) bool) {
		var idx int
		for len(q[idx:]) >= 3 {
			p0, p1, p2 := q[idx], q[idx+1], q[idx+2]

			if idx != 0 {
				p0 = p0.Midpoint(p1)
			}
			if idx+2 < len(q)-1 {
				p2 = p1.Midpoint(p2)
			}

			idx++

			if !yield(QuadBez{p0, p1, p2}) {
				break
			}
		}
	}
}
