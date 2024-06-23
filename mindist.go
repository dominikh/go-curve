package curve

import "math"

//! Minimum distance between two Bézier curves
//!
//! This implements the algorithm in "Computing the minimum distance between
//! two Bézier curves", Chen et al., *Journal of Computational and Applied
//! Mathematics* 229(2009), 294-301

func minDistParam(
	bez1, bez2 []Vec2,
	u, v [2]float64,
	epsilon float64,
	bestAlpha float64,
) [3]float64 {
	if len(bez1) == 0 || len(bez2) == 0 {
		panic("called with empty path")
	}

	n := len(bez1) - 1
	m := len(bez2) - 1
	umin, umax := u[0], u[1]
	vmin, vmax := v[0], v[1]
	umid := (umin + umax) / 2.0
	vmid := (vmin + vmax) / 2.0
	svalues := [4][3]float64{
		{s(umin, vmin, bez1, bez2), umin, vmin},
		{s(umin, vmax, bez1, bez2), umin, vmax},
		{s(umax, vmin, bez1, bez2), umax, vmin},
		{s(umax, vmax, bez1, bez2), umax, vmax},
	}
	alpha := svalues[0][0]
	for _, sval := range svalues {
		alpha = min(alpha, sval[0])
	}
	if alpha > bestAlpha {
		return [3]float64{alpha, umid, vmid}
	}

	if math.Abs(umax-umin) < epsilon || math.Abs(vmax-vmin) < epsilon {
		return [3]float64{alpha, umid, vmid}
	}

	// Property one: D(r>k) > alpha
	isOutside := true
	var minDrk option[float64]
	var minIj option[[2]int]
	for r := range 2 * n {
		for k := range 2 * m {
			dRk := dRk(r, k, bez1, bez2)
			if dRk < alpha {
				isOutside = false
			}
			if !minDrk.isSet || dRk < minDrk.value {
				minDrk.set(dRk)
				minIj.set([2]int{r, k})
			}
		}
	}
	if isOutside {
		return [3]float64{alpha, umid, vmid}
	}

	// Property two: boundary check
	atBoundary0OnBez1 := true
	atBoundary1OnBez1 := true
	atBoundary0OnBez2 := true
	atBoundary1OnBez2 := true
	for i := range 2 * n {
		for j := range 2 * m {
			dij := dRk(i, j, bez1, bez2)
			dkj := dRk(0, j, bez1, bez2)
			if dij < dkj {
				atBoundary0OnBez1 = false
			}
			dkj = dRk(2*n, j, bez1, bez2)
			if dij < dkj {
				atBoundary1OnBez1 = false
			}
			dkj = dRk(i, 0, bez1, bez2)
			if dij < dkj {
				atBoundary0OnBez2 = false
			}
			dkj = dRk(i, 2*n, bez1, bez2)
			if dij < dkj {
				atBoundary1OnBez2 = false
			}
		}
	}
	if atBoundary0OnBez1 && atBoundary0OnBez2 {
		return svalues[0]
	}
	if atBoundary0OnBez1 && atBoundary1OnBez2 {
		return svalues[1]
	}
	if atBoundary1OnBez1 && atBoundary0OnBez2 {
		return svalues[2]
	}
	if atBoundary1OnBez1 && atBoundary1OnBez2 {
		return svalues[3]
	}

	minI, minJ := minIj.unwrap()[0], minIj.unwrap()[1]
	newUmid := umin + (umax-umin)*(float64(minI)/float64(2*n))
	newVmid := vmin + (vmax-vmin)*(float64(minJ)/float64(2*m))

	// Subdivide
	results := [4][3]float64{
		minDistParam(
			bez1,
			bez2,
			[2]float64{umin, newUmid},
			[2]float64{vmin, newVmid},
			epsilon,
			alpha,
		),
		minDistParam(
			bez1,
			bez2,
			[2]float64{umin, newUmid},
			[2]float64{newVmid, vmax},
			epsilon,
			alpha,
		),
		minDistParam(
			bez1,
			bez2,
			[2]float64{newUmid, umax},
			[2]float64{vmin, newVmid},
			epsilon,
			alpha,
		),
		minDistParam(
			bez1,
			bez2,
			[2]float64{newUmid, umax},
			[2]float64{newVmid, vmax},
			epsilon,
			alpha,
		),
	}

	out := results[0]
	for _, res := range results[1:] {
		if math.IsNaN(res[0]) || res[0] < out[0] {
			out = res
		}
	}
	return out
}

func s(u, v float64, bez1, bez2 []Vec2) float64 {
	n := len(bez1) - 1
	m := len(bez2) - 1
	summand := 0.0
	for r := range 2*n + 1 {
		for k := range 2*m + 1 {
			summand +=
				dRk(r, k, bez1, bez2) * basisFunction(2*n, r, u) * basisFunction(2*m, k, v)
		}
	}
	return summand
}

func cRk(r, k int, bez1, bez2 []Vec2) float64 {
	var left Vec2
	{
		n := len(bez1) - 1
		upsilon := min(r, n)
		theta := r - min(n, r)
		items := bez1[:min(upsilon+1, len(bez1))]
		skip := min(theta, len(items))
		items = items[skip:]
		for i, item := range items {
			i += theta
			left = left.Add(item.Mul(float64(choose(n, i)*choose(n, r-i)) / float64(choose(2*n, r))))
		}
	}

	var right Vec2
	{
		m := len(bez2) - 1
		varsigma := min(k, m)
		sigma := k - min(m, k)
		items := bez2[:min(varsigma+1, len(bez2))]
		skip := min(sigma, len(items))
		items = items[skip:]
		for j, item := range items {
			j += skip
			right = right.Add(item.Mul(float64(choose(m, j)*choose(m, k-j)) / float64(choose(2*m, k))))
		}
	}

	return left.Dot(right)
}

func aR(r int, p []Vec2) float64 {
	n := len(p) - 1
	upsilon := min(r, n)
	theta := r - min(n, r)
	var sum float64
	for i := theta; i <= upsilon; i++ {
		dot := p[i].Dot(p[r-i]) // These are bounds checked by the sum limits
		factor := float64(choose(n, i)*choose(n, r-i)) / float64(choose(2*n, r))
		sum += dot * factor
	}
	return sum
}

func dRk(r, k int, bez1, bez2 []Vec2) float64 {
	// In the paper, B_k is used for the second factor, but it's the same thing
	return aR(r, bez1) + aR(k, bez2) - 2.0*cRk(r, k, bez1, bez2)
}

// Bezier basis function
func basisFunction(n, i int, u float64) float64 {
	return float64(choose(n, i)) * math.Pow(1.0-u, float64(n-i)) * math.Pow(u, float64(i))
}

// Binomial co-efficient, but returning zeros for values outside of domain
func choose(n, k int) uint32 {
	if k > n {
		return 0
	}
	p := 1
	bound := n - k
	for i := 1; i <= bound; i++ {
		p *= n
		p /= i
		n -= 1
	}
	return uint32(p)
}
