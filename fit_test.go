package curve

import (
	"fmt"
	"math"
	"testing"
)

func BenchmarkFitToBezpath(b *testing.B) {
	shape := CubicBez{
		P0: Pt(20, 40),
		P1: Pt(40, 80),
		P2: Pt(-40, 40),
		P3: Pt(42, 62),
	}
	for i := range 5 {
		acc := 1.0 / math.Pow(10, float64(i))
		b.Run(fmt.Sprintf("1e-%d", i), func(b *testing.B) {
			offset := NewCubicOffset(shape, -2, acc)
			b.ResetTimer()
			for range b.N {
				for range FitToBezPath(&offset, acc) {
				}
			}
		})
	}
}

func BenchmarkFitToBezpathOpt(b *testing.B) {
	shape := CubicBez{
		P0: Pt(20, 40),
		P1: Pt(40, 80),
		P2: Pt(-40, 40),
		P3: Pt(42, 62),
	}
	for i := range 5 {
		acc := 1.0 / math.Pow(10, float64(i))
		b.Run(fmt.Sprintf("1e-%d", i), func(b *testing.B) {
			offset := NewCubicOffset(shape, -2, acc)
			b.ResetTimer()
			for range b.N {
				FitToBezPathOpt(&offset, acc)
			}
		})
	}
}
