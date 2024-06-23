package curve

import (
	"slices"
	"testing"

	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestQuadSpline(t *testing.T) {
	p1 := Pt(1, 1)
	p2 := Pt(2, 2)
	p3 := Pt(3, 3)
	p5 := Pt(5, 5)
	p8 := Pt(8, 8)
	tests := []struct {
		in  QuadBSpline
		out []QuadBez
	}{
		{make(QuadBSpline, 0), nil},
		{make(QuadBSpline, 1), nil},
		{make(QuadBSpline, 2), nil},
		{QuadBSpline{p1, p2, p3}, []QuadBez{{p1, p2, p3}}},
		{QuadBSpline{p1, p3, p5, p8}, []QuadBez{
			{p1, p3, p3.Midpoint(p5)},
			{p3.Midpoint(p5), p5, p8},
		}},
	}

	for _, tt := range tests {
		got := slices.Collect(tt.in.Quads())
		diff(t, got, tt.out, cmpopts.EquateEmpty())
	}
}
