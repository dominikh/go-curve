package curve

import (
	"slices"
	"testing"
)

func TestSimplifyLinesCorner(t *testing.T) {
	var p BezPath
	p.MoveTo(Pt(1.0, 2.0))
	p.LineTo(Pt(3.0, 4.0))
	p.LineTo(Pt(10.0, 5.0))
	simplified := BezPath(slices.Collect(Simplify(p.Elements(), 1.0, DefaultSimplifyOptions)))
	diff(t, p, simplified)
}
