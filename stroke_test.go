package curve

import (
	"iter"
	"slices"
	"testing"
)

func BenchmarkDash(b *testing.B) {
	shape := Rect{X0: 5, Y0: 5, X1: 20, Y1: 20}.RoundedRect(RoundedRectRadii{2, 2, 2, 2})

	for range b.N {
		StrokePath(
			shape.PathElements(0.001),
			DefaultStroke.
				WithCaps(ButtCap).
				WithDashes(0, []float64{0.1}),
			StrokeOpts{},
			0.01)
	}
}

// Test cases adapted from https://github.com/linebender/vello/pull/388
func TestBrokenStrokes(t *testing.T) {
	brokenCubics := [][4]Point{
		{
			{465.24423, 107.11105},
			{475.50754, 107.11105},
			{475.50754, 107.11105},
			{475.50754, 107.11105},
		},
		{{0.0, -0.01}, {128.0, 128.001}, {128.0, -0.01}, {0.0, 128.001}}, // Near-cusp
		{{0.0, 0.0}, {0.0, -10.0}, {0.0, -10.0}, {0.0, 10.0}},            // Flat line with 180
		{{10.0, 0.0}, {0.0, 0.0}, {20.0, 0.0}, {10.0, 0.0}},              // Flat line with 2 180s
		{{39.0, -39.0}, {40.0, -40.0}, {40.0, -40.0}, {0.0, 0.0}},        // Flat diagonal with 180
		{{40.0, 40.0}, {0.0, 0.0}, {200.0, 200.0}, {0.0, 0.0}},           // Diag w/ an internal 180
		{{0.0, 0.0}, {1e-2, 0.0}, {-1e-2, 0.0}, {0.0, 0.0}},              // Circle
		// Flat line with no turns:
		{
			{400.75, 100.05},
			{400.75, 100.05},
			{100.05, 300.95},
			{100.05, 300.95},
		},
		{{0.5, 0.0}, {0.0, 0.0}, {20.0, 0.0}, {10.0, 0.0}},  // Flat line with 2 180s
		{{10.0, 0.0}, {0.0, 0.0}, {10.0, 0.0}, {10.0, 0.0}}, // Flat line with a 180
	}
	strokeStyle := DefaultStroke.WithWidth(30).WithCaps(ButtCap).WithJoin(MiterJoin)
	for _, cubic := range brokenCubics {
		path := CubicBez{cubic[0], cubic[1], cubic[2], cubic[3]}.PathElements(0.1)
		stroked := BezPath(slices.Collect(StrokePath(path, strokeStyle, StrokeOpts{}, 0.001)))
		if stroked.IsInf() {
			t.Errorf("got infinite stroke for %v", cubic)
		}
	}
}

func TestPathologicalStroke(t *testing.T) {
	curve := CubicBez{
		Pt(602.469, 286.585),
		Pt(641.975, 286.585),
		Pt(562.963, 286.585),
		Pt(562.963, 286.585),
	}
	path := curve.PathElements(0.1)
	strokeStyle := DefaultStroke.WithWidth(1.0)
	stroked := BezPath(slices.Collect(StrokePath(path, strokeStyle, StrokeOpts{}, 0.001)))
	if stroked.IsInf() {
		t.Error("got infinite stroke")
	}
}

func TestDashSequence(t *testing.T) {
	shape := Line{Pt(0.0, 0.0), Pt(21.0, 0.0)}
	dashes := []float64{1.0, 5.0, 2.0, 5.0}
	want := []PathSegment{
		Line{Pt(6.0, 0.0), Pt(8.0, 0.0)}.Seg(),
		Line{Pt(13.0, 0.0), Pt(14.0, 0.0)}.Seg(),
		Line{Pt(19.0, 0.0), Pt(21.0, 0.0)}.Seg(),
		Line{Pt(0.0, 0.0), Pt(1.0, 0.0)}.Seg(),
	}
	it := Segments(Dash(shape.PathElements(0.0), 0.0, dashes))
	got := slices.Collect(iter.Seq[PathSegment](it))
	diff(t, want, got)
}

func TestDashSequenceOffset(t *testing.T) {
	// Same as TestDashSequence, but with a dash offset of 3, which skips the first dash
	// and cuts into the first gap.
	shape := Line{Pt(0.0, 0.0), Pt(21.0, 0.0)}
	dashes := []float64{1.0, 5.0, 2.0, 5.0}
	want := []PathSegment{
		Line{Pt(3.0, 0.0), Pt(5.0, 0.0)}.Seg(),
		Line{Pt(10.0, 0.0), Pt(11.0, 0.0)}.Seg(),
		Line{Pt(16.0, 0.0), Pt(18.0, 0.0)}.Seg(),
	}
	it := Segments(Dash(shape.PathElements(0.0), 3.0, dashes))
	got := slices.Collect(iter.Seq[PathSegment](it))
	diff(t, want, got)
}

func TestStrokeWithMove(t *testing.T) {
	var p BezPath
	p.MoveTo(Pt(0, 0))
	p.LineTo(Pt(110., 506.))
	p.CubicTo(Pt(135.887, 506.953), Pt(169.789, 478.352), Pt(149.129, 502.719))

	for range StrokePath(p.Elements(), DefaultStroke, StrokeOpts{}, 0.1) {
	}
}
