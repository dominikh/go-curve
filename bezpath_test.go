package curve

import (
	"iter"
	"slices"
	"testing"
)

func TestElementsToSegmentsClosePathReferstoLastMove(t *testing.T) {
	last := func(seq iter.Seq[PathSegment]) PathSegment {
		var el PathSegment
		for el = range seq {
		}
		return el
	}
	var p BezPath
	p.MoveTo(Pt(5.0, 5.0))
	p.LineTo(Pt(15.0, 15.0))
	p.MoveTo(Pt(10.0, 10.0))
	p.LineTo(Pt(15.0, 15.0))
	p.ClosePath()

	want := Line{Pt(15, 15), Pt(10, 10)}.Seg()
	diff(t, want, last(p.Segments()))
}

func TestContains(t *testing.T) {
	var path BezPath
	path.MoveTo(Pt(0.0, 0.0))
	path.LineTo(Pt(1.0, 1.0))
	path.LineTo(Pt(2.0, 0.0))
	path.ClosePath()
	if w := path.Winding(Pt(1, 0.5)); w != -1 {
		t.Errorf("got winding %v, want -1", w)
	}
}

func TestReverseUnclosed(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(10, 10)),
			QuadTo(Pt(40, 40), Pt(60, 10)),
			LineTo(Pt(100, 10)),
			CubicTo(Pt(125, 10), Pt(150, 50), Pt(125, 60)),
		},
		[]PathElement{
			MoveTo(Pt(125, 60)),
			CubicTo(Pt(150, 50), Pt(125, 10), Pt(100, 10)),
			LineTo(Pt(60, 10)),
			QuadTo(Pt(40, 40), Pt(10, 10)),
		},
	)
}

func TestReverseClosedTriangle(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(100, 100)),
			LineTo(Pt(150, 200)),
			LineTo(Pt(50, 200)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(50, 200)),
			LineTo(Pt(150, 200)),
			LineTo(Pt(100, 100)),
			ClosePath(),
		},
	)
}

func TestReverseClosedShape(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(125, 100)),
			QuadTo(Pt(200, 150), Pt(175, 300)),
			CubicTo(Pt(150, 150), Pt(50, 150), Pt(25, 300)),
			QuadTo(Pt(0, 150), Pt(75, 100)),
			LineTo(Pt(100, 50)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(100, 50)),
			LineTo(Pt(75, 100)),
			QuadTo(Pt(0, 150), Pt(25, 300)),
			CubicTo(Pt(50, 150), Pt(150, 150), Pt(175, 300)),
			QuadTo(Pt(200, 150), Pt(125, 100)),
			ClosePath(),
		},
	)
}

func TestReverseMultipleSubpaths(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(10, 10)),
			QuadTo(Pt(40, 40), Pt(60, 10)),
			LineTo(Pt(100, 10)),
			CubicTo(Pt(125, 10), Pt(150, 50), Pt(125, 60)),
			MoveTo(Pt(100, 100)),
			LineTo(Pt(150, 200)),
			LineTo(Pt(50, 200)),
			ClosePath(),
			MoveTo(Pt(125, 100)),
			QuadTo(Pt(200, 150), Pt(175, 300)),
			CubicTo(Pt(150, 150), Pt(50, 150), Pt(25, 300)),
			QuadTo(Pt(0, 150), Pt(75, 100)),
			LineTo(Pt(100, 50)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(125, 60)),
			CubicTo(Pt(150, 50), Pt(125, 10), Pt(100, 10)),
			LineTo(Pt(60, 10)),
			QuadTo(Pt(40, 40), Pt(10, 10)),
			MoveTo(Pt(50, 200)),
			LineTo(Pt(150, 200)),
			LineTo(Pt(100, 100)),
			ClosePath(),
			MoveTo(Pt(100, 50)),
			LineTo(Pt(75, 100)),
			QuadTo(Pt(0, 150), Pt(25, 300)),
			CubicTo(Pt(50, 150), Pt(150, 150), Pt(175, 300)),
			QuadTo(Pt(200, 150), Pt(125, 100)),
			ClosePath(),
		},
	)
}

func TestReverseMultipleMoves(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(2.0, 2.0)),
			MoveTo(Pt(3.0, 3.0)),
			ClosePath(),
			MoveTo(Pt(4.0, 4.0)),
		},
		[]PathElement{
			MoveTo(Pt(2.0, 2.0)),
			MoveTo(Pt(3.0, 3.0)),
			ClosePath(),
			MoveTo(Pt(4.0, 4.0)),
		},
	)
}

func TestControlBox(t *testing.T) {
	// a sort of map ping looking thing drawn with a single cubic
	// cbox is wildly different than tight box
	var p BezPath
	p.MoveTo(Pt(200, 300))
	p.CubicTo(Pt(50, 50), Pt(350, 50), Pt(200, 300))
	want := Rect{50, 50, 350, 300}
	diff(t, p.ControlBox(), want)
	if a, b := p.ControlBox().Area(), p.BoundingBox().Area(); a < b {
		t.Errorf("control box is smaller than bounding box: %v < %v", a, b)
	}
}

func TestGetSegment(t *testing.T) {
	// Segment(i) should produce the same results as Segments()[i-1]
	var circle BezPath = slices.Collect(Circle{Pt(10.0, 10.0), 2.0}.PathElements(0.1))
	segments := slices.Collect(circle.Segments())
	var getSegments []PathSegment
	for i := 1; ; i++ {
		seg, ok := circle.Segment(i)
		if !ok {
			break
		}
		getSegments = append(getSegments, seg)
	}
	diff(t, segments, getSegments)
}

// The following are direct port of fonttools'
// reverseContourPen_test.py::test_reverse_pen, adapted to rust, excluding
// test cases that don't apply because we don't implement
// outputImpliedClosingLine=False.
// https://github.com/fonttools/fonttools/blob/85c80be/Tests/pens/reverseContourPen_test.py#L6-L467

func TestReverseClosedLastLineNotOnMove(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			LineTo(Pt(1.0, 1.0)),
			LineTo(Pt(2.0, 2.0)),
			LineTo(Pt(3.0, 3.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(3.0, 3.0)),
			LineTo(Pt(2.0, 2.0)),
			LineTo(Pt(1.0, 1.0)),
			LineTo(Pt(0.0, 0.0)), // closing line NOT implied
			ClosePath(),
		},
	)
}

func TestReverseClosedLastLineOverlapsMove(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			LineTo(Pt(1.0, 1.0)),
			LineTo(Pt(2.0, 2.0)),
			LineTo(Pt(0.0, 0.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			LineTo(Pt(2.0, 2.0)),
			LineTo(Pt(1.0, 1.0)),
			LineTo(Pt(0.0, 0.0)), // closing line NOT implied
			ClosePath(),
		},
	)
}

func TestReverseClosedDuplicateLineFollowingMove(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			LineTo(Pt(0.0, 0.0)),
			LineTo(Pt(1.0, 1.0)),
			LineTo(Pt(2.0, 2.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(2.0, 2.0)),
			LineTo(Pt(1.0, 1.0)),
			LineTo(Pt(0.0, 0.0)), // duplicate line retained
			LineTo(Pt(0.0, 0.0)),
			ClosePath(),
		},
	)
}

func TestReverseClosedTwoLines(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			LineTo(Pt(1.0, 1.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(1.0, 1.0)),
			LineTo(Pt(0.0, 0.0)), // closing line NOT implied
			ClosePath(),
		},
	)
}

func TestReverseClosedLastCurveOverlapsMove(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			CubicTo(Pt(1.0, 1.0), Pt(2.0, 2.0), Pt(3.0, 3.0)),
			CubicTo(Pt(4.0, 4.0), Pt(5.0, 5.0), Pt(0.0, 0.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)), // no extra lineTo added here
			CubicTo(Pt(5.0, 5.0), Pt(4.0, 4.0), Pt(3.0, 3.0)),
			CubicTo(Pt(2.0, 2.0), Pt(1.0, 1.0), Pt(0.0, 0.0)),
			ClosePath(),
		},
	)
}

func TestReverseClosedLastCurveNotOnMove(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			CubicTo(Pt(1.0, 1.0), Pt(2.0, 2.0), Pt(3.0, 3.0)),
			CubicTo(Pt(4.0, 4.0), Pt(5.0, 5.0), Pt(6.0, 6.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(6.0, 6.0)), // the previously implied line
			CubicTo(Pt(5.0, 5.0), Pt(4.0, 4.0), Pt(3.0, 3.0)),
			CubicTo(Pt(2.0, 2.0), Pt(1.0, 1.0), Pt(0.0, 0.0)),
			ClosePath(),
		},
	)
}

func TestReverseClosedLineCurveLine(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			LineTo(Pt(1.0, 1.0)), // this line...
			CubicTo(Pt(2.0, 2.0), Pt(3.0, 3.0), Pt(4.0, 4.0)),
			CubicTo(Pt(5.0, 5.0), Pt(6.0, 6.0), Pt(7.0, 7.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(7.0, 7.0)),
			CubicTo(Pt(6.0, 6.0), Pt(5.0, 5.0), Pt(4.0, 4.0)),
			CubicTo(Pt(3.0, 3.0), Pt(2.0, 2.0), Pt(1.0, 1.0)),
			LineTo(Pt(0.0, 0.0)), // ... does NOT become implied
			ClosePath(),
		},
	)
}

func TestReverseClosedLastQuadOverlapsMove(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			QuadTo(Pt(1.0, 1.0), Pt(2.0, 2.0)),
			QuadTo(Pt(3.0, 3.0), Pt(0.0, 0.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)), // no extra lineTo added here
			QuadTo(Pt(3.0, 3.0), Pt(2.0, 2.0)),
			QuadTo(Pt(1.0, 1.0), Pt(0.0, 0.0)),
			ClosePath(),
		},
	)
}

func TestReverseClosedLastQuadNotOnMove(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			QuadTo(Pt(1.0, 1.0), Pt(2.0, 2.0)),
			QuadTo(Pt(3.0, 3.0), Pt(4.0, 4.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(4.0, 4.0)), // the previously implied line
			QuadTo(Pt(3.0, 3.0), Pt(2.0, 2.0)),
			QuadTo(Pt(1.0, 1.0), Pt(0.0, 0.0)),
			ClosePath(),
		},
	)
}

func TestReverseClosedLineQuadLine(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			LineTo(Pt(1.0, 1.0)), // this line...
			QuadTo(Pt(2.0, 2.0), Pt(3.0, 3.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(3.0, 3.0)),
			QuadTo(Pt(2.0, 2.0), Pt(1.0, 1.0)),
			LineTo(Pt(0.0, 0.0)), // ... does NOT become implied
			ClosePath(),
		},
	)
}

func TestReverseEmpty(t *testing.T) {
	reverseHelper(t, []PathElement{}, []PathElement{})
}

func TestReverseSinglePoint(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{MoveTo(Pt(0.0, 0.0))},
		[]PathElement{MoveTo(Pt(0.0, 0.0))},
	)
}

func TestReverseSinglePointClosed(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{MoveTo(Pt(0.0, 0.0)), ClosePath()},
		[]PathElement{MoveTo(Pt(0.0, 0.0)), ClosePath()},
	)
}

func TestReverseSingleLineOpen(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			LineTo(Pt(1.0, 1.0)),
		},
		[]PathElement{
			MoveTo(Pt(1.0, 1.0)),
			LineTo(Pt(0.0, 0.0)),
		},
	)
}

func TestReverseSingleCurveOpen(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			CubicTo(Pt(1.0, 1.0), Pt(2.0, 2.0), Pt(3.0, 3.0)),
		},
		[]PathElement{
			MoveTo(Pt(3.0, 3.0)),
			CubicTo(Pt(2.0, 2.0), Pt(1.0, 1.0), Pt(0.0, 0.0)),
		},
	)
}

func TestReverseCurveLineOpen(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			CubicTo(Pt(1.0, 1.0), Pt(2.0, 2.0), Pt(3.0, 3.0)),
			LineTo(Pt(4.0, 4.0)),
		},
		[]PathElement{
			MoveTo(Pt(4.0, 4.0)),
			LineTo(Pt(3.0, 3.0)),
			CubicTo(Pt(2.0, 2.0), Pt(1.0, 1.0), Pt(0.0, 0.0)),
		},
	)
}

func TestReverseLineCurveOpen(t *testing.T) {
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			LineTo(Pt(1.0, 1.0)),
			CubicTo(Pt(2.0, 2.0), Pt(3.0, 3.0), Pt(4.0, 4.0)),
		},
		[]PathElement{
			MoveTo(Pt(4.0, 4.0)),
			CubicTo(Pt(3.0, 3.0), Pt(2.0, 2.0), Pt(1.0, 1.0)),
			LineTo(Pt(0.0, 0.0)),
		},
	)
}

func TestReverseDuplicatePointAfterMove(t *testing.T) {
	// Test case from: https://github.com/googlei18n/cu2qu/issues/51#issue-179370514
	// Simplified to only use atomic QuadTo (no QuadSplines).
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(848.0, 348.0)),
			LineTo(Pt(848.0, 348.0)),
			QuadTo(Pt(848.0, 526.0), Pt(449.0, 704.0)),
			QuadTo(Pt(848.0, 171.0), Pt(848.0, 348.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(848.0, 348.0)),
			QuadTo(Pt(848.0, 171.0), Pt(449.0, 704.0)),
			QuadTo(Pt(848.0, 526.0), Pt(848.0, 348.0)),
			LineTo(Pt(848.0, 348.0)),
			ClosePath(),
		},
	)
}

func TestReverseDuplicatePointAtEnd(t *testing.T) {
	// Test case from: https://github.com/googlefonts/fontmake/issues/572
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 651.0)),
			LineTo(Pt(0.0, 101.0)),
			LineTo(Pt(0.0, 101.0)),
			LineTo(Pt(0.0, 651.0)),
			LineTo(Pt(0.0, 651.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(0.0, 651.0)),
			LineTo(Pt(0.0, 651.0)),
			LineTo(Pt(0.0, 101.0)),
			LineTo(Pt(0.0, 101.0)),
			LineTo(Pt(0.0, 651.0)),
			ClosePath(),
		},
	)
}

func TestReverseLines(t *testing.T) {
	// https://github.com/fonttools/fonttools/blob/bf265ce49e0cae6f032420a4c80c31d8e16285b8/Tests/pens/reverseContourPen_test.py#L7
	reverseHelper(
		t,
		[]PathElement{
			MoveTo(Pt(0.0, 0.0)),
			LineTo(Pt(1.0, 1.0)),
			LineTo(Pt(2.0, 2.0)),
			LineTo(Pt(3.0, 3.0)),
			ClosePath(),
		},
		[]PathElement{
			MoveTo(Pt(3, 3)),
			LineTo(Pt(2, 2)),
			LineTo(Pt(1, 1)),
			LineTo(Pt(0, 0)),
			ClosePath(),
		},
	)
}

func reverseHelper(t *testing.T, contour, want []PathElement) {
	t.Helper()

	var got []PathElement = BezPath(contour).ReverseSubpaths()
	diff(t, got, want)
}
