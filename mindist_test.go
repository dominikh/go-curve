package curve

import (
	"testing"

	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestChoose(t *testing.T) {
	diff(t, choose(6, 0), uint32(1))
	diff(t, choose(6, 1), uint32(6))
	diff(t, choose(6, 2), uint32(15))
}

func TestD_rk(t *testing.T) {
	bez1 := []Vec2{
		Vec(129.0, 139.0),
		Vec(190.0, 139.0),
		Vec(201.0, 364.0),
		Vec(90.0, 364.0),
	}
	bez2 := []Vec2{
		Vec(309.0, 159.0),
		Vec(178.0, 159.0),
		Vec(215.0, 408.0),
		Vec(309.0, 408.0),
	}
	b := aR(1, bez2)
	diff(t, b, 80283.0, cmpopts.EquateApprox(0, 0.005))
	d := dRk(0, 1, bez1, bez2)
	diff(t, d, 9220.0, cmpopts.EquateApprox(0, 0.005))
}

func TestMinDist(t *testing.T) {
	bez1 := CubicBez{
		Pt(129.0, 139.0),
		Pt(190.0, 139.0),
		Pt(201.0, 364.0),
		Pt(90.0, 364.0),
	}.Seg()
	bez2 := CubicBez{
		Pt(309.0, 159.0),
		Pt(178.0, 159.0),
		Pt(215.0, 408.0),
		Pt(309.0, 408.0),
	}.Seg()
	mindist := bez1.MinDist(bez2, 0.001)
	diff(t, mindist.Distance, 50.9966, cmpopts.EquateApprox(0, 0.5))
}

func TestMinDistOverflow(t *testing.T) {
	bez1 := CubicBez{
		Pt(232.0, 126.0),
		Pt(134.0, 126.0),
		Pt(139.0, 232.0),
		Pt(141.0, 301.0),
	}.Seg()
	bez2 := Line{Pt(359.0, 416.0), Pt(367.0, 755.0)}.Seg()
	mindist := bez1.MinDist(bez2, 0.001)
	diff(t, mindist.Distance, 246.4731222669117, cmpopts.EquateApprox(0, 0.5))
}

func TestMinDistOutOfOrder(t *testing.T) {
	bez1 := CubicBez{
		Pt(287.0, 182.0),
		Pt(346.0, 277.0),
		Pt(356.0, 299.0),
		Pt(359.0, 416.0),
	}.Seg()
	bez2 := Line{Pt(141.0, 301.0), Pt(152.0, 709.0)}.Seg()
	mindist1 := bez1.MinDist(bez2, 0.5)
	mindist2 := bez2.MinDist(bez1, 0.5)
	diff(t, mindist1.Distance, mindist2.Distance, cmpopts.EquateApprox(0, 0.5))
}
