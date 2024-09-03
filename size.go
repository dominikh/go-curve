package curve

import (
	"fmt"
	"math"
)

type Size struct {
	Width  float64
	Height float64
}

// Sz returns the size x×y.
func Sz(w, h float64) Size {
	return Size{
		Width:  w,
		Height: h,
	}
}

func (sz Size) String() string {
	return fmt.Sprintf("%g×%g", sz.Width, sz.Height)
}

func (sz Size) AsVec2() Vec2 {
	return Vec2{
		X: sz.Width,
		Y: sz.Height,
	}
}

func (sz Size) MaxSide() float64 {
	return max(sz.Width, sz.Height)
}

func (sz Size) MinSide() float64 {
	return min(sz.Width, sz.Height)
}

func (sz Size) Area() float64 {
	return sz.Width * sz.Height
}

func (sz Size) Clamp(minSz, maxSz Size) Size {
	w := min(max(sz.Width, minSz.Width), maxSz.Width)
	h := min(max(sz.Height, minSz.Height), maxSz.Height)
	return Size{
		Width:  w,
		Height: h,
	}
}

func (sz Size) Splat() (w float64, h float64) {
	return sz.Width, sz.Height
}

// Round returns a new size with width and height rounded to the nearest integers.
func (sz Size) Round() Size {
	return Size{
		Width:  math.Round(sz.Width),
		Height: math.Round(sz.Height),
	}
}

// Ceil returns a new size with width and height rounded up to the nearest integers.
func (sz Size) Ceil() Size {
	return Size{
		Width:  math.Ceil(sz.Width),
		Height: math.Ceil(sz.Height),
	}
}

// Floor returns a new size with width and height rounded down to the nearest integers.
func (sz Size) Floor() Size {
	return Size{
		Width:  math.Floor(sz.Width),
		Height: math.Floor(sz.Height),
	}
}

// Expand returns a new size with width and height rounded away from zero to the
// nearest integers.
func (sz Size) Expand() Size {
	return Size{
		Width:  expand(sz.Width),
		Height: expand(sz.Height),
	}
}

// Trunc returns a new size with width and height rounded towards zero to the nearest
// integers.
func (sz Size) Trunc() Size {
	return Size{
		Width:  math.Trunc(sz.Width),
		Height: math.Trunc(sz.Height),
	}
}

func (sz Size) AspectRatio() float64 {
	// XXX why height / width and not width / height?
	return sz.Height / sz.Width
}

// IsInf reports whether at least one of width and height is infinite.
func (sz Size) IsInf() bool {
	return math.IsInf(sz.Width, 0) || math.IsInf(sz.Height, 0)
}

// IsNaN reports whether at least one of width and height is NaN.
func (sz Size) IsNaN() bool {
	return math.IsNaN(sz.Width) || math.IsNaN(sz.Height)
}

// Scale multiplies sz by f.
func (sz Size) Scale(f float64) Size {
	return Size{
		Width:  sz.Width * f,
		Height: sz.Height * f,
	}
}

func (sz Size) Add(v Vec2) Size {
	return Size{
		Width:  sz.Width + v.X,
		Height: sz.Height + v.Y,
	}
}
