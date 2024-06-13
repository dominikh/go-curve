package curve

import (
	"testing"

	"github.com/google/go-cmp/cmp"
)

func diff(t *testing.T, want, got any, opts ...cmp.Option) {
	t.Helper()
	if d := cmp.Diff(want, got, opts...); d != "" {
		t.Error(d)
	}
}
