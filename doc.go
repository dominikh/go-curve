// Package curve provides primitives and routines for 2D shapes, curves, and
// paths. It was designed to serve the needs of 2D graphics applications, but it
// is intended to be general enough to be useful for other applications.
//
// # Kurbo
//
// This package is a manual, idiomatic Go port of the [kurbo] Rust crate. Kurbo
// contains several novel approaches to problems in 2D curve design,
// particularly stroke expansion. This package's development will closely follow
// kurbo's.
//
// We may introduce functionality of our own, but this will likely be limited to
// simple features.
//
// # Features
//
// We provide the following notable features:
//
//   - Curve fitting (see [FitToBezPath])
//   - Stroke expansion (see [StrokePath])
//   - Dashing (see [Dash])
//   - Curve simplification (see [Simplify])
//   - Affine transformations (see [Affine])
//   - Approximating cubic Béziers with quadratic splines (see [CubicBez.ApproxQuadSpline])
//   - Flattening paths to lines (see [Flatten])
//
// # Shapes, parametric curves, and paths
//
// The three core primitives of this package are shapes, curves, and paths.
//
// [Shape] describes geometric shapes that have a perimeter and a bounding box,
// and that can be converted to a series of path elements. [ClosedShape]
// describes those shapes that additionally have a closed outline (such as
// rectangles, but not lines) and thus have an area and can compute a point's
// [winding number].
//
// This package includes the following shapes:
//   - [Arc]
//   - [Circle]
//   - [CircleSegment]
//   - [CubicBez]
//   - [Ellipse]
//   - [Line]
//   - [QuadBez]
//   - [Rect]
//   - [RoundedRect]
//
// [ParametricCurve] describes parametrized curves. These curves can be
// evaluated at t ∈ [0, 1] and return pairs of (x, y) values, commonly
// interpreted as points in a 2D Cartesian coordinate system. The simplest
// parametric curve is the [Line], whose evaluation is a simple linear
// interpolation between its start and end points. More complex curves are the
// quadratic and cubic Béziers, for example.
//
// [Arclener] is an optional interface implemented by curves that can compute
// their length. [SignedAreaer] is an optional interface implemented by curves
// that can compute their signed area. [ArclenSolver] is an optional interface
// implemented by curves that can efficiently solve for t given an arc length.
//
// [FittableCurve] is closely related to [ParametricCurve] and describes
// parametrized curves in a way that allows efficient curve fitting. Using curve
// fitting, arbitrary curves can be approximated with paths. See the example of
// [FitToBezPathOpt] for a custom curve that implements a quartic function and
// that gets approximated by cubic Béziers.
//
// This package includes the following curves:
//   - [Line]
//   - [CubicBez]
//   - [PathSegment] (as it is a wrapper for [Line], [CubicBez], and [QuadBez])
//   - [QuadBez]
//
// # Bézier paths
//
// Bézier paths consist of lines and quadratic and cubic Béziers. All shapes can
// represent themselves as Bézier paths, and more complex curves can be
// approximated by them.
//
// [BezPath] represents Bézier paths as a slice of path elements and provides
// methods for building paths as well as for any path manipulation that needs
// access to the collection of path elements. The two other representations are
// iter.Seq[PathElement] and iter.Seq[PathSegment], see the Iterators section
// for more on that.
//
// # Path elements and segments
//
// This package provides two representations for paths: [PathElement] and
// [PathSegment]. Path elements are akin to drawing commands in graphics APIs
// like PostScript, consisting of pen moves ([MoveTo]) and various drawing
// commands ([LineTo], [QuadTo], etc.) Each command moves the current position
// of the pen, which acts as the start position of the following drawing
// command.
//
// Segments, on the other hand, are self-contained descriptions of a portion of
// the path, containing explicit start points.
//
// Using [Elements] and [Segments], you can freely convert between the two
// representations.
//
// The package supports [MoveTo], [LineTo], [QuadTo], [CubicTo], and [ClosePath]
// line elements, and correspondingly has path segments for lines, quadratic
// Béziers, and cubic Béziers.
//
// # Cubic Béziers
//
// The primary curve used by this package is the cubic Bézier. Functions that
// convert shapes to paths do so by approximating them with cubic Béziers (and
// lines). Similarly, curve fitting approximates curves with cubic Béziers (and
// lines). This is done because cubic Béziers hit the sweet spot of
// expressiveness and computability.
//
// However, we also provide functions for working with quadratic Béziers,
// including approximating cubic Béziers with quadratic Béziers. This
// functionality is particularly useful for font design, as TrueType fonts
// (unlike OpenType) use quadratic Béziers.
//
// # Iterators
//
// Many functions in this package can operate on individual path elements or
// segments at a time and don't need random access. Similarly, many functions
// don't have to remember or modify the sequence of path elements or segments
// they produce. In those cases, they accept and return iterators, to avoid
// having to allocate slices.
//
// Functions that cannot work with iterators directly will instead accept or
// return slices, to make it clear that they allocate. You can use
// [slices.Collect] to turn iterators into slices, and [slices.Values] to turn
// slices into iterators.
//
// A notable example are [FitToBezPath] and [FitToBezPathOpt]. The former
// returns an iterator, as it can produce elements one at a time. The latter has
// to do backtracking to find the optimal solution, and thus returns the slice
// it had to build up.
//
// Many functions, such as [Dash] or [Simplify] can be considered adapters that
// map from one sequence to another.
//
// # Literature
//
// This package makes use of the following ideas:
//   - [A Primer on Bézier Curves]
//   - [Algorithm 1010: Boosting Efficiency in Solving Quartic Equations with No Compromise in Accuracy] by Orellana and De Michele
//   - [An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality] by Oliveira and Takahashi
//   - [Approximate a circle with cubic Bézier curves] by Spencer Mortensen
//   - [Calculating Area of Closed Curves in ℝ²]
//   - [Fitting cubic Bézier curves] by Raph Levien
//   - [Flattening quadratic Béziers] by Raph Levien
//   - [Green's theorem]
//   - [How to solve a cubic equation, revisited] by Christoph Peters
//   - [Parallel curves of cubic Béziers] by Raph Levien
//   - [Simplifying Bézier paths] by Raph Levien
//   - https://github.com/Pomax/BezierInfo-2/issues/238
//
// [A Primer on Bézier Curves]: https://pomax.github.io/bezierinfo/
// [Algorithm 1010: Boosting Efficiency in Solving Quartic Equations with No Compromise in Accuracy]: https://cristiano-de-michele.netlify.app/publication/orellana-2020/orellana-2020.pdf
// [An Enhancement of the Bisection Method Average Performance Preserving Minmax Optimality]: https://dl.acm.org/doi/10.1145/3423597
// [Approximate a circle with cubic Bézier curves]: https://spencermortensen.com/articles/bezier-circle/
// [Fitting cubic Bézier curves]: https://raphlinus.github.io/curves/2021/03/11/bezier-fitting.html
// [Flattening quadratic Béziers]: https://raphlinus.github.io/graphics/curves/2019/12/23/flatten-quadbez.html
// [Green's theorem]: https://en.wikipedia.org/wiki/Green%27s_theorem
// [How to solve a cubic equation, revisited]: https://momentsingraphics.de/CubicRoots.html
// [Parallel curves of cubic Béziers]: https://raphlinus.github.io/curves/2022/09/09/parallel-beziers.html
// [Simplifying Bézier paths]: https://raphlinus.github.io/curves/2023/04/18/bezpath-simplify.html
// [kurbo]: https://github.com/linebender/kurbo
// [winding number]: https://en.wikipedia.org/wiki/Winding_number
// [Calculating Area of Closed Curves in ℝ²]: http://ich.deanmcnamee.com/graphics/2016/03/30/CurveArea.html
package curve
