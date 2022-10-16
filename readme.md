# Analytical vector rendering
Just some experiments with computing the coverage of vectors by clipping vectors into closed-form shapes (lines, bezier quadratics and cubics) and calculating the coverage using [Green's Theorem](https://en.wikipedia.org/wiki/Green%27s_theorem). No triangulation is performed.

## Running examples
Ensure Rust is installed and run the examples via:
```
cargo run --release --example=render
cargo run --release --example=shape
cargo run --release --example=point_test
```

## Status
The 'render' example demonstrates per-pixel projection onto the shape and clipping, producing an alpha amount. The green pixels show the errors currently present in the algorithm, where the area is not computed correctly (likely due to precission issues). The algorithm is designed to be ran on a Fragment shader of a GPU and would hopefully be quick enough to support rendering dynamic vectors in 3D with subpixel quality.
