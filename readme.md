# Analytical vector rendering
Just some experiments with computing the coverage of vectors by clipping vectors into closed-form shapes (lines, bezier quadratics and cubics) and calculating the coverage using [Green's Theorem](https://en.wikipedia.org/wiki/Green%27s_theorem). No triangulation is performed.

## Running examples
Ensure Rust is installed and run the examples via:
```
cargo run --release --example=render
cargo run --release --example=blur
cargo run --release --example=shape
cargo run --release --example=point_test
```

## Rendering
The 'render' example demonstrates per-pixel projection onto the shape and clipping, producing an alpha amount. The algorithm is designed to be ran on a Fragment shader of a GPU (coming soon) and will hopefully be quick enough to support rendering dynamic vectors in 3D with subpixel quality.

## Blur
High-quality blurs can be applied easily by changing the size of each projected pixel with little to no extra cost. The blurs use the same analytical coverage computation and are very high quality. The shape of each projected pixel can be changed to obtain different blur effects.
