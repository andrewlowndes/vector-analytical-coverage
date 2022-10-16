use gen_iter::gen_iter;

pub type Pixel = (i32, i32);
pub type Pixel3 = (i32, i32, i32);
pub type PixelAA = (i32, i32, u8);

/**
 * Semi-automated port of easyFilter to Rust using generators/iterators
 * @author Zingl Alois
 * @date 22.08.2016
 * @version 1.2
*/
pub fn plot_line(mut x0: i32, mut y0: i32, x1: i32, y1: i32) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        let dx = (x1 - x0).abs();
        let sx = if x0 < x1 { 1 } else { -(1) };
        let dy = -(y1 - y0).abs();
        let sy = if y0 < y1 { 1 } else { -(1) };
        let mut err = dx + dy;
        // error value e_xy
        let mut e2;
        loop {
            yield (x0, y0);
            e2 = 2 * err;
            if e2 >= dy {
                // e_xy+e_x > 0
                if x0 == x1 {
                    break;
                }
                err += dy;
                x0 += sx
            }
            if e2 > dx {
                continue;
            }
            // e_xy+e_y < 0
            if y0 == y1 {
                break;
            }
            err += dx;
            y0 += sy
        }
    })
}

pub fn plot_line3d(
    mut x0: i32,
    mut y0: i32,
    mut z0: i32,
    mut x1: i32,
    mut y1: i32,
    mut z1: i32,
) -> impl Iterator<Item = Pixel3> {
    gen_iter!(move {
        let dx = (x1 - x0).abs();
        let sx = if x0 < x1 { 1 } else { -(1) };
        let dy = (y1 - y0).abs();
        let sy = if y0 < y1 { 1 } else { -(1) };
        let dz = (z1 - z0).abs();
        let sz = if z0 < z1 { 1 } else { -(1) };
        let dm = if dx > dy && dx > dz {
            dx
        } else if dy > dz {
            dy
        } else {
            dz
        };
        let mut i = dm;
        // max diff
        z1 = dm / 2;
        y1 = z1;
        x1 = y1;
        // error offset
        loop {
            yield (x0, y0, z0);
            let fresh0 = i;
            i -= 1;
            if fresh0 == 0 {
                break;
            }
            x1 -= dx;
            if x1 < 0 {
                x1 += dm;
                x0 += sx
            }
            y1 -= dy;
            if y1 < 0 {
                y1 += dm;
                y0 += sy
            }
            z1 -= dz;
            if z1 < 0 {
                z1 += dm;
                z0 += sz
            }
        }
    })
}

pub fn plot_ellipse(xm: i32, ym: i32, a: i32, b: i32) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        let mut x = -a;
        let mut y = 0;
        // II. quadrant from bottom left to top right
        let mut e2: i64 = b as i64 * b as i64;
        // error of 1.step
        let mut err: i64 = x as i64 * (2_i64 * e2 + x as i64) + e2;
        loop {
            yield (xm - x, ym + y);
            // I. Quadrant
            yield (xm + x, ym + y);
            // II. Quadrant
            yield (xm + x, ym - y);
            // III. Quadrant
            yield (xm - x, ym - y);
            // IV. Quadrant
            e2 = 2_i64 * err;
            if e2 >= (x * 2 + 1) as i64 * b as i64 * b as i64 {
                // e_xy+e_x > 0
                x += 1;
                err += (x * 2 + 1) as i64 * b as i64 * b as i64
            }
            if e2 <= (y * 2 + 1) as i64 * a as i64 * a as i64 {
                // e_xy+e_y < 0
                y += 1;
                err += (y * 2 + 1) as i64 * a as i64 * a as i64
            }
            if x > 0 {
                break;
            }
        }
        loop {
            let fresh1 = y;
            y += 1;
            if fresh1 >= b {
                break;
            }
            // too early stop of flat ellipses a=1,
            yield (xm, ym + y);
            // -> finish tip of ellipse
            yield (xm, ym - y);
        }
    })
}

pub fn plot_optimized_ellipse(xm: i32, ym: i32, a: i32, b: i32) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        let mut x: i64 = -a as i64;
        let mut y: i64 = 0_i64;
        // II. quadrant from bottom left to top right
        let mut e2: i64 = b as i64;
        let mut dx: i64 = (1_i64 + 2_i64 * x) * e2 * e2;
        // error increment
        let mut dy: i64 = x * x;
        // error of 1.step
        let mut err: i64 = dx + dy;
        loop {
            yield ((xm as i64 - x) as i32, (ym as i64 + y) as i32);
            // y step
            yield ((xm as i64 + x) as i32, (ym as i64 + y) as i32);
            yield ((xm as i64 + x) as i32, (ym as i64 - y) as i32);
            yield ((xm as i64 - x) as i32, (ym as i64 - y) as i32);
            e2 = 2_i64 * err;
            if e2 >= dx {
                x += 1;
                dx += 2_i64 * b as i64 * b as i64;
                err += dx
            }
            // I. Quadrant
            // II. Quadrant
            // III. Quadrant
            // IV. Quadrant
            // x step
            if e2 <= dy {
                y += 1;
                dy += 2_i64 * a as i64 * a as i64;
                err += dy
            }
            if x > 0_i64 {
                break;
            }
        }
        loop {
            let fresh2 = y;
            y += 1;
            if fresh2 >= b as i64 {
                break;
            }
            // too early stop for flat ellipses with a=1,
            yield (xm, (ym as i64 + y) as i32);
            // -> finish tip of ellipse
            yield (xm, (ym as i64 - y) as i32);
        }
    })
}

pub fn plot_circle(xm: i32, ym: i32, mut r: i32) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        let mut x = -r;
        let mut y = 0;
        let mut err = 2 - 2 * r;
        loop {
            // bottom left to top right
            yield (xm - x, ym + y);
            // -> x-step now
            yield (xm - y, ym - x);
            yield (xm + x, ym - y);
            yield (xm + y, ym + x);
            r = err;
            if r <= y {
                y += 1;
                err += y * 2 + 1
            }
            // I. Quadrant +x +y
            // II. Quadrant -x +y
            // III. Quadrant -x -y
            // IV. Quadrant +x -y
            // e_xy+e_y < 0
            if r > x || err > y {
                // e_xy+e_x > 0 or no 2nd y-step
                x += 1;
                err += x * 2 + 1
            }
            if x >= 0 {
                break;
            }
        }
    })
}

pub fn plot_ellipse_rect(
    mut x0: i32,
    mut y0: i32,
    mut x1: i32,
    mut y1: i32,
) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        // rectangular parameter enclosing the ellipse
        let mut a: i64 = (x1 - x0).abs() as i64;
        let b: i64 = (y1 - y0).abs() as i64;
        let mut b1: i64 = b & 1_i64;
        // diameter
        let mut dx: f64 = 4_f64 * (1.0f64 - a as f64) * b as f64 * b as f64;
        let mut dy: f64 = (4_i64 * (b1 + 1_i64) * a * a) as f64;
        // error increment
        let mut err: f64 = dx + dy + (b1 * a * a) as f64;
        let mut e2: f64;
        // error of 1.step
        if x0 > x1 {
            x0 = x1;
            x1 = (x1 as i64 + a) as i32
        }
        // if called with swapped points
        if y0 > y1 {
            y0 = y1
        }
        // .. exchange them
        y0 = (y0 as i64 + (b + 1_i64) / 2_i64) as i32;
        y1 = (y0 as i64 - b1) as i32;
        // starting pixel
        a = 8_i64 * a * a;
        b1 = 8_i64 * b * b;
        loop {
            yield (x1, y0);
            // x step
            yield (x0, y0);
            yield (x0, y1);
            yield (x1, y1);
            e2 = 2_f64 * err;
            if e2 <= dy {
                y0 += 1;
                y1 -= 1;
                dy += a as f64;
                err += dy
            }
            if e2 >= dx || 2_f64 * err > dy {
                x0 += 1;
                x1 -= 1;
                dx += b1 as f64;
                err += dx
            }
            if x0 > x1 {
                break;
            }
        }
        while (y0 - y1) as i64 <= b {
            // I. Quadrant
            // II. Quadrant
            // III. Quadrant
            // IV. Quadrant
            // y step
            // too early stop of flat ellipses a=1
            yield (x0 - 1, y0);
            // -> finish tip of ellipse
            let fresh3 = y0;
            y0 += 1;
            yield (x1 + 1, fresh3);
            yield (x0 - 1, y1);
            let fresh4 = y1;
            y1 -= 1;
            yield (x1 + 1, fresh4);
        }
    })
}

pub fn plot_quad_bezier_seg(
    mut x0: i32,
    mut y0: i32,
    x1: i32,
    mut y1: i32,
    mut x2: i32,
    mut y2: i32,
) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        // plot a limited quadratic Bezier segment
        let mut sx = x2 - x1;
        let mut sy = y2 - y1;
        let mut xx: i64 = (x0 - x1) as i64;
        let mut yy: i64 = (y0 - y1) as i64;
        let mut xy: i64;
        // relative values for checks
        let mut dx: f64;
        let mut dy: f64;
        let mut err: f64;
        let mut cur: f64 = (xx * sy as i64 - yy * sx as i64) as f64;
        // curvature
        assert!(xx * sx as i64 <= 0_i64 && yy * sy as i64 <= 0_i64);

        // sign of gradient must not change
        if sx as i64 * sx as i64 + sy as i64 * sy as i64 > xx * xx + yy * yy {
            // begin with longer part
            x2 = x0;
            x0 = sx + x1;
            y2 = y0;
            y0 = sy + y1;
            cur = -cur
            // swap P0 P2
        }
        if cur != 0 as f64 {
            // no straight line
            xx += sx as i64;
            sx = if x0 < x2 { 1 } else { -(1) };
            xx *= sx as i64;
            // gradient negates -> algorithm fails
            yy += sy as i64;
            sy = if y0 < y2 { 1 } else { -(1) };
            yy *= sy as i64;
            xy = 2_i64 * xx * yy;
            xx *= xx;
            yy *= yy;
            if (cur * sx as f64 * sy as f64) < 0 as f64 {
                // x step direction
                // y step direction
                // differences 2nd degree
                // negated curvature?
                xx = -xx;
                yy = -yy;
                xy = -xy;
                cur = -cur
            }
            dx = 4.0f64 * sy as f64 * cur * (x1 - x0) as f64 + xx as f64 - xy as f64;
            dy = 4.0f64 * sx as f64 * cur * (y0 - y1) as f64 + yy as f64 - xy as f64;
            xx += xx;
            yy += yy;
            err = dx + dy + xy as f64;
            loop
            // differences 1st degree
            // error 1st step
            {
                yield (x0, y0);
                // y step
                if x0 == x2 && y0 == y2 {
                    return;
                }
                y1 = (2_f64 * err < dx) as i32;
                if 2_f64 * err > dy {
                    x0 += sx;
                    dx -= xy as f64;
                    dy += yy as f64;
                    err += dy
                }
                // plot curve
                // last pixel -> curve finished
                // save value for test of y step
                // x step
                if y1 != 0 {
                    y0 += sy;
                    dy -= xy as f64;
                    dx += xx as f64;
                    err += dx
                }
                if !(dy < 0 as f64 && dx > 0 as f64) {
                    break;
                }
            }
        }
        // plot remaining part to end
        for point in plot_line(x0, y0, x2, y2) {
            yield point;
        }
    })
}

pub fn plot_quad_bezier(
    mut x0: i32,
    mut y0: i32,
    mut x1: i32,
    mut y1: i32,
    mut x2: i32,
    mut y2: i32,
) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        // plot any quadratic Bezier curve
        let mut x = x0 - x1;
        let mut y = y0 - y1;
        let mut t: f64 = (x0 - 2 * x1 + x2) as f64;
        let mut r: f64;
        if x as i64 * (x2 - x1) as i64 > 0_i64 {
            // horizontal cut at P4?
            if y as i64 * (y2 - y1) as i64 > 0_i64 {
                // vertical cut at P6 too?
                if ((y0 - 2 * y1 + y2) as f64 / t * x as f64).abs() > (y).abs() as f64 {
                    // which first?
                    x0 = x2;
                    x2 = x + x1;
                    y0 = y2;
                    y2 = y + y1
                    // swap points
                }
            }
            // P0 = P4, P1 = P8
            t = (x0 - x1) as f64 / t;
            r = (1_f64 - t) * ((1_f64 - t) * y0 as f64 + 2.0f64 * t * y1 as f64)
                + t * t * y2 as f64;
            t = (x0 * x2 - x1 * x1) as f64 * t / (x0 - x1) as f64;
            x = (t + 0.5f64).floor() as i32;
            y = (r + 0.5f64).floor() as i32;
            r = (y1 - y0) as f64 * (t - x0 as f64) / (x1 - x0) as f64 + y0 as f64;
            for point in plot_quad_bezier_seg(x0, y0, x, (r + 0.5f64).floor() as i32, x, y) {
                yield point;
            }
            r = (y1 - y2) as f64 * (t - x2 as f64) / (x1 - x2) as f64 + y2 as f64;
            x1 = x;
            x0 = x1;
            y0 = y;
            y1 = (r + 0.5f64).floor() as i32
        }
        if (y0 - y1) as i64 * (y2 - y1) as i64 > 0_i64 {
            // now horizontal cut at P4 comes first
            // By(t=P4)
            // gradient dP4/dx=0
            // intersect P3 | P0 P1
            // intersect P4 | P1 P2
            // vertical cut at P6?
            t = (y0 - 2 * y1 + y2) as f64;
            t = (y0 - y1) as f64 / t;
            r = (1_f64 - t) * ((1_f64 - t) * x0 as f64 + 2.0f64 * t * x1 as f64)
                + t * t * x2 as f64;
            // P0 = P6, P1 = P7
            t = (y0 * y2 - y1 * y1) as f64 * t / (y0 - y1) as f64;
            x = (r + 0.5f64).floor() as i32;
            y = (t + 0.5f64).floor() as i32;
            r = (x1 - x0) as f64 * (t - y0 as f64) / (y1 - y0) as f64 + x0 as f64;
            for point in plot_quad_bezier_seg(x0, y0, (r + 0.5f64).floor() as i32, y, x, y) {
                yield point;
            }
            r = (x1 - x2) as f64 * (t - y2 as f64) / (y1 - y2) as f64 + x2 as f64;
            x0 = x;
            x1 = (r + 0.5f64).floor() as i32;
            y1 = y;
            y0 = y1
        }
        for point in plot_quad_bezier_seg(x0, y0, x1, y1, x2, y2) {
            yield point;
        }
        // Bx(t=P6)
        // gradient dP6/dy=0
        // intersect P6 | P0 P1
        // intersect P7 | P1 P2
        // remaining part
    })
}

pub fn plot_quad_rational_bezier_seg(
    mut x0: i32,
    mut y0: i32,
    mut x1: i32,
    mut y1: i32,
    mut x2: i32,
    mut y2: i32,
    mut w: f32,
) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        // plot a limited rational Bezier segment, squared weight
        let mut sx = x2 - x1;
        let mut sy = y2 - y1;
        // relative values for checks
        let mut dx: f64 = (x0 - x2) as f64;
        let mut dy: f64 = (y0 - y2) as f64;
        let mut xx: f64 = (x0 - x1) as f64;
        let mut yy: f64 = (y0 - y1) as f64;
        let mut xy: f64 = xx * sy as f64 + yy * sx as f64;
        let mut cur: f64 = xx * sy as f64 - yy * sx as f64;
        let mut err: f64;
        // curvature
        assert!(xx * sx as f64 <= 0.0f64 && yy * sy as f64 <= 0.0f64);

        // sign of gradient must not change
        if cur != 0.0f64 && w as f64 > 0.0f64 {
            // no straight line
            if (sx as i64 * sx as i64 + sy as i64 * sy as i64) as f64 > xx * xx + yy * yy {
                // begin with longer part
                x2 = x0;
                x0 = (x0 as f64 - dx) as i32;
                y2 = y0;
                y0 = (y0 as f64 - dy) as i32;
                cur = -cur
                // swap P0 P2
            }
            xx = 2.0f64 * (4.0f64 * w as f64 * sx as f64 * xx + dx * dx);
            // gradient negates -> algorithm fails
            yy = 2.0f64 * (4.0f64 * w as f64 * sy as f64 * yy + dy * dy);
            sx = if x0 < x2 { 1 } else { -(1) };
            sy = if y0 < y2 { 1 } else { -(1) };
            xy = -2.0f64 * sx as f64 * sy as f64 * (2.0f64 * w as f64 * xy + dx * dy);
            if (cur * sx as f64 * sy as f64) < 0.0f64 {
                // differences 2nd degree
                // x step direction
                // y step direction
                // negated curvature?
                xx = -xx;
                yy = -yy;
                xy = -xy;
                cur = -cur
            }
            dx = 4.0f64 * w as f64 * (x1 - x0) as f64 * sy as f64 * cur + xx / 2.0f64 + xy;
            dy = 4.0f64 * w as f64 * (y0 - y1) as f64 * sx as f64 * cur + yy / 2.0f64 + xy;
            if (w as f64) < 0.5f64 && (dy > xy || dx < xy) {
                // differences 1st degree
                // flat ellipse, algorithm fails
                cur = (w as f64 + 1.0f64) / 2.0f64;
                w = (w as f64).sqrt() as f32;
                xy = 1.0f64 / (w as f64 + 1.0f64);
                sx = ((x0 as f64 + 2.0f64 * w as f64 * x1 as f64 + x2 as f64) * xy / 2.0f64 + 0.5f64)
                    .floor() as i32;
                // subdivide curve in half
                sy = ((y0 as f64 + 2.0f64 * w as f64 * y1 as f64 + y2 as f64) * xy / 2.0f64 + 0.5f64)
                    .floor() as i32;
                dx = ((w * x1 as f32 + x0 as f32) as f64 * xy + 0.5f64).floor();
                dy = ((y1 as f32 * w + y0 as f32) as f64 * xy + 0.5f64).floor();
                let next_points: Box<dyn Iterator<Item = Pixel>> = Box::new(plot_quad_rational_bezier_seg(x0, y0, dx as i32, dy as i32, sx, sy, cur as f32));
                for point in next_points {
                    yield point;
                }
                // plot separately
                dx = ((w * x1 as f32 + x2 as f32) as f64 * xy + 0.5f64).floor();
                dy = ((y1 as f32 * w + y2 as f32) as f64 * xy + 0.5f64).floor();
                let next_points: Box<dyn Iterator<Item = Pixel>> = Box::new(plot_quad_rational_bezier_seg(sx, sy, dx as i32, dy as i32, x2, y2, cur as f32));
                for point in next_points {
                    yield point;
                }
                return;
            }
            // error 1.step
            err = dx + dy - xy;
            loop {
                yield (x0, y0);
                // x step
                if x0 == x2 && y0 == y2 {
                    return;
                }
                x1 = (2_f64 * err > dy) as i32;
                y1 = (2_f64 * (err + yy) < -dy) as i32;
                if 2_f64 * err < dx || y1 != 0 {
                    y0 += sy;
                    dy += xy;
                    dx += xx;
                    err += dx
                }
                // plot curve
                // last pixel -> curve finished
                // save value for test of x step
                // y step
                if 2_f64 * err > dx || x1 != 0 {
                    x0 += sx;
                    dx += xy;
                    dy += yy;
                    err += dy
                }
                if !(dy <= xy && dx >= xy) {
                    break;
                }
            }
        }
        // plot remaining needle to end
        for point in plot_line(x0, y0, x2, y2) {
            yield point;
        }
    })
}

pub fn plot_quad_rational_bezier(
    mut x0: i32,
    mut y0: i32,
    mut x1: i32,
    mut y1: i32,
    mut x2: i32,
    mut y2: i32,
    mut w: f32,
) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        // plot any quadratic rational Bezier curve
        let mut x = x0 - 2 * x1 + x2;
        let mut y = y0 - 2 * y1 + y2;
        let mut xx: f64 = (x0 - x1) as f64;
        let mut yy: f64 = (y0 - y1) as f64;
        let mut ww: f64;
        let mut t: f64;
        let mut q: f64;
        assert!(w as f64 >= 0.0f64);
        if xx * (x2 - x1) as f64 > 0 as f64 {
            // horizontal cut at P4?
            if yy * (y2 - y1) as f64 > 0 as f64 {
                // vertical cut at P6 too?
                if (xx * y as f64).abs() > (yy * x as f64).abs() {
                    // which first?
                    x0 = x2;
                    x2 = (xx + x1 as f64) as i32;
                    y0 = y2;
                    y2 = (yy + y1 as f64) as i32;
                    // swap points
                }
            }
            // P0 = P4, P1 = P8
            if x0 == x2 || w as f64 == 1.0f64 {
                t = (x0 - x1) as f64 / x as f64
            } else {
                // now horizontal cut at P4 comes first
                // non-rational or rational case
                q = (4.0f64 * w as f64 * w as f64 * (x0 - x1) as f64 * (x2 - x1) as f64
                    + ((x2 - x0) as i64 * (x2 - x0) as i64) as f64)
                    .sqrt();
                if x1 < x0 {
                    q = -q
                }
                t = (2.0f64 * w as f64 * (x0 - x1) as f64 - x0 as f64 + x2 as f64 + q)
                    / (2.0f64 * (1.0f64 - w as f64) * (x2 - x0) as f64)
                // t at P4
            }
            q = 1.0f64 / (2.0f64 * t * (1.0f64 - t) * (w as f64 - 1.0f64) + 1.0f64);
            xx = (t * t * (x0 as f64 - 2.0f64 * w as f64 * x1 as f64 + x2 as f64)
                + 2.0f64 * t * (w * x1 as f32 - x0 as f32) as f64
                + x0 as f64)
                * q;
            yy = (t * t * (y0 as f64 - 2.0f64 * w as f64 * y1 as f64 + y2 as f64)
                + 2.0f64 * t * (w * y1 as f32 - y0 as f32) as f64
                + y0 as f64)
                * q;
            ww = t * (w as f64 - 1.0f64) + 1.0f64;
            ww = ww * ww * q;
            w = (((1.0f64 - t) * (w as f64 - 1.0f64) + 1.0f64) * (q).sqrt()) as f32;
            x = (xx + 0.5f64).floor() as i32;
            y = (yy + 0.5f64).floor() as i32;
            yy = (xx - x0 as f64) * (y1 - y0) as f64 / (x1 - x0) as f64 + y0 as f64;
            for point in plot_quad_rational_bezier_seg(x0, y0, x, (yy + 0.5f64).floor() as i32, x, y, ww as f32) {
                yield point;
            }
            yy = (xx - x2 as f64) * (y1 - y2) as f64 / (x1 - x2) as f64 + y2 as f64;
            y1 = (yy + 0.5f64).floor() as i32;
            x1 = x;
            x0 = x1;
            y0 = y
        }
        if (y0 - y1) as i64 * (y2 - y1) as i64 > 0_i64 {
            // sub-divide at t
            // = P4
            // squared weight P3
            // weight P8
            // P4
            // intersect P3 | P0 P1
            // intersect P4 | P1 P2
            // vertical cut at P6?
            if y0 == y2 || w as f64 == 1.0f64 {
                t = (y0 - y1) as f64 / (y0 as f64 - 2.0f64 * y1 as f64 + y2 as f64)
            } else {
                // non-rational or rational case
                q = (4.0f64 * w as f64 * w as f64 * (y0 - y1) as f64 * (y2 - y1) as f64
                    + ((y2 - y0) as i64 * (y2 - y0) as i64) as f64)
                    .sqrt();
                if y1 < y0 {
                    q = -q
                }
                t = (2.0f64 * w as f64 * (y0 - y1) as f64 - y0 as f64 + y2 as f64 + q)
                    / (2.0f64 * (1.0f64 - w as f64) * (y2 - y0) as f64)
                // t at P6
            }
            q = 1.0f64 / (2.0f64 * t * (1.0f64 - t) * (w as f64 - 1.0f64) + 1.0f64);
            // P0 = P6, P1 = P7
            xx = (t * t * (x0 as f64 - 2.0f64 * w as f64 * x1 as f64 + x2 as f64)
                + 2.0f64 * t * (w * x1 as f32 - x0 as f32) as f64
                + x0 as f64)
                * q;
            yy = (t * t * (y0 as f64 - 2.0f64 * w as f64 * y1 as f64 + y2 as f64)
                + 2.0f64 * t * (w * y1 as f32 - y0 as f32) as f64
                + y0 as f64)
                * q;
            ww = t * (w as f64 - 1.0f64) + 1.0f64;
            ww = ww * ww * q;
            w = (((1.0f64 - t) * (w as f64 - 1.0f64) + 1.0f64) * (q).sqrt()) as f32;
            x = (xx + 0.5f64).floor() as i32;
            y = (yy + 0.5f64).floor() as i32;
            xx = (x1 - x0) as f64 * (yy - y0 as f64) / (y1 - y0) as f64 + x0 as f64;
            for point in plot_quad_rational_bezier_seg(x0, y0, (xx + 0.5f64).floor() as i32, y, x, y, ww as f32) {
                yield point;
            }
            xx = (x1 - x2) as f64 * (yy - y2 as f64) / (y1 - y2) as f64 + x2 as f64;
            x1 = (xx + 0.5f64).floor() as i32;
            x0 = x;
            y1 = y;
            y0 = y1
        }
        for point in plot_quad_rational_bezier_seg(x0, y0, x1, y1, x2, y2, w * w) {
            yield point;
        }
        // sub-divide at t
        // = P6
        // squared weight P5
        // weight P7
        // P6
        // intersect P6 | P0 P1
        // intersect P7 | P1 P2
        // remaining
    })
}

pub fn plot_rotated_ellipse_rect(
    x0: i32,
    y0: i32,
    x1: i32,
    y1: i32,
    zd: i64,
) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        // rectangle enclosing the ellipse, integer rotation angle
        let mut xd = x1 - x0;
        let mut yd = y1 - y0;
        let mut w: f32 = (xd as i64 * yd as i64) as f32;
        if zd == 0_i64 {
            for point in plot_ellipse_rect(x0, y0, x1, y1) {
                yield point;
            }
            return;
        }
        // looks nicer
        if w as f64 != 0.0f64 {
            w = (w - zd as f32) / (w + w)
        }
        // squared weight of P1
        assert!(w as f64 <= 1.0f64 && w as f64 >= 0.0f64);
        xd = ((xd as f32 * w).floor() as f64 + 0.5f64) as i32;
        yd = ((yd as f32 * w).floor() as f64 + 0.5f64) as i32;
        // snap xe,ye to int
        for point in plot_quad_rational_bezier_seg(x0, y0 + yd, x0, y0, x0 + xd, y0, (1.0f64 - w as f64) as f32) {
            yield point;
        }
        for point in plot_quad_rational_bezier_seg(x0, y0 + yd, x0, y1, x1 - xd, y1, w) {
            yield point;
        }
        for point in plot_quad_rational_bezier_seg(x1, y1 - yd, x1, y1, x1 - xd, y1, (1.0f64 - w as f64) as f32) {
            yield point;
        }
        for point in plot_quad_rational_bezier_seg(x1, y1 - yd, x1, y0, x0 + xd, y0, w) {
            yield point;
        }
    })
}

pub fn plot_rotated_ellipse(
    x: i32,
    y: i32,
    mut a: i32,
    mut b: i32,
    angle: f32,
) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        // plot ellipse rotated by angle (radian)
        let mut xd: f32 = (a as i64 * a as i64) as f32;
        let mut yd: f32 = (b as i64 * b as i64) as f32;
        let s: f32 = (angle as f64).sin() as f32;
        let mut zd: f32 = (xd - yd) * s;
        // ellipse rotation
        xd = ((xd - zd * s).sqrt() as f64) as f32;
        yd = ((yd + zd * s).sqrt() as f64) as f32;
        // surrounding rectangle
        a = (xd as f64 + 0.5f64) as i32;
        b = (yd as f64 + 0.5f64) as i32;
        zd = zd * a as f32 * b as f32 / (xd * yd);
        // scale to integer
        for point in plot_rotated_ellipse_rect(
            x - a,
            y - b,
            x + a,
            y + b,
            ((4_f32 * zd) as f64 * (angle as f64).cos()) as i64,
        ) {
            yield point;
        }
    })
}

pub fn plot_cubic_bezier_seg(
    mut x0: i32,
    mut y0: i32,
    mut x1: f32,
    mut y1: f32,
    mut x2: f32,
    y2: f32,
    mut x3: i32,
    mut y3: i32,
) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        // plot limited cubic Bezier segment
        let mut f;
        let mut fx;
        let mut fy;
        let mut leg = 1;
        let mut sx = if x0 < x3 { 1 } else { -(1) };
        let mut sy = if y0 < y3 { 1 } else { -(1) };
        // step direction
        let xc: f32 = -((x0 as f32 + x1 - x2 - x3 as f32) as f64).abs() as f32;
        let xa: f32 = xc - (4 * sx) as f32 * (x1 - x2);
        let mut xb: f32 = sx as f32 * (x0 as f32 - x1 - x2 + x3 as f32);
        let yc: f32 = -((y0 as f32 + y1 - y2 - y3 as f32) as f64).abs() as f32;
        let ya: f32 = yc - (4 * sy) as f32 * (y1 - y2);
        let mut yb: f32 = sy as f32 * (y0 as f32 - y1 - y2 + y3 as f32);
        let mut ab: f64;
        let mut ac: f64;
        let mut bc: f64;
        let mut cb: f64;
        let mut xx: f64;
        let mut xy: f64;
        let mut yy: f64;
        let mut dx: f64;
        let mut dy: f64;
        let mut ex: f64;
        let mut pxy: f64;
        let ep: f64 = 0.01f64;
        // check for curve restrains
        // slope P0-P1 == P2-P3    and  (P0-P3 == P1-P2      or   no slope change)
        assert!(
            (((x1 - x0 as f32) * (x2 - x3 as f32)) as f64) < ep
                && ((((x3 - x0) as f32 * (x1 - x2)) as f64) < ep
                    || ((xb * xb) as f64) < (xa * xc) as f64 + ep)
        );
        assert!(
            (((y1 - y0 as f32) * (y2 - y3 as f32)) as f64) < ep
                && ((((y3 - y0) as f32 * (y1 - y2)) as f64) < ep
                    || ((yb * yb) as f64) < (ya * yc) as f64 + ep)
        );
        if xa == 0 as f32 && ya == 0 as f32 {
            // quadratic Bezier
            sx = (((3_f32 * x1 - x0 as f32 + 1_f32) / 2_f32) as f64).floor() as i32;
            sy = (((3_f32 * y1 - y0 as f32 + 1_f32) / 2_f32) as f64).floor() as i32;
            // new midpoint
            for point in plot_quad_bezier_seg(x0, y0, sx, sy, x3, y3) {
                yield point;
            }
            return;
        }
        x1 = (x1 - x0 as f32) * (x1 - x0 as f32) + (y1 - y0 as f32) * (y1 - y0 as f32) + 1_f32;
        // line lengths
        x2 = (x2 - x3 as f32) * (x2 - x3 as f32) + (y2 - y3 as f32) * (y2 - y3 as f32) + 1_f32;
        loop {
            // loop over both ends
            ab = (xa * yb - xb * ya) as f64;
            ac = (xa * yc - xc * ya) as f64;
            bc = (xb * yc - xc * yb) as f64;
            ex = ab * (ab + ac - 3_f64 * bc) + ac * ac;
            // P0 part of self-intersection loop?
            f = if ex > 0 as f64 {
                1.0 as i32
            } else {
                ((1_f32 + 1024_f32 / x1) as f64).sqrt() as i32
            };
            // calculate resolution
            ab *= f as f64;
            ac *= f as f64;
            bc *= f as f64;
            ex *= (f * f) as f64;
            // increase resolution
            xy = 9_f64 * (ab + ac + bc) / 8_f64;
            cb = (8_f32 * (xa - ya)) as f64;
            // init differences of 1st degree
            dx = 27_f64
                * (8_f64 * ab * (yb * yb - ya * yc) as f64 + ex * (ya + 2_f32 * yb + yc) as f64)
                / 64_f64
                - (ya * ya) as f64 * (xy - ya as f64);
            dy = 27_f64
                * (8_f64 * ab * (xb * xb - xa * xc) as f64 - ex * (xa + 2_f32 * xb + xc) as f64)
                / 64_f64
                - (xa * xa) as f64 * (xy + xa as f64);
            // init differences of 2nd degree
            xx = 3_f64
                * (3_f64 * ab * (3_f32 * yb * yb - ya * ya - 2_f32 * ya * yc) as f64
                    - ya as f64 * (3_f64 * ac * (ya + yb) as f64 + ya as f64 * cb))
                / 4_f64;
            yy = 3_f64
                * (3_f64 * ab * (3_f32 * xb * xb - xa * xa - 2_f32 * xa * xc) as f64
                    - xa as f64 * (3_f64 * ac * (xa + xb) as f64 + xa as f64 * cb))
                / 4_f64;
            xy = (xa * ya) as f64 * (6_f64 * ab + 6_f64 * ac - 3_f64 * bc + cb);
            ac = (ya * ya) as f64;
            cb = (xa * xa) as f64;
            xy = 3_f64
                * (xy + (9 * f) as f64 * (cb * yb as f64 * yc as f64 - (xb * xc) as f64 * ac)
                    - (18_f32 * xb * yb) as f64 * ab)
                / 8_f64;
            if ex < 0 as f64 {
                // negate values if inside self-intersection loop
                dx = -dx;
                dy = -dy;
                xx = -xx;
                yy = -yy;
                xy = -xy;
                ac = -ac;
                cb = -cb
            }
            // init differences of 3rd degree
            ab = (6_f32 * ya) as f64 * ac;
            ac *= (-(6) as f32 * xa) as f64;
            bc = (6_f32 * ya) as f64 * cb;
            cb *= (-(6) as f32 * xa) as f64;
            dx += xy;
            ex = dx + dy;
            dy += xy;
            // error of 1st step
            pxy = xy;
            fy = f;
            fx = fy;
            's_201: while x0 != x3 && y0 != y3 {
                yield (x0, y0);
                // pixel ahead valid
                // plot curve
                loop
                // move sub-steps of one pixel
                {
                    if dx > pxy || dy < pxy {
                        break 's_201;
                    }
                    // confusing values
                    y1 = (2_f64 * ex - dy) as f32;
                    // save value for test of y step
                    if 2_f64 * ex >= dx {
                        // x sub-step
                        fx -= 1;
                        dx += xx;
                        ex += dx;
                        xy += ac;
                        dy += xy;
                        yy += bc;
                        xx += ab
                    }
                    if y1 <= 0 as f32 {
                        // y sub-step
                        fy -= 1;
                        dy += yy;
                        ex += dy;
                        xy += bc;
                        dx += xy;
                        xx += ac;
                        yy += cb
                    }
                    if !(fx > 0 && fy > 0) {
                        break;
                    }
                }
                // pixel complete?
                if 2 * fx <= f {
                    x0 += sx;
                    fx += f
                }
                // x step
                if 2 * fy <= f {
                    y0 += sy;
                    fy += f
                }
                if pxy == xy as f64 && dx < 0 as f64 && dy > 0 as f64 {
                    pxy = ep;
                }
            }
            xx = x0 as f64;
            x0 = x3;
            x3 = xx as i32;
            sx = -sx;
            xb = -xb;
            // y step
            // swap legs
            yy = y0 as f64;
            y0 = y3;
            y3 = yy as i32;
            sy = -sy;
            yb = -yb;
            x1 = x2;
            let fresh5 = leg;
            leg -= 1;
            if fresh5 == 0 {
                break;
            }
        }
        // try other end
        // remaining part in case of cusp or crunode
        for point in plot_line(x0, y0, x3, y3) {
            yield point;
        }
    })
}

pub fn plot_cubic_bezier(
    mut x0: i32,
    mut y0: i32,
    x1: i32,
    y1: i32,
    x2: i32,
    y2: i32,
    mut x3: i32,
    mut y3: i32,
) -> impl Iterator<Item = Pixel> {
    gen_iter!(move {
        // plot any cubic Bezier curve
        let mut n = 0;
        let mut i;
        let xc: i64 = (x0 + x1 - x2 - x3) as i64;
        let xa: i64 = xc - (4 * (x1 - x2)) as i64;
        let xb: i64 = (x0 - x1 - x2 + x3) as i64;
        let xd: i64 = xb + (4 * (x1 + x2)) as i64;
        let yc: i64 = (y0 + y1 - y2 - y3) as i64;
        let ya: i64 = yc - (4 * (y1 - y2)) as i64;
        let yb: i64 = (y0 - y1 - y2 + y3) as i64;
        let yd: i64 = yb + (4 * (y1 + y2)) as i64;
        let mut fx0: f32 = x0 as f32;
        let mut fx1: f32;
        let mut fx2: f32;
        let mut fx3: f32;
        let mut fy0: f32 = y0 as f32;
        let mut fy1: f32;
        let mut fy2: f32;
        let mut fy3: f32;
        let mut t1: f64 = (xb * xb - xa * xc) as f64;
        let mut t2: f64;
        let mut t: [f64; 5] = [0.; 5];
        // sub-divide curve at gradient sign changes
        if xa == 0_i64 {
            // horizontal
            if (xc).abs() < 2 * (xb).abs() {
                let fresh6 = n;
                n += 1;
                t[fresh6 as usize] = xc as f64 / (2.0f64 * xb as f64)
            }
            // one change
        } else if t1 > 0.0f64 {
            // two changes
            t2 = (t1).sqrt();
            t1 = (xb as f64 - t2) / xa as f64;
            if (t1).abs() < 1.0f64 {
                let fresh7 = n;
                n += 1;
                t[fresh7 as usize] = t1
            }
            t1 = (xb as f64 + t2) / xa as f64;
            if (t1).abs() < 1.0f64 {
                let fresh8 = n;
                n += 1;
                t[fresh8 as usize] = t1
            }
        }
        t1 = (yb * yb - ya * yc) as f64;
        if ya == 0_i64 {
            // vertical
            if (yc).abs() < 2 * (yb).abs() {
                let fresh9 = n;
                n += 1;
                t[fresh9 as usize] = yc as f64 / (2.0f64 * yb as f64)
            }
            // one change
        } else if t1 > 0.0f64 {
            // two changes
            t2 = (t1).sqrt();
            t1 = (yb as f64 - t2) / ya as f64;
            if (t1).abs() < 1.0f64 {
                let fresh10 = n;
                n += 1;
                t[fresh10 as usize] = t1
            }
            t1 = (yb as f64 + t2) / ya as f64;
            if (t1).abs() < 1.0f64 {
                let fresh11 = n;
                n += 1;
                t[fresh11] = t1
            }
        }
        i = 1;
        while i < n {
            // bubble sort of 4 points
            t1 = t[(i - 1) as usize];
            if t1 > t[i] {
                t[(i - 1) as usize] = t[i];
                t[i] = t1;
                i = 0
            }
            i += 1
        }
        t1 = -1.0f64;
        t[n as usize] = 1.0f64;
        // begin / end point
        i = 0;
        while i <= n {
            // plot each segment separately
            t2 = t[i];
            // sub-divide at t[i-1], t[i]
            fx1 = ((t1 * (t1 * xb as f64 - (2_i64 * xc) as f64)
                - t2 * (t1 * (t1 * xa as f64 - (2_i64 * xb) as f64) + xc as f64)
                + xd as f64)
                / 8_f64
                - fx0 as f64) as f32;
            fy1 = ((t1 * (t1 * yb as f64 - (2_i64 * yc) as f64)
                - t2 * (t1 * (t1 * ya as f64 - (2_i64 * yb) as f64) + yc as f64)
                + yd as f64)
                / 8_f64
                - fy0 as f64) as f32;
            fx2 = ((t2 * (t2 * xb as f64 - (2_i64 * xc) as f64)
                - t1 * (t2 * (t2 * xa as f64 - (2_i64 * xb) as f64) + xc as f64)
                + xd as f64)
                / 8_f64
                - fx0 as f64) as f32;
            fy2 = ((t2 * (t2 * yb as f64 - (2_i64 * yc) as f64)
                - t1 * (t2 * (t2 * ya as f64 - (2_i64 * yb) as f64) + yc as f64)
                + yd as f64)
                / 8_f64
                - fy0 as f64) as f32;
            fx3 = ((t2 * (t2 * ((3_i64 * xb) as f64 - t2 * xa as f64) - (3_i64 * xc) as f64)
                + xd as f64)
                / 8_f64) as f32;
            fx0 -= fx3;
            fy3 = ((t2 * (t2 * ((3_i64 * yb) as f64 - t2 * ya as f64) - (3_i64 * yc) as f64)
                + yd as f64)
                / 8_f64) as f32;
            fy0 -= fy3;
            x3 = (fx3 as f64 + 0.5f64).floor() as i32;
            y3 = (fy3 as f64 + 0.5f64).floor() as i32;
            // scale bounds to int
            if fx0 as f64 != 0.0f64 {
                fx0 = (x0 - x3) as f32 / fx0;
                fx1 *= fx0;
                fx2 *= fx0
            }
            if fy0 as f64 != 0.0f64 {
                fy0 = (y0 - y3) as f32 / fy0;
                fy1 *= fy0;
                fy2 *= fy0
            }
            if x0 != x3 || y0 != y3 {
                // segment t1 - t2
                for point in plot_cubic_bezier_seg(
                    x0,
                    y0,
                    x0 as f32 + fx1,
                    y0 as f32 + fy1,
                    x0 as f32 + fx2,
                    y0 as f32 + fy2,
                    x3,
                    y3,
                ) {
                    yield point;
                }
            }
            x0 = x3;
            y0 = y3;
            fx0 = fx3;
            fy0 = fy3;
            t1 = t2;
            i += 1
        }
    })
}

pub fn plot_line_aa(mut x0: i32, mut y0: i32, x1: i32, y1: i32) -> impl Iterator<Item = PixelAA> {
    gen_iter!(move {
        // draw a black (0) anti-aliased line on white (255) background
        let sx = if x0 < x1 { 1 } else { -(1) };
        let sy = if y0 < y1 { 1 } else { -(1) };
        let mut x2;
        let mut dx: i64 = (x1 - x0).abs() as i64;
        let mut dy: i64 = (y1 - y0).abs() as i64;
        let mut err: i64 = dx * dx + dy * dy;
        let mut e2: i64 = if err == 0_i64 {
            1_f64
        } else {
            (0xffff7f_i64 as f64) / (err as f64).sqrt()
        } as i64;
        // multiplication factor
        dx *= e2;
        dy *= e2;
        // error value e_xy
        err = dx - dy;
        // pixel loop
        loop {
            yield (x0, y0, (((err - dx + dy).abs()) >> 16) as u8);
            e2 = err;
            x2 = x0;
            if 2_i64 * e2 >= -dx {
                // x step
                if x0 == x1 {
                    break;
                }
                if e2 + dy < 0xff0000_i64 {
                    yield (x0, y0 + sy, ((e2 + dy) >> 16) as u8);
                }
                err -= dy;
                x0 += sx
            }
            if 2_i64 * e2 > dy {
                continue;
            }
            // y step
            if y0 == y1 {
                break;
            }
            if dx - e2 < 0xff0000_i64 {
                yield (x2 + sx, y0, ((dx - e2) >> 16) as u8);
            }
            err += dx;
            y0 += sy
        }
    })
}

pub fn plot_circle_aa(xm: i32, ym: i32, mut r: i32) -> impl Iterator<Item = PixelAA> {
    gen_iter!(move {
        // draw a black anti-aliased circle on white background
        let mut x = -r;
        let mut y = 0;
        // II. quadrant from bottom left to top right
        let mut i;
        let mut x2;
        let mut e2;
        let mut err = 2 - 2 * r;
        // error of 1.step
        r = 1 - err;
        loop {
            i = 255 * (err - 2 * (x + y).abs() - 2) / r;
            // get blend value of pixel
            yield (xm - x, ym + y, i as u8);
            // I. Quadrant
            yield (xm - y, ym - x, i as u8);
            // II. Quadrant
            yield (xm + x, ym - y, i as u8);
            // III. Quadrant
            yield (xm + y, ym + x, i as u8);
            // IV. Quadrant
            e2 = err;
            x2 = x;
            // remember values
            if err + y > 0 {
                // x step
                i = 255 * (err - 2 * x - 1) / r;
                // outward pixel
                if i < 256 {
                    yield (xm - x, ym + y + 1, i as u8);
                    yield (xm - y - 1, ym - x, i as u8);
                    yield (xm + x, ym - y - 1, i as u8);
                    yield (xm + y + 1, ym + x, i as u8);
                }
                x += 1;
                err += x * 2 + 1
            }
            if e2 + x2 <= 0 {
                // y step
                i = 255 * (2 * y + 3 - e2) / r;
                // inward pixel
                if i < 256 {
                    yield (xm - x2 - 1, ym + y, i as u8);
                    yield (xm - y, ym - x2 - 1, i as u8);
                    yield (xm + x2 + 1, ym - y, i as u8);
                    yield (xm + y, ym + x2 + 1, i as u8);
                }
                y += 1;
                err += y * 2 + 1
            }
            if x >= 0 {
                break;
            }
        }
    })
}

pub fn plot_ellipse_rect_aa(
    mut x0: i32,
    mut y0: i32,
    mut x1: i32,
    mut y1: i32,
) -> impl Iterator<Item = PixelAA> {
    gen_iter!(move {
        // draw a black anti-aliased rectangular ellipse on white background
        let mut a: i64 = (x1 - x0).abs() as i64;
        let b: i64 = (y1 - y0).abs() as i64;
        let mut b1: i64 = b & 1_i64;
        // diameter
        let mut dx: f32 = (4_f64 * (a as f64 - 1.0f64) * b as f64 * b as f64) as f32;
        let mut dy: f32 = (4_i64 * (b1 + 1_i64) * a * a) as f32;
        // error increment
        let mut ed: f32;
        let mut i: f32;
        let mut err: f32 = (b1 * a * a) as f32 - dx + dy;
        // error of 1.step
        let mut f: bool;
        if a == 0_i64 || b == 0_i64 {
            for point in plot_line(x0, y0, x1, y1) {
                yield (point.0, point.1, 255);
            }
            return;
        }
        if x0 > x1 {
            x0 = x1;
            x1 = (x1 as i64 + a) as i32
        }
        // if called with swapped points
        if y0 > y1 {
            y0 = y1
        }
        // .. exchange them
        y0 = (y0 as i64 + (b + 1_i64) / 2_i64) as i32;
        y1 = (y0 as i64 - b1) as i32;
        // starting pixel
        a = 8_i64 * a * a;
        b1 = 8_i64 * b * b;
        loop {
            /* approximate ed=(dx*dx+dy*dy).sqrt() */
            i = (dx as f64).min(dy as f64) as f32;
            ed = (dx as f64).max(dy as f64) as f32;
            if y0 == y1 + 1 && err > dy && a > b1 {
                ed = (255_f64 * 4.0f64 / a as f64) as f32
            } else {
                // x error increment
                // x-tip
                ed = 255_f32 / (ed + 2_f32 * ed * i * i / (4_f32 * ed * ed + i * i))
            }
            // approximation
            i = (ed as f64 * ((err + dx - dy) as f64).abs()) as f32;
            // get intensity value by pixel error
            yield (x0, y0, i as u8);
            yield (x0, y1, i as u8);
            yield (x1, y0, i as u8);
            yield (x1, y1, i as u8);
            f = 2_f32 * err + dy >= 0 as f32;
            if f {
                // x step, remember condition
                if x0 >= x1 {
                    break;
                }
                i = ed * (err + dx);
                if i < 255_f32 {
                    yield (x0, y0 + 1, i as u8);
                    yield (x0, y1 - 1, i as u8);
                    yield (x1, y0 + 1, i as u8);
                    yield (x1, y1 - 1, i as u8);
                }
                // do error increment later since values are still needed
            }
            if 2_f32 * err <= dx {
                // y step
                i = ed * (dy - err);
                if i < 255_f32 {
                    yield (x0 + 1, y0, i as u8);
                    yield (x1 - 1, y0, i as u8);
                    yield (x0 + 1, y1, i as u8);
                    yield (x1 - 1, y1, i as u8);
                }
                y0 += 1;
                y1 -= 1;
                dy += a as f32;
                err += dy
            }
            if f {
                x0 += 1;
                x1 -= 1;
                dx -= b1 as f32;
                err -= dx
            }
        }
        x0 -= 1;
        let fresh12 = x1;
        x1 += 1;
        if x0 == fresh12 {
            // too early stop of flat ellipses
            while ((y0 - y1) as i64) < b {
                i = ((255 * 4) as f64 * ((err + dx) as f64).abs() / b1 as f64) as f32;
                // -> finish tip of ellipse
                y0 += 1;
                yield (x0, y0, i as u8);
                yield (x1, y0, i as u8);
                y1 -= 1;
                yield (x0, y1, i as u8);
                yield (x1, y1, i as u8);
                dy += a as f32;
                err += dy
            }
        };
    })
}

pub fn plot_quad_bezier_seg_aa(
    mut x0: i32,
    mut y0: i32,
    mut x1: i32,
    y1: i32,
    mut x2: i32,
    mut y2: i32,
) -> impl Iterator<Item = PixelAA> {
    gen_iter!(move {
        // draw an limited anti-aliased quadratic Bezier segment
        let mut sx = x2 - x1;
        let mut sy = y2 - y1;
        let mut xx: i64 = (x0 - x1) as i64;
        let mut yy: i64 = (y0 - y1) as i64;
        let mut xy: i64;
        // relative values for checks
        let mut dx: f64;
        let mut dy: f64;
        let mut err: f64;
        let mut ed: f64;
        let mut cur: f64 = (xx * sy as i64 - yy * sx as i64) as f64;
        // curvature
        assert!(xx * sx as i64 <= 0_i64 && yy * sy as i64 <= 0_i64);

        // sign of gradient must not change
        if sx as i64 * sx as i64 + sy as i64 * sy as i64 > xx * xx + yy * yy {
            // begin with longer part
            x2 = x0;
            x0 = sx + x1;
            y2 = y0;
            y0 = sy + y1;
            cur = -cur
            // swap P0 P2
        }
        if cur != 0 as f64 {
            // no straight line
            xx += sx as i64;
            sx = if x0 < x2 { 1 } else { -(1) };
            xx *= sx as i64;
            // gradient negates -> close curves
            yy += sy as i64;
            sy = if y0 < y2 { 1 } else { -(1) };
            yy *= sy as i64;
            xy = 2_i64 * xx * yy;
            xx *= xx;
            yy *= yy;
            if (cur * sx as f64 * sy as f64) < 0 as f64 {
                // x step direction
                // y step direction
                // differences 2nd degree
                // negated curvature?
                xx = -xx;
                yy = -yy;
                xy = -xy;
                cur = -cur
            }
            dx = 4.0f64 * sy as f64 * (x1 - x0) as f64 * cur + xx as f64 - xy as f64;
            dy = 4.0f64 * sx as f64 * (y0 - y1) as f64 * cur + yy as f64 - xy as f64;
            xx += xx;
            yy += yy;
            err = dx + dy + xy as f64;
            loop
            // differences 1st degree
            // error 1st step
            {
                cur = (dx + xy as f64).min(-xy as f64 - dy);
                ed = (dx + xy as f64).max(-xy as f64 - dy);
                // approximate error distance
                ed += 2_f64 * ed * cur * cur / (4_f64 * ed * ed + cur * cur);
                yield (x0, y0, (255_f64 * (err - dx - dy - xy as f64).abs() / ed) as u8);
                // plot curve
                if x0 == x2 || y0 == y2 {
                    break;
                }
                // last pixel -> curve finished
                x1 = x0;
                cur = dx - err;
                let y1 = 2_f64 * err + dy < 0 as f64;
                if 2_f64 * err + dx > 0 as f64 {
                    // x step
                    if err - dy < ed {
                        yield (x0, y0 + sy, (255_f64 * (err - dy).abs() / ed) as u8);
                    }
                    x0 += sx;
                    dx -= xy as f64;
                    dy += yy as f64;
                    err += dy
                }
                if y1 {
                    // y step
                    if cur < ed {
                        yield (x1 + sx, y0, (255_f64 * (cur).abs() / ed) as u8);
                    }
                    y0 += sy;
                    dy -= xy as f64;
                    dx += xx as f64;
                    err += dx
                }
                if dy >= dx {
                    break;
                }
            }
        }
        for point in plot_line_aa(x0, y0, x2, y2) {
            yield point;
        }
        // plot remaining needle to end
    })
}

pub fn plot_quad_rational_bezier_seg_aa(
    mut x0: i32,
    mut y0: i32,
    mut x1: i32,
    y1: i32,
    mut x2: i32,
    mut y2: i32,
    mut w: f32,
) -> impl Iterator<Item = PixelAA> {
    gen_iter!(move {
        // draw an anti-aliased rational quadratic Bezier segment, squared weight
        let mut sx = x2 - x1;
        let mut sy = y2 - y1;
        // relative values for checks
        let mut dx: f64 = (x0 - x2) as f64;
        let mut dy: f64 = (y0 - y2) as f64;
        let mut xx: f64 = (x0 - x1) as f64;
        let mut yy: f64 = (y0 - y1) as f64;
        let mut xy: f64 = xx * sy as f64 + yy * sx as f64;
        let mut cur: f64 = xx * sy as f64 - yy * sx as f64;
        let mut err: f64;
        let mut ed: f64;
        // curvature
        let mut f: bool;
        assert!(xx * sx as f64 <= 0.0f64 && yy * sy as f64 <= 0.0f64);
        // sign of gradient must not change
        if cur != 0.0f64 && w as f64 > 0.0f64 {
            // no straight line
            if (sx as i64 * sx as i64 + sy as i64 * sy as i64) as f64 > xx * xx + yy * yy {
                // begin with longer part
                x2 = x0;
                x0 = (x0 as f64 - dx) as i32;
                y2 = y0;
                y0 = (y0 as f64 - dy) as i32;
                cur = -cur
                // swap P0 P2
            }
            xx = 2.0f64 * (4.0f64 * w as f64 * sx as f64 * xx + dx * dx);
            // gradient negates -> algorithm fails
            yy = 2.0f64 * (4.0f64 * w as f64 * sy as f64 * yy + dy * dy);
            sx = if x0 < x2 { 1 } else { -(1) };
            sy = if y0 < y2 { 1 } else { -(1) };
            xy = -2.0f64 * sx as f64 * sy as f64 * (2.0f64 * w as f64 * xy + dx * dy);
            if (cur * sx as f64 * sy as f64) < 0 as f64 {
                // differences 2nd degree
                // x step direction
                // y step direction
                // negated curvature?
                xx = -xx;
                yy = -yy;
                cur = -cur;
                xy = -xy
            }
            dx = 4.0f64 * w as f64 * (x1 - x0) as f64 * sy as f64 * cur + xx / 2.0f64 + xy;
            dy = 4.0f64 * w as f64 * (y0 - y1) as f64 * sx as f64 * cur + yy / 2.0f64 + xy;
            if (w as f64) < 0.5f64 && dy > dx {
                // differences 1st degree
                // flat ellipse, algorithm fails
                cur = (w as f64 + 1.0f64) / 2.0f64;
                w = (w as f64).sqrt() as f32;
                xy = 1.0f64 / (w as f64 + 1.0f64);
                sx = ((x0 as f64 + 2.0f64 * w as f64 * x1 as f64 + x2 as f64) * xy / 2.0f64 + 0.5f64)
                    .floor() as i32;
                // subdivide curve in half
                sy = ((y0 as f64 + 2.0f64 * w as f64 * y1 as f64 + y2 as f64) * xy / 2.0f64 + 0.5f64)
                    .floor() as i32;
                dx = ((w * x1 as f32 + x0 as f32) as f64 * xy + 0.5f64).floor();
                dy = ((y1 as f32 * w + y0 as f32) as f64 * xy + 0.5f64).floor();
                let next_points: Box<dyn Iterator<Item = PixelAA>> = Box::new(plot_quad_rational_bezier_seg_aa(x0, y0, dx as i32, dy as i32, sx, sy, cur as f32));
                for point in next_points {
                    yield point;
                }
                // plot apart
                dx = ((w * x1 as f32 + x2 as f32) as f64 * xy + 0.5f64).floor();
                dy = ((y1 as f32 * w + y2 as f32) as f64 * xy + 0.5f64).floor();
                let next_points: Box<dyn Iterator<Item = PixelAA>> = Box::new(plot_quad_rational_bezier_seg_aa(
                    sx, sy, dx as i32, dy as i32, x2, y2, cur as f32,
                ));
                for point in next_points {
                    yield point;
                }
                return;
            }
            // error 1st step
            err = dx + dy - xy;
            // pixel loop
            loop {
                cur = (dx - xy).min(xy - dy);
                ed = (dx - xy).max(xy - dy);
                ed += 2_f64 * ed * cur * cur / (4.0f64 * ed * ed + cur * cur);
                // y step
                x1 = (255_f64 * (err - dx - dy + xy).abs() / ed) as i32;
                if x1 < 256 {
                    yield (x0, y0, x1 as u8);
                }
                f = 2_f64 * err + dy < 0 as f64;
                if f {
                    // approximate error distance
                    // get blend value by pixel error
                    // plot curve
                    // y step
                    if y0 == y2 {
                        return;
                    }
                    // last pixel -> curve finished
                    if dx - err < ed {
                        yield (x0 + sx, y0, (255_f64 * (dx - err).abs() / ed) as u8);
                    }
                }
                if 2_f64 * err + dx > 0 as f64 {
                    // x step
                    if x0 == x2 {
                        return;
                    }
                    // last pixel -> curve finished
                    if err - dy < ed {
                        yield (x0, y0 + sy, (255_f64 * (err - dy).abs() / ed) as u8);
                    }
                    x0 += sx;
                    dx += xy;
                    dy += yy;
                    err += dy
                }
                if f {
                    y0 += sy;
                    dy += xy;
                    dx += xx;
                    err += dx
                }
                if dy >= dx {
                    break;
                }
            }
        }
        // plot remaining needle to end
        for point in plot_line_aa(x0, y0, x2, y2) {
            yield point;
        }
    })
}

pub fn plot_cubic_bezier_seg_aa(
    mut x0: i32,
    mut y0: i32,
    mut x1: f32,
    mut y1: f32,
    mut x2: f32,
    mut y2: f32,
    mut x3: i32,
    mut y3: i32,
) -> impl Iterator<Item = PixelAA> {
    gen_iter!(move {
        // plot limited anti-aliased cubic Bezier segment
        let mut f;
        let mut fx;
        let mut fy;
        let mut leg = 1;
        let mut sx = if x0 < x3 { 1 } else { -(1) };
        let mut sy = if y0 < y3 { 1 } else { -(1) };
        // step direction
        let xc: f32 = -((x0 as f32 + x1 - x2 - x3 as f32) as f64).abs() as f32;
        let xa: f32 = xc - (4 * sx) as f32 * (x1 - x2);
        let mut xb: f32 = sx as f32 * (x0 as f32 - x1 - x2 + x3 as f32);
        let yc: f32 = -((y0 as f32 + y1 - y2 - y3 as f32) as f64).abs() as f32;
        let ya: f32 = yc - (4 * sy) as f32 * (y1 - y2);
        let mut yb: f32 = sy as f32 * (y0 as f32 - y1 - y2 + y3 as f32);
        let mut ab: f64;
        let mut ac: f64;
        let mut bc: f64;
        let mut ba: f64;
        let mut xx: f64;
        let mut xy: f64;
        let mut yy: f64;
        let mut dx: f64;
        let mut dy: f64;
        let mut ex: f64;
        let mut px: f64;
        let mut py: f64;
        let mut ed: f64;
        let mut ip: f64;
        let ep: f64 = 0.01f64;
        // check for curve restrains
        // slope P0-P1 == P2-P3     and  (P0-P3 == P1-P2      or  no slope change)
        assert!(
            (((x1 - x0 as f32) * (x2 - x3 as f32)) as f64) < ep
                && ((((x3 - x0) as f32 * (x1 - x2)) as f64) < ep
                    || ((xb * xb) as f64) < (xa * xc) as f64 + ep)
        );
        assert!(
            (((y1 - y0 as f32) * (y2 - y3 as f32)) as f64) < ep
                && ((((y3 - y0) as f32 * (y1 - y2)) as f64) < ep
                    || ((yb * yb) as f64) < (ya * yc) as f64 + ep)
        );
        if xa == 0 as f32 && ya == 0 as f32 {
            // quadratic Bezier
            sx = (((3_f32 * x1 - x0 as f32 + 1_f32) / 2_f32) as f64).floor() as i32;
            sy = (((3_f32 * y1 - y0 as f32 + 1_f32) / 2_f32) as f64).floor() as i32;
            // new midpoint
            for point in plot_quad_bezier_seg_aa(x0, y0, sx, sy, x3, y3) {
                yield point;
            }
            return;
        }
        x1 = (x1 - x0 as f32) * (x1 - x0 as f32) + (y1 - y0 as f32) * (y1 - y0 as f32) + 1_f32;
        // line lengths
        x2 = (x2 - x3 as f32) * (x2 - x3 as f32) + (y2 - y3 as f32) * (y2 - y3 as f32) + 1_f32;
        's_56: loop {
            // loop over both ends
            ab = (xa * yb - xb * ya) as f64;
            ac = (xa * yc - xc * ya) as f64;
            bc = (xb * yc - xc * yb) as f64;
            ip = 4_f64 * ab * bc - ac * ac;
            // self intersection loop at all?
            ex = ab * (ab + ac - 3_f64 * bc) + ac * ac;
            // P0 part of self-intersection loop?
            f = if ex > 0 as f64 {
                1_i32
            } else {
                ((1_f32 + 1024_f32 / x1) as f64).sqrt() as i32
            };
            // calculate resolution
            ab *= f as f64;
            ac *= f as f64;
            bc *= f as f64;
            ex *= (f * f) as f64;
            // increase resolution
            xy = 9_f64 * (ab + ac + bc) / 8_f64;
            ba = (8_f32 * (xa - ya)) as f64;
            // init differences of 1st degree
            dx = 27_f64
                * (8_f64 * ab * (yb * yb - ya * yc) as f64 + ex * (ya + 2_f32 * yb + yc) as f64)
                / 64_f64
                - (ya * ya) as f64 * (xy - ya as f64);
            dy = 27_f64
                * (8_f64 * ab * (xb * xb - xa * xc) as f64 - ex * (xa + 2_f32 * xb + xc) as f64)
                / 64_f64
                - (xa * xa) as f64 * (xy + xa as f64);
            // init differences of 2nd degree
            xx = 3_f64
                * (3_f64 * ab * (3_f32 * yb * yb - ya * ya - 2_f32 * ya * yc) as f64
                    - ya as f64 * (3_f64 * ac * (ya + yb) as f64 + ya as f64 * ba))
                / 4_f64;
            yy = 3_f64
                * (3_f64 * ab * (3_f32 * xb * xb - xa * xa - 2_f32 * xa * xc) as f64
                    - xa as f64 * (3_f64 * ac * (xa + xb) as f64 + xa as f64 * ba))
                / 4_f64;
            xy = (xa * ya) as f64 * (6_f64 * ab + 6_f64 * ac - 3_f64 * bc + ba);
            ac = (ya * ya) as f64;
            ba = (xa * xa) as f64;
            xy = 3_f64
                * (xy + (9 * f) as f64 * (ba * yb as f64 * yc as f64 - (xb * xc) as f64 * ac)
                    - (18_f32 * xb * yb) as f64 * ab)
                / 8_f64;
            if ex < 0 as f64 {
                // negate values if inside self-intersection loop
                dx = -dx;
                dy = -dy;
                xx = -xx;
                yy = -yy;
                xy = -xy;
                ac = -ac;
                ba = -ba
            }
            // init differences of 3rd degree
            ab = (6_f32 * ya) as f64 * ac;
            ac *= (-(6) as f32 * xa) as f64;
            bc = (6_f32 * ya) as f64 * ba;
            ba *= (-(6) as f32 * xa) as f64;
            dx += xy;
            ex = dx + dy;
            dy += xy;
            // error of 1st step
            fy = f;
            fx = fy;
            's_205: loop {
                if !(x0 != x3 && y0 != y3) {
                    break 's_56;
                }
                y1 = (xy - dx).abs().min((dy - xy).abs()) as f32;
                ed = (xy - dx).abs().max((dy - xy).abs());
                // approximate error distance
                ed = f as f64
                    * (ed
                        + 2_f64 * ed * y1 as f64 * y1 as f64
                            / (4_f64 * ed * ed + (y1 * y1) as f64));
                y1 = (255_f64
                    * (ex - (f - fx + 1) as f64 * dx - (f - fy + 1) as f64 * dy + f as f64 * xy).abs()
                    / ed) as f32;
                if y1 < 256_f32 {
                    yield (x0, y0, y1 as u8);
                }
                // plot curve
                px = (ex - (f - fx + 1) as f64 * dx + (fy - 1) as f64 * dy).abs();
                // pixel intensity x move
                py = (ex + (fx - 1) as f64 * dx - (f - fy + 1) as f64 * dy).abs();
                // pixel intensity y move
                y2 = y0 as f32;
                loop
                // move sub-steps of one pixel
                {
                    if ip >= -ep {
                        // intersection possible? -> check..
                        if dx + xx > xy || dy + yy < xy {
                            break 's_205;
                        }
                    }
                    // two x or y steps
                    y1 = (2_f64 * ex + dx) as f32;
                    // save value for test of y step
                    if 2_f64 * ex + dy > 0 as f64 {
                        // x sub-step
                        fx -= 1;
                        dx += xx;
                        ex += dx;
                        xy += ac;
                        dy += xy;
                        yy += bc;
                        xx += ab
                    } else if y1 > 0 as f32 {
                        break 's_205;
                    }
                    // tiny nearly cusp
                    if y1 <= 0 as f32 {
                        // y sub-step
                        fy -= 1;
                        dy += yy;
                        ex += dy;
                        xy += bc;
                        dx += xy;
                        xx += ac;
                        yy += ba
                    }
                    if !(fx > 0 && fy > 0) {
                        break;
                    }
                }
                // pixel complete?
                if 2 * fy <= f {
                    // x+ anti-aliasing pixel
                    if py < ed {
                        yield (x0 + sx, y0, (255_f64 * py / ed) as u8);
                    }
                    // y step
                    y0 += sy;
                    fy += f
                }
                if 2 * fx <= f {
                    // plot curve
                    // y+ anti-aliasing pixel
                    if px < ed {
                        yield (x0, (y2 + sy as f32) as i32, (255_f64 * px / ed) as u8);
                    }
                    // x step
                    x0 += sx;
                    fx += f
                }
            }
            // plot curve
            // finish curve by line
            if 2_f64 * ex < dy && 2 * fy <= f + 2 {
                // round x+ approximation pixel
                if py < ed {
                    yield (x0 + sx, y0, (255_f64 * py / ed) as u8);
                }
                // plot curve
                y0 += sy
            }
            if 2_f64 * ex > dx && 2 * fx <= f + 2 {
                // round y+ approximation pixel
                if px < ed {
                    yield (x0, (y2 + sy as f32) as i32, (255_f64 * px / ed) as u8);
                }
                // plot curve
                x0 += sx
            }
            xx = x0 as f64;
            x0 = x3;
            x3 = xx as i32;
            sx = -sx;
            xb = -xb;
            // swap legs
            yy = y0 as f64;
            y0 = y3;
            y3 = yy as i32;
            sy = -sy;
            yb = -yb;
            x1 = x2;
            let fresh13 = leg;
            leg -= 1;
            if fresh13 == 0 {
                break;
            }
        }
        // try other end
        for point in plot_line_aa(x0, y0, x3, y3) {
            yield point;
        }
        // remaining part in case of cusp or crunode
    })
}

pub fn plot_line_width(
    mut x0: i32,
    mut y0: i32,
    x1: i32,
    y1: i32,
    mut wd: f32,
) -> impl Iterator<Item = PixelAA> {
    gen_iter!(move {
        // plot an anti-aliased line of width wd
        let dx = (x1 - x0).abs();
        let sx = if x0 < x1 { 1 } else { -(1) };
        let dy = (y1 - y0).abs();
        let sy = if y0 < y1 { 1 } else { -(1) };
        let mut err = dx - dy;
        let mut e2;
        let mut x2;
        let mut y2;
        // error value e_xy
        let ed: f32 = if dx + dy == 0 {
            1_f64
        } else {
            ((dx as f32 * dx as f32 + dy as f32 * dy as f32) as f64).sqrt()
        } as f32;
        wd = (wd + 1_f32) / 2_f32;
        loop {
            // pixel loop
            yield (
                x0,
                y0,
                0.0_f64.max((255_f32 * ((err - dx + dy).abs() as f32 / ed - wd + 1_f32)) as f64) as u8,
            );
            e2 = err;
            x2 = x0;
            if 2 * e2 >= -dx {
                // x step
                e2 += dy;
                y2 = y0;
                while (e2 as f32) < ed * wd && (y1 != y2 || dx > dy) {
                    y2 += sy;
                    yield (
                        x0,
                        y2,
                        0.0_f64.max((255_f32 * ((e2).abs() as f32 / ed - wd + 1_f32)) as f64) as u8,
                    );
                    e2 += dx
                }
                if x0 == x1 {
                    break;
                }
                e2 = err;
                err -= dy;
                x0 += sx
            }
            if 2 * e2 > dy {
                continue;
            }
            // y step
            e2 = dx - e2;
            while (e2 as f32) < ed * wd && (x1 != x2 || dx < dy) {
                x2 += sx;
                yield (
                    x2,
                    y0,
                    0.0_f64.max((255_f32 * ((e2).abs() as f32 / ed - wd + 1_f32)) as f64) as u8,
                );
                e2 += dy
            }
            if y0 == y1 {
                break;
            }
            err += dx;
            y0 += sy
        }
    })
}

pub fn plot_quad_spline<'a>(
    x: &'a mut [i32],
    y: &'a mut [i32],
) -> impl Iterator<Item = Pixel> + 'a {
    gen_iter!(move {
        let n = x.len();
        // plot quadratic spline, destroys input arrays x,y
        let mut mi: f32 = 1_f32;
        let mut m: [f32; 6] = [0.; 6];
        // diagonal constants of matrix
        let mut i: usize = 0;
        let mut x0;
        let mut y0;
        let mut x1;
        let mut y1;
        let mut x2 = x[n as usize];
        let mut y2 = y[n as usize];
        // need at least 3 points P[0]..P[n]
        x0 =
            8 * x[1] -
                2 * x[0_usize];
        x[1] = x0;
        // first row of matrix
        y0 =
            8 * y[1] -
                2 * y[0_usize];
        y[1] = y0;
        while i < n {
            // forward sweep
            if (i - 2) < 6 {
                mi = (1.0f64 / (6.0f64 - mi as f64)) as f32;
                m[(i - 2) as usize] = mi
            }
            x0 =
                ((8.0 * x[i] as f32
                        - x0 as f32 * mi) as f64 + 0.5f64).floor() as i32;
            x[i] = x0;
            // store yi
            y0 =
                ((8.0 * y[i] as f32
                        - y0 as f32 * mi) as f64 + 0.5f64).floor() as i32;
            y[i] = y0;
            i += 1
        }
        x1 =
            ((x0 - 2 * x2) as f64 /
                    (5.0f64 - mi as f64) + 0.5f64).floor() as i32;
        // correction last row
        y1 =
            ((y0 - 2 * y2) as f64 /
                    (5.0f64 - mi as f64) + 0.5f64).floor() as i32;
        i = n - 2;
        while i > 0 {
            // back substitution
            if i <= 6 { mi = m[(i - 1) as usize] }
            x0 =
                ((x[i] - x1) as f64 * mi as
                        f64 + 0.5f64).floor() as i32;
            // next corner
            y0 =
                ((y[i] - y1) as f64 * mi as
                        f64 + 0.5f64).floor() as i32;
            for pixel in plot_quad_bezier((x0 + x1) / 2, (y0 + y1) / 2, x1, y1, x2, y2) {
                yield pixel;
            }
            x2 = (x0 + x1) / 2;
            x1 = x0;
            y2 = (y0 + y1) / 2;
            y1 = y0;
            i -= 1
        }
        for pixel in plot_quad_bezier(x[0_usize], y[0_usize], x1, y1, x2, y2) {
            yield pixel;
        }
    })
}

pub fn plot_cubic_spline<'a>(
    x: &'a mut [i32],
    y: &'a mut [i32],
) -> impl Iterator<Item = Pixel> + 'a {
    gen_iter!(move {
        let n = x.len();
        // plot cubic spline, destroys input arrays x,y
        let mut mi: f32 = 0.25f64 as f32;
        let mut m: [f32; 6] = [0.; 6];
        // diagonal constants of matrix
        let mut x3 = x[(n - 1) as usize];
        let mut y3 = y[(n - 1) as usize];
        let mut x4 = x[n as usize];
        let mut y4 = y[n as usize];
        let mut i;
        let mut x0;
        let mut y0;
        let mut x1;
        let mut y1;
        let mut x2;
        let mut y2;

        assert!(n > 2);

        // need at least 4 points P[0]..P[n]
        x0 =
            12 * x[1] -
                3 * x[0_usize];
        x[1] = x0;
        // first row of matrix
        y0 =
            12 * y[1] -
                3 * y[0_usize];
        y[1] = y0;
        i = 2;
        while i < n {
            // foreward sweep
            if (i - 2) < 6 {
                mi = (0.25f64 / (2.0f64 - mi as f64)) as f32;
                m[(i - 2) as usize] = mi
            }
            x0 =
                ((12.0 * x[i] as
                        f32 -
                        (2 * x0) as f32 * mi) as
                        f64 + 0.5f64).floor() as i32;
            x[i] = x0;
            y0 =
                ((12.0 * y[i] as
                        f32 -
                        (2 * y0) as f32 * mi) as
                        f64 + 0.5f64).floor() as i32;
            y[i] = y0;
            i += 1
        }
        x2 =
            (((x0 - 3 * x4) as f32 /
                    (7_f32 -
                            4_f32 * mi)) as
                    f64 + 0.5f64).floor() as i32;
        // correct last row
        y2 =
            (((y0 - 3 * y4) as f32 /
                    (7_f32 -
                            4_f32 * mi)) as
                    f64 + 0.5f64).floor() as i32;
        for pixel in plot_cubic_bezier(x3, y3, (x2 + x4) / 2, (y2 + y4) / 2, x4, y4, x4, y4) {
            yield pixel;
        }
        if (n - 3) < 6 {
            mi = m[(n - 3) as usize]
        }
        x1 =
            (((x[(n - 2) as usize] -
                        2 * x2) as f32 * mi) as
                    f64 + 0.5f64).floor() as i32;
        y1 =
            (((y[(n - 2) as usize] -
                        2 * y2) as f32 * mi) as
                    f64 + 0.5f64).floor() as i32;
        i = n - 3;
        while i > 0 {
            // back substitution
            if i <= 6 { mi = m[(i - 1) as usize] }
            x0 =
                (((x[i] - 2 * x1) as
                        f32 * mi) as f64 + 0.5f64).floor() as
                    i32;
            y0 =
                (((y[i] - 2 * y1) as
                        f32 * mi) as f64 + 0.5f64).floor() as
                    i32;
            x4 =
                ((x0 + 4 * x1 + x2 + 3) as
                        f64 / 6.0f64).floor() as i32;
            // reconstruct P[i]
            y4 =
                ((y0 + 4 * y1 + y2 + 3) as
                        f64 / 6.0f64).floor() as i32;
            for pixel in plot_cubic_bezier(x4, y4,
                            (((2 * x1 + x2) /
                                    3) as f64 +
                                    0.5f64).floor() as i32,
                            (((2 * y1 + y2) /
                                    3) as f64 +
                                    0.5f64).floor() as i32,
                            (((x1 + 2 * x2) /
                                    3) as f64 +
                                    0.5f64).floor() as i32,
                            (((y1 + 2 * y2) /
                                    3) as f64 +
                                    0.5f64).floor() as i32, x3, y3) {
                yield pixel;
            }
            x3 = x4;
            y3 = y4;
            x2 = x1;
            y2 = y1;
            x1 = x0;
            y1 = y0;
            i -= 1
        }
        x0 = x[0_usize];
        x4 =
            ((3 * x0 + 7 * x1 +
                    2 * x2 + 6) as f64
                    / 12.0f64).floor() as i32;
        // reconstruct P[1]
        y0 = y[0_usize];
        y4 =
            ((3 * y0 + 7 * y1 +
                    2 * y2 + 6) as f64
                    / 12.0f64).floor() as i32;
        for pixel in plot_cubic_bezier(x4, y4,
                        (((2 * x1 + x2) / 3) as
                                f64 + 0.5f64).floor() as i32,
                        (((2 * y1 + y2) / 3) as
                                f64 + 0.5f64).floor() as i32,
                        (((x1 + 2 * x2) / 3) as
                                f64 + 0.5f64).floor() as i32,
                        (((y1 + 2 * y2) / 3) as
                                f64 + 0.5f64).floor() as i32, x3,
                        y3) {
            yield pixel;
        }
        for pixel in plot_cubic_bezier(x0, y0, x0, y0, (x0 + x1) / 2, (y0 + y1) / 2, x4, y4) {
            yield pixel;
        }
    })
}
