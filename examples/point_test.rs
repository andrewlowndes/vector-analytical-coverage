use glm::vec2;
use minifb::{Key, Window, WindowOptions};
use std::time::Duration;
use svg::{
    data::cubic_font::cubic_font_shape,
    shapes::{point_in_shape, shape_aabb},
    slice2d::{rgb, Slice2d},
};

const WIDTH: usize = 800;
const HEIGHT: usize = 600;

fn main() {
    let mut buffer = Slice2d {
        data: vec![0; WIDTH * HEIGHT],
        width: WIDTH,
        height: HEIGHT,
    };

    let mut window = Window::new(
        "Point test example - ESC to exit",
        WIDTH,
        HEIGHT,
        WindowOptions::default(),
    )
    .expect("No window");
    window.limit_update_rate(Some(Duration::from_micros(16600)));

    let mut shape = cubic_font_shape();

    //scale and move the shape so we can see it
    let scale = 400.0;
    let translate = vec2(0.0, 400.0);

    for line in &mut shape.lines {
        line.a = (line.a * scale) + translate;
        line.b = (line.b * scale) + translate;
    }

    for quadratic in &mut shape.quadratics {
        quadratic.a = (quadratic.a * scale) + translate;
        quadratic.b = (quadratic.b * scale) + translate;
        quadratic.c = (quadratic.c * scale) + translate;
    }

    for cubic in &mut shape.cubics {
        cubic.a = (cubic.a * scale) + translate;
        cubic.b = (cubic.b * scale) + translate;
        cubic.c = (cubic.c * scale) + translate;
        cubic.d = (cubic.d * scale) + translate;
    }

    let mut offline_buffer = buffer.data.clone();

    let black = rgb(0, 0, 0);
    let red = rgb(255, 0, 0);

    //rasterize each pixel and store in an offline buffer
    let aabb = shape_aabb(&shape);

    let mut i = 0;
    for y in 0..HEIGHT {
        for x in 0..WIDTH {
            if point_in_shape(&shape, &aabb, &vec2(x as f32, y as f32)) {
                offline_buffer[i] = red;
            } else {
                offline_buffer[i] = black;
            }

            i += 1;
        }
    }

    while window.is_open() && !window.is_key_down(Key::Escape) {
        buffer.data.fill(0);

        //at run-time just copy our slice over to the framebuffer
        buffer.data[..].copy_from_slice(&offline_buffer);

        window
            .update_with_buffer(&buffer.data, WIDTH, HEIGHT)
            .unwrap();
    }
}
