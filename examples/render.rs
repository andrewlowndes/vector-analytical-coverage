use glm::vec2;
use minifb::{Key, Window, WindowOptions};
use std::time::Duration;
use svg::{
    shapes::{Line, in_range, clip_shape, shape_area, Shape},
    slice2d::{rgb, Slice2d}, data::cubic_font::cubic_font_shape,
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
        "Render example - ESC to exit",
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

    let mut offline_buffer = buffer.clone();

    let red = rgb(255, 0, 0);
    let green = rgb(0, 255, 0);
    let blue = rgb(0, 0, 255);

    //rasterize each pixel and store in an offline buffer
    let mut i = 0;
    'outer: for y in 0..HEIGHT {
        for x in 0..WIDTH {
            let clipping_vertices = &[
                vec2(x as f32, y as f32),
                vec2(x as f32, y as f32 + 1.0),
                vec2(x as f32 + 1.0, y as f32 + 1.0),
                vec2(x as f32 + 1.0, y as f32),
            ];

            let clipping_lines = vec![
                Line(clipping_vertices[0], clipping_vertices[1]),
                Line(clipping_vertices[1], clipping_vertices[2]),
                Line(clipping_vertices[2], clipping_vertices[3]),
                Line(clipping_vertices[3], clipping_vertices[0]),
            ];

            let clipping_shape = Shape(clipping_lines, vec![], vec![]);
            
            let clipped_shape = clip_shape(&shape, &clipping_shape.lines);

            let clipped_area = -shape_area(&clipped_shape);
            let alpha = (clipped_area * 255.0).floor() as u8;

            if !in_range(clipped_area, 0.0, 1.0) {
                /*
                dbg!(&clipped_shape, &clipped_area, &alpha);
                offline_buffer.data.fill(0);
                offline_buffer.shape(&shape, red);
                offline_buffer.shape(&clipping_shape, green);
                offline_buffer.shape(&clipped_shape, blue);
                break 'outer;
                */

                offline_buffer.data[i] = green;
                i += 1;
                continue;
            }

            offline_buffer.data[i] = rgb(alpha, 0, 0);
            i += 1;
        }
    }

    while window.is_open() && !window.is_key_down(Key::Escape) {
        buffer.data.fill(0);

        //at run-time just copy our slice over to the framebuffer
        buffer.data[..].copy_from_slice(&offline_buffer.data);

        window
            .update_with_buffer(&buffer.data, WIDTH, HEIGHT)
            .unwrap();
    }
}
