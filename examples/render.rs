use glm::{vec2, Vec2, GenNumVec};
use minifb::{Key, Window, WindowOptions};
use std::{time::Duration, f32::EPSILON};
use svg::{
    data::cubic_font::cubic_font_shape,
    shapes::{clip_shape, in_range, shape_area, Line, Shape, shape_aabb, transform_shape, aabb_aabb_intersect, cubic_aabb, line_cubic_intersect, line_cubic_intersect_debug, point_in_shape},
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
    transform_shape(&mut shape, scale, translate);
    let aabb = shape_aabb(&shape);

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

            let mut clipping_shape = Shape(clipping_lines, vec![], vec![]);

            let (mut clipped_shape, intersections) = clip_shape(&shape, &clipping_shape.lines);

            let clipped_area = -shape_area(&clipped_shape);
            let alpha = (clipped_area * 255.0).floor() as u8;

            if clipped_area > 1.1 {
                //if the area exceeds the 0-1 range (allow for a margin) then there is an error so show it
                /*
                dbg!(&intersections);

                let clipping_aabb = shape_aabb(&clipping_shape);
                for joint in &shape.cubics {
                    if !aabb_aabb_intersect(&clipping_aabb, &cubic_aabb(joint)) {
                        continue;
                    }

                    for (_clip_index, clip_joint) in clipping_shape.lines.iter().enumerate() {
                        dbg!(joint, line_cubic_intersect_debug(clip_joint, joint));
                    }
                }

                for (_clip_index, clip_joint) in clipping_shape.lines.iter().enumerate() {
                    dbg!(point_in_shape(&shape, &aabb, &clip_joint.a, true));
                }
    
                let clipped_aabb = shape_aabb(&clipped_shape);
                let clipped_square = clipped_aabb.max - clipped_aabb.min;
                let clipped_square_area = clipped_square.x * clipped_square.y;

                //zoom into the clipping shape so we can see what is going on
                let clipping_aabb = shape_aabb(&clipping_shape);
                let scale = 200.0;
                let translate = (-clipping_aabb.min * scale) + 100.0;

                transform_shape(&mut shape, scale, translate);
                transform_shape(&mut clipping_shape, scale, translate);
                transform_shape(&mut clipped_shape, scale, translate);
                
                dbg!(&clipped_shape, &clipped_area, &alpha, &clipped_square_area);
                offline_buffer.data.fill(0);

                offline_buffer.shape(&shape, red);
                offline_buffer.shape(&clipping_shape, green);
                offline_buffer.shape(&clipped_shape, blue);
                break 'outer;
                */
                
                //draw the pixel trying to shade as green and then carry on (good for showing number of incorrect pixels)
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
