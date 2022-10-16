use glm::vec2;
use minifb::{Key, Window, WindowOptions};
use std::time::Duration;
use svg::{
    shapes::{line_quadratic_intersect, Line, Quadratic},
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
        "Cubic example - ESC to exit",
        WIDTH,
        HEIGHT,
        WindowOptions::default(),
    )
    .expect("No window");
    window.limit_update_rate(Some(Duration::from_micros(16600)));

    //let line1 = Line(vec2(50.0, 120.0), vec2(300.0, 360.0));
    let line2 = Line(vec2(100.0, 340.0), vec2(400.0, 20.0));
    let quadratic = Quadratic(vec2(150.0, 120.0), vec2(100.0, 300.0), vec2(300.0, 360.0));

    let red = rgb(255, 0, 0);
    let green = rgb(0, 255, 0);
    let blue = rgb(0, 0, 255);

    let point_radius = 2;
    let point_diameter = point_radius * 2;

    while window.is_open() && !window.is_key_down(Key::Escape) {
        buffer.data.fill(0);

        //buffer.line(line1.a.x as i32, line1.a.y as i32, line1.b.x as i32, line1.b.y as i32, red);
        buffer.quadratic(
            quadratic.a.x as i32,
            quadratic.a.y as i32,
            quadratic.b.x as i32,
            quadratic.b.y as i32,
            quadratic.c.x as i32,
            quadratic.c.y as i32,
            red,
        );
        buffer.line(
            line2.a.x as i32,
            line2.a.y as i32,
            line2.b.x as i32,
            line2.b.y as i32,
            green,
        );

        for intersection in line_quadratic_intersect(&line2, &quadratic) {
            buffer.rectangle(
                intersection.pos.x as i32 - point_radius,
                intersection.pos.y as i32 - point_radius,
                point_diameter as usize,
                point_diameter as usize,
                blue,
            );
        }

        window
            .update_with_buffer(&buffer.data, WIDTH, HEIGHT)
            .unwrap();
    }
}
