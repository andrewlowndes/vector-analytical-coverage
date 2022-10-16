use glm::vec2;
use minifb::{Key, MouseButton, MouseMode, Window, WindowOptions};
use serde::{Deserialize, Serialize};
use std::{
    fs::{read_to_string, File},
    io::Write,
    time::Duration,
};
use svg::{
    data::cubic_font::cubic_font_shape,
    shapes::{clip_shape, Line, Shape},
    slice2d::{rgb, Slice2d},
};

const WIDTH: usize = 800;
const HEIGHT: usize = 600;

#[derive(Clone, Default, PartialEq, Serialize, Deserialize, Debug)]
struct Settings {
    pub delta_x: f32,
    pub delta_y: f32,
}

const SETTINGS_PATH: &str = "settings.json";

fn main() -> anyhow::Result<()> {
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

    let mut shape = cubic_font_shape();

    let red = rgb(255, 0, 0);
    let green = rgb(0, 255, 0);
    let blue = rgb(0, 0, 255);

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

    let mut settings: Settings = read_to_string(SETTINGS_PATH)
        .map_err(anyhow::Error::new)
        .and_then(|str| serde_json::from_str(&str).map_err(anyhow::Error::new))
        .unwrap_or_default();

    let clipping_vertices = vec![
        vec2(200.0 + settings.delta_x, 100.0 + settings.delta_y),
        vec2(200.0 + settings.delta_x, 250.0 + settings.delta_y),
        vec2(350.0 + settings.delta_x, 250.0 + settings.delta_y),
        vec2(350.0 + settings.delta_x, 100.0 + settings.delta_y),
    ];

    let clipping_lines = vec![
        Line(clipping_vertices[0], clipping_vertices[1]),
        Line(clipping_vertices[1], clipping_vertices[2]),
        Line(clipping_vertices[2], clipping_vertices[3]),
        Line(clipping_vertices[3], clipping_vertices[0]),
    ];

    let mut clipping_shape = Shape(clipping_lines, vec![], vec![]);

    let mut mouse_dragging = false;
    let mut last_mouse_pos = (0.0, 0.0);

    while window.is_open() && !window.is_key_down(Key::Escape) {
        //logic
        if window.get_mouse_down(MouseButton::Left) {
            if let Some(mouse_pos) = window.get_mouse_pos(MouseMode::Discard) {
                if mouse_dragging {
                    let delta_x = mouse_pos.0 - last_mouse_pos.0;
                    let delta_y = mouse_pos.1 - last_mouse_pos.1;

                    for line in clipping_shape.lines.iter_mut() {
                        line.a.x += delta_x;
                        line.a.y += delta_y;
                        line.b.x += delta_x;
                        line.b.y += delta_y;
                    }
                } else {
                    mouse_dragging = true;
                }

                last_mouse_pos = mouse_pos;
            }
        } else {
            mouse_dragging = false;
        }

        //draw
        buffer.data.fill(0);

        //draw the outline of the shapes for now for reference
        buffer.shape(&shape, red);

        //draw the clipping shape
        buffer.shape(&clipping_shape, green);

        //intersect a clipping shape and produce a new shape and then overlay
        let clipped_shape = clip_shape(&shape, &clipping_shape.lines);
        buffer.shape(&clipped_shape, blue);

        window.update_with_buffer(&buffer.data, WIDTH, HEIGHT)?
    }

    //store the coords in a json file so we can re-use them next time it opens
    settings.delta_x = clipping_shape.lines[0].a.x - 200.0;
    settings.delta_y = clipping_shape.lines[0].a.y - 100.0;

    File::create(SETTINGS_PATH)?.write_all(&serde_json::to_vec(&settings)?)?;

    Ok(())
}
