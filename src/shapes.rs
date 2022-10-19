use glm::{abs, acos, atan2, cos, distance, max, min, mix_s, pow, sign, sin, sqrt, vec2, Vec2};
use itertools::Itertools;
use std::{cmp::Ordering, collections::VecDeque};

#[derive(Clone, PartialEq, Debug)]
pub struct Aabb {
    pub min: Vec2,
    pub max: Vec2,
}

#[allow(non_snake_case)]
pub fn Aabb(min: Vec2, max: Vec2) -> Aabb {
    Aabb { min, max }
}

#[derive(Clone, PartialEq, Debug)]
pub struct Line {
    pub a: Vec2,
    pub b: Vec2,
}

#[allow(non_snake_case)]
pub fn Line(a: Vec2, b: Vec2) -> Line {
    Line { a, b }
}

#[derive(Clone, PartialEq, Debug)]
pub struct Quadratic {
    pub a: Vec2,
    pub b: Vec2,
    pub c: Vec2,
}

#[allow(non_snake_case)]
pub fn Quadratic(a: Vec2, b: Vec2, c: Vec2) -> Quadratic {
    Quadratic { a, b, c }
}

#[derive(Clone, PartialEq, Debug)]
pub struct Cubic {
    pub a: Vec2,
    pub b: Vec2,
    pub c: Vec2,
    pub d: Vec2,
}

#[allow(non_snake_case)]
pub fn Cubic(a: Vec2, b: Vec2, c: Vec2, d: Vec2) -> Cubic {
    Cubic { a, b, c, d }
}

#[derive(Clone, PartialEq, Debug)]
pub struct Shape {
    pub lines: Vec<Line>,
    pub quadratics: Vec<Quadratic>,
    pub cubics: Vec<Cubic>,
}

#[allow(non_snake_case)]
pub fn Shape(lines: Vec<Line>, quadratics: Vec<Quadratic>, cubics: Vec<Cubic>) -> Shape {
    Shape {
        lines,
        quadratics,
        cubics,
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Debug)]
pub enum IntersectLocation {
    OnEdge,
    Outside,
    OutsideToInside,
    Inside,
    InsideToOutside,
}

impl IntersectLocation {
    fn from_winding(val: f32) -> Self {
        match val.total_cmp(&0.0) {
            Ordering::Less => IntersectLocation::OutsideToInside,
            Ordering::Greater => IntersectLocation::InsideToOutside,
            Ordering::Equal => IntersectLocation::OnEdge,
        }
    }
}

#[derive(Clone, PartialEq, Debug)]
pub struct Intersection {
    pub t: f32,
    pub line_t: f32,
    pub pos: Vec2,
    pub location: IntersectLocation,
}

#[allow(non_snake_case)]
pub fn Intersection(t: f32, line_t: f32, pos: Vec2, location: IntersectLocation) -> Intersection {
    Intersection {
        t,
        line_t,
        pos,
        location,
    }
}

const EPSILON: f32 = f32::EPSILON; //this will need to be tinkered based on the precision of the float
const TAU: f32 = std::f32::consts::TAU;

pub fn transform_vec(vec: &mut Vec2, scale: f32, translate: Vec2) {
    *vec = *vec * scale + translate;
}

pub fn transform_shape(shape: &mut Shape, scale: f32, translate: Vec2) {
    for line in &mut shape.lines {
        transform_vec(&mut line.a, scale, translate);
        transform_vec(&mut line.b, scale, translate);
    }
    for quadratic in &mut shape.quadratics {
        transform_vec(&mut quadratic.a, scale, translate);
        transform_vec(&mut quadratic.b, scale, translate);
        transform_vec(&mut quadratic.c, scale, translate);
    }
    for cubic in &mut shape.cubics {
        transform_vec(&mut cubic.a, scale, translate);
        transform_vec(&mut cubic.b, scale, translate);
        transform_vec(&mut cubic.c, scale, translate);
        transform_vec(&mut cubic.d, scale, translate);
    }
}

pub fn in_range(val: f32, min_val: f32, max_val: f32) -> bool {
    (min_val..=max_val).contains(&val)
}

pub fn det2(a: Vec2, b: Vec2) -> f32 {
    a.x * b.y - b.x * a.y
}

pub fn line_dir(line: &Line) -> Vec2 {
    line.b - line.a
}

pub fn line_angle(line: &Line) -> f32 {
    let dir = line.b - line.a;
    -atan2(dir.x, dir.y)
}

pub fn line_size(line: &Line) -> f32 {
    distance(line.a, line.b)
}

pub fn line_aabb(line: &Line) -> Aabb {
    Aabb(min(line.a, line.b), max(line.a, line.b))
}

pub fn quadratic_aabb(test: &Quadratic) -> Aabb {
    Aabb(
        min(min(test.a, test.b), test.c),
        max(max(test.a, test.b), test.c),
    )
}

pub fn cubic_aabb(test: &Cubic) -> Aabb {
    Aabb(
        min(min(min(test.a, test.b), test.c), test.d),
        max(max(max(test.a, test.b), test.c), test.d),
    )
}

pub fn lines_aabb(lines: &[Line]) -> Aabb {
    let mut aabb = Aabb(vec2(f32::MAX, f32::MAX), vec2(f32::MIN, f32::MIN));

    for Line { a, b } in lines {
        aabb.min.x = min(a.x, min(b.x, aabb.min.x));
        aabb.min.y = min(a.y, min(b.y, aabb.min.y));
        aabb.max.x = max(a.x, max(b.x, aabb.max.x));
        aabb.max.y = max(a.y, max(b.y, aabb.max.y));
    }

    aabb
}

pub fn shape_aabb(shape: &Shape) -> Aabb {
    let mut aabb = lines_aabb(&shape.lines);

    for Quadratic { a, b, c } in &shape.quadratics {
        aabb.min.x = min(a.x, min(b.x, min(c.x, aabb.min.x)));
        aabb.min.y = min(a.y, min(b.y, min(c.y, aabb.min.y)));
        aabb.max.x = max(a.x, max(b.x, max(c.x, aabb.max.x)));
        aabb.max.y = max(a.y, max(b.y, max(c.y, aabb.max.y)));
    }

    for Cubic { a, b, c, d } in &shape.cubics {
        aabb.min.x = min(a.x, min(b.x, min(c.x, min(d.x, aabb.min.x))));
        aabb.min.y = min(a.y, min(b.y, min(c.y, min(d.y, aabb.min.y))));
        aabb.max.x = max(a.x, max(b.x, max(c.x, max(d.x, aabb.max.x))));
        aabb.max.y = max(a.y, max(b.y, max(c.y, max(d.y, aabb.max.y))));
    }

    aabb
}

//http://ich.deanmcnamee.com/graphics/2016/03/30/CurveArea.html
pub fn line_area(line: &Line) -> f32 {
    det2(line.a, line.b)
}

pub fn quadratic_area(quadratic: &Quadratic) -> f32 {
    (quadratic.c.x * (-quadratic.a.y - 2.0 * quadratic.b.y)
        + 2.0 * quadratic.b.x * (quadratic.c.y - quadratic.a.y)
        + quadratic.a.x * (2.0 * quadratic.b.y + quadratic.c.y))
        / 3.0
}

pub fn cubic_area(cubic: &Cubic) -> f32 {
    (cubic.d.x * (-cubic.a.y - 3.0 * cubic.b.y - 6.0 * cubic.c.y)
        - 3.0 * cubic.c.x * (cubic.a.y + cubic.b.y - 2.0 * cubic.d.y)
        + 3.0 * cubic.b.x * (-2.0 * cubic.a.y + cubic.c.y + cubic.d.y)
        + cubic.a.x * (6.0 * cubic.b.y + 3.0 * cubic.c.y + cubic.d.y))
        / 10.0
}

pub fn shape_area(shape: &Shape) -> f32 {
    let mut area = 0.0;

    for line in &shape.lines {
        area += line_area(line);
    }

    for quadratic in &shape.quadratics {
        area += quadratic_area(quadratic);
    }

    for cubic in &shape.cubics {
        area += cubic_area(cubic);
    }

    area / 2.0
}

pub fn point_in_lines(lines: &[Line], point: &Vec2) -> bool {
    //check the point is on the positive side of the joint
    //for now just use one algorithm
    
    lines
        .iter()
        .all(|joint| det2(*point - joint.a, line_dir(joint)) >= 0.0)
    
    //point_in_shape(&Shape(lines.to_vec(), vec![], vec![]), &lines_aabb(&lines), point)
}

pub fn point_in_aabb(aabb: &Aabb, point: &Vec2) -> bool {
    point.x <= aabb.max.x && point.x >= aabb.min.x && point.y <= aabb.max.y && point.y >= aabb.min.y
}

pub fn aabb_aabb_intersect(a: &Aabb, b: &Aabb) -> bool {
    a.min.x <= b.max.x && a.max.x >= b.min.x && a.min.y <= b.max.y && a.max.y >= b.min.y
}

pub fn cuberoot(x: f32) -> f32 {
    sign(x) * pow(abs(x), 1.0 / 3.0)
}

//utilities
pub fn linear(a: f32, b: f32) -> Vec<f32> {
    if abs(a) < EPSILON {
        return vec![];
    }

    vec![-b / a]
}

pub fn quadratic(a: f32, b: f32, c: f32) -> Vec<f32> {
    if abs(a) < EPSILON {
        return linear(b, c);
    }

    let nom = sqrt(b * b - (4.0 * a * c));
    let denom = 2.0 * a;

    vec![(-b + nom) / denom, (-b - nom) / denom]
}


pub fn cubic(a: f32, b: f32, c: f32, d: f32) -> Vec<f32> {
    if abs(a) < EPSILON {
        return quadratic(b, c, d);
    }

    // Convert to depressed cubic t^3+pt+q = 0 (subst x = t - b/3a)
    let p = (3.0 * a * c - b * b) / (3.0 * a * a);
    let q = (2.0 * b * b * b - 9.0 * a * b * c + 27.0 * a * a * d) / (27.0 * a * a * a);

    let roots: Vec<f32>;

    if abs(p) < EPSILON {
        // p = 0 -> t^3 = -q -> t = -q^1/3
        roots = vec![cuberoot(-q)];
    } else if abs(q) < EPSILON {
        // q = 0 -> t^3 + pt = 0 -> t(t^2+p)=0
        if p < 0.0 {
            roots = vec![0.0, sqrt(-p), -sqrt(-p)];
        } else {
            roots = vec![0.0];
        }
    } else {
        let temp_val = q * q / 4.0 + p * p * p / 27.0;

        if abs(temp_val) < EPSILON {
            // D = 0 -> two roots
            roots = vec![-1.5 * q / p, 3.0 * q / p];
        } else if temp_val > 0.0 {
            // Only one real root
            let u = cuberoot(-q / 2.0 - sqrt(temp_val));
            roots = vec![u - p / (3.0 * u)];
        } else {
            // D < 0, three roots, but needs to use complex numbers/trigonometric solution
            let u = 2.0 * sqrt(-p / 3.0);
            let t = acos(3.0 * q / p / u) / 3.0; // D < 0 implies p < 0 and acos argument in [-1..1]
            let k = TAU / 3.0;
            roots = vec![u * cos(t), u * cos(t - k), u * cos(t - 2.0 * k)];
        }
    }

    // Convert back from depressed cubic
    let divisor = b / (3.0 * a);
    return roots.iter().map(|val| val - divisor).collect::<Vec<_>>();
}

pub fn rotate_point(point: Vec2, angle: f32) -> Vec2 {
    vec2(
        point.x * cos(angle) - point.y * sin(angle),
        point.x * sin(angle) + point.y * cos(angle),
    )
}

const ROUNDING_T: f32 = 0.0001;

pub fn clean_t(t: f32) -> f32 {
    //t = parseFloat(t.toFixed(ROUNDING_PRECISION));
    //console.log(t);

    //allow a little overlap at the extremes
    if abs(t) < ROUNDING_T {
        0.0
    } else if abs(t - 1.0) < ROUNDING_T {
        1.0
    } else {
        t
    }
}

pub fn line_line_intersect(line: &Line, test: &Line) -> Vec<Intersection> {
    //align the two lines and then solve the equation
    let angle = line_angle(line);
    let test_a = rotate_point(test.a - line.a, angle);
    let test_b = rotate_point(test.b - line.a, angle);

    let a = test_b.y - test_a.y;
    let b = test_a.y;

    let derivative = sign(a);
    let line_mag = line_size(line);

    let ts = if abs(test_a.y) < EPSILON && abs(test_b.y) < EPSILON {
        //parallel lines, 2 or 0 intersections based on the entry and exit points
        let test_mag = line_size(test);
        let line_a = rotate_point(line.a - test.a, angle);
        let line_b = rotate_point(line.b - test.a, angle);

        vec![0.0, 1.0, line_a.x / test_mag, line_b.x / test_mag]
    } else {
        linear(a, b).into_iter().map(clean_t).collect::<Vec<_>>()
    };

    ts.into_iter()
        .filter(|t| in_range(*t, 0.0, 1.0))
        .filter_map(|t| {
            let line_pos = mix_s(test_a, test_b, t);
            let line_t = clean_t(line_pos.x / line_mag);

            if in_range(line_t, 0.0, 1.0) {
                let pos = mix_s(line.a, line.b, line_t);

                Some(Intersection(
                    t,
                    line_t,
                    pos,
                    IntersectLocation::from_winding(derivative),
                ))
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
}

pub fn line_quadratic_intersect(line: &Line, test: &Quadratic) -> Vec<Intersection> {
    //align the cubic curve to the line and then solve the equation
    let angle = line_angle(line);
    let test_a = rotate_point(test.a - line.a, angle);
    let test_b = rotate_point(test.b - line.a, angle);
    let test_c = rotate_point(test.c - line.a, angle);

    let a = test_a.y - 2.0 * test_b.y + test_c.y;
    let b = 2.0 * (test_b.y - test_a.y);
    let c = test_a.y;

    let derivative_a = 2.0 * a;
    let derivative_b = b;

    let line_mag = line_size(line);

    quadratic(a, b, c)
        .into_iter()
        .map(clean_t)
        .filter(|t| in_range(*t, 0.0, 1.0))
        .filter_map(|t| {
            let line_pos = mix_s(mix_s(test_a, test_b, t), mix_s(test_b, test_c, t), t);

            let line_t = clean_t(line_pos.x / line_mag);

            if in_range(line_t, 0.0, 1.0) {
                let pos = mix_s(line.a, line.b, line_t);

                Some(Intersection(
                    t,
                    line_t,
                    pos,
                    IntersectLocation::from_winding(sign(derivative_a * t + derivative_b)),
                ))
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
}

pub fn line_cubic_intersect(line: &Line, test: &Cubic) -> Vec<Intersection> {
    //align the cubic curve to the line and then solve the equation
    let angle = line_angle(line);
    let test_a = rotate_point(test.a - line.a, angle);
    let test_b = rotate_point(test.b - line.a, angle);
    let test_c = rotate_point(test.c - line.a, angle);
    let test_d = rotate_point(test.d - line.a, angle);

    let a = -test_a.y + 3.0 * (test_b.y - test_c.y) + test_d.y;
    let b = 3.0 * (test_a.y - 2.0 * test_b.y + test_c.y);
    let c = 3.0 * (-test_a.y + test_b.y);
    let d = test_a.y;

    let derivative_a = 3.0 * a;
    let derivative_b = 2.0 * b;
    let derivative_c = c;

    //rather than re-project the points back out
    let ts = cubic(a, b, c, d).into_iter().map(clean_t);

    let line_mag = line_size(line);

    ts.into_iter()
        .filter(|t| in_range(*t, 0.0, 1.0))
        .filter_map(|t| {
            let mid = mix_s(test_b, test_c, t);

            let line_pos = mix_s(
                mix_s(mix_s(test_a, test_b, t), mid, t),
                mix_s(mid, mix_s(test_c, test_d, t), t),
                t,
            );

            let line_t = clean_t(line_pos.x / line_mag);

            if in_range(line_t, 0.0, 1.0) {
                let pos = mix_s(line.a, line.b, line_t);

                Some(Intersection(
                    t,
                    line_t,
                    pos,
                    IntersectLocation::from_winding(sign(
                        derivative_a * t * t + derivative_b * t + derivative_c,
                    )),
                ))
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
}

pub fn line_cubic_intersect_debug(line: &Line, test: &Cubic) -> Vec<Intersection> {
    //align the cubic curve to the line and then solve the equation
    let angle = line_angle(line);
    let test_a = rotate_point(test.a - line.a, angle);
    let test_b = rotate_point(test.b - line.a, angle);
    let test_c = rotate_point(test.c - line.a, angle);
    let test_d = rotate_point(test.d - line.a, angle);

    let a = -test_a.y + 3.0 * (test_b.y - test_c.y) + test_d.y;
    let b = 3.0 * (test_a.y - 2.0 * test_b.y + test_c.y);
    let c = 3.0 * (-test_a.y + test_b.y);
    let d = test_a.y;

    let derivative_a = 3.0 * a;
    let derivative_b = 2.0 * b;
    let derivative_c = c;

    //rather than re-project the points back out
    let ts = cubic(a, b, c, d).into_iter().map(clean_t);

    let line_mag = line_size(line);

    ts.into_iter()
        .filter(|t| in_range(*t, 0.0, 1.0))
        .filter_map(|t| {
            let mid = mix_s(test_b, test_c, t);

            let line_pos = mix_s(
                mix_s(mix_s(test_a, test_b, t), mid, t),
                mix_s(mid, mix_s(test_c, test_d, t), t),
                t,
            );

            let line_t = clean_t(line_pos.x / line_mag);
            dbg!(&line_t);
        

            if in_range(line_t, 0.0, 1.0) {
                let pos = mix_s(line.a, line.b, line_t);

                Some(Intersection(
                    t,
                    line_t,
                    pos,
                    IntersectLocation::from_winding(sign(
                        derivative_a * t * t + derivative_b * t + derivative_c,
                    )),
                ))
            } else {
                None
            }
        })
        .collect::<Vec<_>>()
}

pub fn line_quadratic_split(
    joint: &Quadratic,
    from_t: f32,
    to_t: f32,
    start_point: &Vec2,
    end_point: &Vec2,
) -> Quadratic {
    let mid_point = mix_s(
        mix_s(joint.a, joint.b, from_t),
        mix_s(joint.b, joint.c, from_t),
        to_t,
    );

    Quadratic(*start_point, mid_point, *end_point)
}

pub fn line_cubic_split(
    joint: &Cubic,
    from_t: f32,
    to_t: f32,
    start_point: &Vec2,
    end_point: &Vec2,
) -> Cubic {
    let mid_point1 = mix_s(joint.a, joint.b, from_t);
    let mid_point2 = mix_s(joint.b, joint.c, from_t);
    let mid_point3 = mix_s(joint.c, joint.d, from_t);

    let control_point1 = mix_s(
        mix_s(mid_point1, mid_point2, from_t),
        mix_s(mid_point2, mid_point3, from_t),
        to_t,
    );
    let control_point2 = mix_s(
        mix_s(mid_point1, mid_point2, to_t),
        mix_s(mid_point2, mid_point3, to_t),
        to_t,
    );

    Cubic(*start_point, control_point1, control_point2, *end_point)
}

pub fn point_in_shape(shape: &Shape, aabb: &Aabb, point: &Vec2, debug: bool) -> bool {
    if !point_in_aabb(aabb, point) {
        return false;
    }

    let scanline = Line(*point, vec2(aabb.min.x, point.y));

    //rather than sum the windings, just find the closest line and determine the 'direction' of the line
    let mut max_x = f32::MIN;
    let mut max_intersection: Option<Intersection> = None;

    for joint in &shape.lines {
        let joint_aabb = line_aabb(joint);

        if point.y >= joint_aabb.min.y && point.y <= joint_aabb.max.y {
            for intersection in line_line_intersect(&scanline, joint) {
                if intersection.location != IntersectLocation::OnEdge && in_range(intersection.pos.x, max_x, point.x)
                {
                    max_x = intersection.pos.x;
                    max_intersection = Some(intersection);
                }
            }
        }
    }

    for joint in &shape.quadratics {
        let joint_aabb = quadratic_aabb(joint);

        if point.y >= joint_aabb.min.y && point.y <= joint_aabb.max.y {
            for intersection in line_quadratic_intersect(&scanline, joint) {
                if intersection.location != IntersectLocation::OnEdge && in_range(intersection.pos.x, max_x, point.x)
                {
                    max_x = intersection.pos.x;
                    max_intersection = Some(intersection);
                }
            }
        }
    }

    for joint in &shape.cubics {
        let joint_aabb = cubic_aabb(joint);

        if point.y >= joint_aabb.min.y && point.y <= joint_aabb.max.y {
            for intersection in line_cubic_intersect(&scanline, joint) {
                if intersection.location != IntersectLocation::OnEdge && in_range(intersection.pos.x, max_x, point.x)
                {
                    max_x = intersection.pos.x;
                    max_intersection = Some(intersection);
                }
            }
        }
    }

    if debug {
        dbg!(&max_intersection, point);
    }

    max_intersection
        .map(|intersection| intersection.location == IntersectLocation::OutsideToInside || intersection.pos == *point)
        .unwrap_or_default()
}

pub fn clip_shape(shape: &Shape, clipping_lines: &[Line]) -> (Shape, VecDeque<Vec<Intersection>>) {
    let clipping_lines_aabb = lines_aabb(clipping_lines);
    let aabb = shape_aabb(shape);

    if !aabb_aabb_intersect(&clipping_lines_aabb, &aabb) {
        return (Shape(vec![], vec![], vec![]), VecDeque::from(vec![]));
    }

    //TODO: replace, should not use any complex structures
    let mut clip_joint_intersections_vec: Vec<Vec<Intersection>> = vec![vec![]; clipping_lines.len()];
    clip_joint_intersections_vec.fill_with(Default::default);

    let mut clip_joint_intersections = VecDeque::from(clip_joint_intersections_vec);

    let mut lines = vec![];
    for joint in &shape.lines {
        //test aaabb of the edge against the clipping shape (completely inside or outside)
        if !aabb_aabb_intersect(&clipping_lines_aabb, &line_aabb(joint)) {
            continue;
        }

        //test edge vertices are inside or outside clipping shape
        //kinda pointless as ineffective for curves so just get all of the intersections
        let mut intersections = vec![];
        for (clip_index, clip_joint) in clipping_lines.iter().enumerate() {
            let mut clip_intersections = line_line_intersect(clip_joint, joint); //.filter(intersection => intersection.winding != 0);
            intersections.append(&mut clip_intersections.clone()); //TODO: refactor to avoid clone
            clip_joint_intersections[clip_index].append(&mut clip_intersections);
        }

        //split the curve up at it's intersections in the correct order
        intersections.sort_by(|a, b| a.t.total_cmp(&b.t));

        //if the intersection point is going outside
        let location = if point_in_lines(clipping_lines, &joint.a) {
            IntersectLocation::Inside
        } else {
            IntersectLocation::Outside
        };

        let mut prev_intersection = Intersection(0.0, 0.0, joint.a, location);
        for intersection in intersections {
            if [
                IntersectLocation::OutsideToInside,
                IntersectLocation::Inside,
            ]
            .contains(&prev_intersection.location)
                && [
                    IntersectLocation::InsideToOutside,
                    IntersectLocation::Inside,
                ]
                .contains(&intersection.location)
            {
                lines.push(Line(prev_intersection.pos, intersection.pos));
            }
            prev_intersection = intersection;
        }

        //repeat for the last point in the joint too
        if [
            IntersectLocation::OutsideToInside,
            IntersectLocation::Inside,
        ]
        .contains(&prev_intersection.location)
            && point_in_lines(clipping_lines, &joint.b)
        {
            lines.push(Line(prev_intersection.pos, joint.b));
        }
    }

    //again for quadratics
    let mut quadratics = vec![];
    for joint in &shape.quadratics {
        //test aaabb of the edge against the clipping shape (completely inside or outside)
        if !aabb_aabb_intersect(&clipping_lines_aabb, &quadratic_aabb(joint)) {
            continue;
        }

        //test edge vertices are inside or outside clipping shape
        //kinda pointless as ineffective for curves so just get all of the intersections
        let mut intersections = vec![];
        for (clip_index, clip_joint) in clipping_lines.iter().enumerate() {
            let mut clip_intersections = line_quadratic_intersect(clip_joint, joint); //.filter(intersection => intersection.winding != 0);
            intersections.append(&mut clip_intersections.clone()); //TODO: refactor to avoid clone
            clip_joint_intersections[clip_index].append(&mut clip_intersections);
        }

        //split the curve up at it's intersections in the correct order
        intersections.sort_by(|a, b| a.t.total_cmp(&b.t));

        //if the intersection point is going outside
        let location = if point_in_lines(clipping_lines, &joint.a) {
            IntersectLocation::Inside
        } else {
            IntersectLocation::Outside
        };

        let mut prev_intersection = Intersection(0.0, 0.0, joint.a, location);
        for intersection in intersections {
            if [
                IntersectLocation::OutsideToInside,
                IntersectLocation::Inside,
            ]
            .contains(&prev_intersection.location)
                && [
                    IntersectLocation::InsideToOutside,
                    IntersectLocation::Inside,
                ]
                .contains(&intersection.location)
            {
                quadratics.push(line_quadratic_split(
                    joint,
                    prev_intersection.t,
                    intersection.t,
                    &prev_intersection.pos,
                    &intersection.pos,
                ));
            }
            prev_intersection = intersection;
        }

        //repeat for the last point in the joint too
        if [
            IntersectLocation::OutsideToInside,
            IntersectLocation::Inside,
        ]
        .contains(&prev_intersection.location)
            && point_in_lines(clipping_lines, &joint.c)
        {
            quadratics.push(line_quadratic_split(
                joint,
                prev_intersection.t,
                1.0,
                &prev_intersection.pos,
                &joint.c,
            ));
        }
    }

    //again for cubics
    let mut cubics = vec![];
    for joint in &shape.cubics {
        if !aabb_aabb_intersect(&clipping_lines_aabb, &cubic_aabb(joint)) {
            continue;
        }

        let mut intersections = vec![];
        for (clip_index, clip_joint) in clipping_lines.iter().enumerate() {
            let mut clip_intersections = line_cubic_intersect(clip_joint, joint); //.filter(intersection => intersection.winding != 0);
            intersections.append(&mut clip_intersections.clone()); //TODO: refactor to avoid clone
            clip_joint_intersections[clip_index].append(&mut clip_intersections);
        }

        intersections.sort_by(|a, b| a.t.total_cmp(&b.t));

        let location = if point_in_lines(clipping_lines, &joint.a) {
            IntersectLocation::Inside
        } else {
            IntersectLocation::Outside
        };

        let mut prev_intersection = Intersection(0.0, 0.0, joint.a, location);
        for intersection in intersections {
            if [
                IntersectLocation::OutsideToInside,
                IntersectLocation::Inside,
            ]
            .contains(&prev_intersection.location)
                && [
                    IntersectLocation::InsideToOutside,
                    IntersectLocation::Inside,
                ]
                .contains(&intersection.location)
            {
                cubics.push(line_cubic_split(
                    joint,
                    prev_intersection.t,
                    intersection.t,
                    &prev_intersection.pos,
                    &intersection.pos,
                ));
            }
            prev_intersection = intersection;
        }

        if [
            IntersectLocation::OutsideToInside,
            IntersectLocation::Inside,
        ]
        .contains(&prev_intersection.location)
            && point_in_lines(clipping_lines, &joint.d)
        {
            cubics.push(line_cubic_split(
                joint,
                prev_intersection.t,
                1.0,
                &prev_intersection.pos,
                &joint.d,
            ));
        }
    }

    //if we do not have intersections then our shape is either contained (1) or outside (0)
    if lines.is_empty() && quadratics.is_empty() && cubics.is_empty() {
        if clipping_lines.iter().all(|clip_joint| point_in_shape(shape, &aabb, &clip_joint.a, false)) {
            lines.append(&mut clipping_lines.to_vec());
        }
    } else {
        //determine if the clipping polygon vertices are inside the shape we want to clip
        //we also need to include the clipping shape joints where the shape goes outside the clipping shape
        for clip_joint in clipping_lines {
            //we now do the same splitting process with the clipping shape but re-using the intersection tests we have already done
            let mut intersections = clip_joint_intersections.pop_front().unwrap();

            //need to ensure intersections are sorted too (vertex intersections have the same line_t and opposite windings)
            //intersections = intersections.filter(intersection => intersection.t != 0 && intersection.t != 1 && intersection.line_t != 0 && intersection.line_t != 1);
            intersections.sort_by(|a, b| {
                a.line_t
                    .total_cmp(&b.line_t)
                    .then_with(|| a.location.cmp(&b.location))
            });

            intersections = intersections.into_iter().unique_by(|item| ((item.line_t * 10000.0).floor() as u32, item.location)).collect::<Vec<_>>();

            let location = if point_in_shape(shape, &aabb, &clip_joint.a, false) {
                IntersectLocation::Inside
            } else {
                IntersectLocation::Outside
            };

            let mut prev_intersection = &Intersection(0.0, 0.0, clip_joint.a, location);
            for intersection in intersections.iter() {
                if [
                    IntersectLocation::Inside,
                    IntersectLocation::InsideToOutside,
                ]
                .contains(&prev_intersection.location)
                    && [
                        IntersectLocation::Inside,
                        IntersectLocation::OutsideToInside,
                    ]
                    .contains(&intersection.location)
                {
                    lines.push(Line(prev_intersection.pos, intersection.pos));
                }
                prev_intersection = intersection;
            }

            //repeat for the last point in the joint too
            if [
                IntersectLocation::Inside,
                IntersectLocation::InsideToOutside,
            ]
            .contains(&prev_intersection.location)
                && point_in_shape(shape, &aabb, &clip_joint.b, false)
            {
                lines.push(Line(prev_intersection.pos, clip_joint.b));
            }

            //push the transformed intersections back on the end so we can debug
            clip_joint_intersections.push_back(intersections);
        }
    }

    (Shape(lines, quadratics, cubics), clip_joint_intersections)
}
