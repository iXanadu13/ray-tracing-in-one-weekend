use indicatif::ProgressIterator;
use itertools::Itertools;
use std::fs;
use glam::DVec3 as vec3;

const MAX_VALUE: u8 = 255;

const ASPECT_RATIO: f64 = 16.0 / 9.0;
const IMAGE_WIDTH: u32 = 400;
const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;

const VIEWPORT_HEIGHT: f64 = 2.0;
const VIEWPORT_WIDTH: f64 = VIEWPORT_HEIGHT * (IMAGE_WIDTH as f64/IMAGE_HEIGHT as f64);

// Calculate the vectors across the horizontal and down the vertical viewport edges.
const VIEWPORT_U: vec3 = 
    vec3::new(VIEWPORT_WIDTH, 0., 0.);
const VIEWPORT_V: vec3 = 
    vec3::new(0., -VIEWPORT_HEIGHT, 0.);

const CAMERA_CENTER:vec3 = vec3::new(0.,0.,0.);
const FOCAL_LENGTH: f64 = 1.;

fn main() -> Result<(), Box<dyn std::error::Error>>{
    // Calculate the image height, and ensure that it's at least 1.
    if IMAGE_HEIGHT < 1 {
        panic!("image height is at least 1");
    }
    // Calculate the horizontal and vertical delta vectors from pixel to pixel.
    let pixel_delta_u: vec3 = 
        VIEWPORT_U / IMAGE_WIDTH as f64;
    let pixel_delta_v: vec3 = 
        VIEWPORT_V / IMAGE_HEIGHT as f64;
    // Calculate the location of the upper left pixel.
    let viewport_upper_left: vec3 = CAMERA_CENTER
        - vec3::new(0., 0., FOCAL_LENGTH)
        - VIEWPORT_U / 2.
        - VIEWPORT_V / 2.;
    let pixel00_loc: vec3 = viewport_upper_left
        + 0.5 * (pixel_delta_u + pixel_delta_v);

    let pixels = (0..IMAGE_HEIGHT)
        .cartesian_product(0..IMAGE_WIDTH)
        .progress_count(
            IMAGE_HEIGHT as u64 * IMAGE_WIDTH as u64
        )
        .map(|(y, x)|{
            let pixel_center = pixel00_loc
                + (x as f64 * pixel_delta_u)
                + (y as f64 * pixel_delta_v);
            let ray_direction = 
                pixel_center - CAMERA_CENTER;
            let ray = Ray{
                origin: CAMERA_CENTER,
                direction: ray_direction
            };
            let pixel_color = ray.color() * 255.0;
            format!{
                "{} {} {}",
                pixel_color.x, pixel_color.y, pixel_color.z
            }
        })
        .join("\n");
    
    fs::write("output.ppm", format!(
"P3
{IMAGE_WIDTH} {IMAGE_HEIGHT}
{MAX_VALUE}
{pixels}
"
    ))?;

    Ok(())
}

struct Ray{
    origin: vec3,
    direction: vec3,
}
impl Ray{
    fn at(&self, t: f64) -> vec3 {
        self.origin + t * self.direction
    }
    #[allow(non_snake_case)]
    fn color(&self) -> vec3 {
        let t = hit_sphere(&vec3::new(0., 0., -1.), 0.5, self);
        if t > 0.0 {
            let N = (self.at(t) - vec3::new(0., 0., -1.)).normalize();
            return 0.5 * (N + 1.0);
        }
        let unit_direction: vec3 = 
            self.direction.normalize();
        let a = 0.5 * (unit_direction.y + 1.0);
        // 线性插值
        // blendedValue=(1−a)⋅startValue+a⋅endValue
        return (1.0 - a) * vec3::new(1.0, 1.0, 1.0)
            + a * vec3::new(0.5, 0.7, 1.0);
    }
}

fn hit_sphere(center: &vec3, radius: f64, ray: &Ray) -> f64 {
    let oc: vec3 = ray.origin - *center;  // we use Q - C instead
    let a = ray.direction.length_squared();
    let half_b = ray.direction.dot(oc);
    let c = oc.length_squared() - radius * radius;
    let discriminant = half_b * half_b - a * c;
    if discriminant < 0. {
        -1.
    }
    else{
        (- half_b - discriminant.sqrt()) / a
    }
}