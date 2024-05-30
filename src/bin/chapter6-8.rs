use indicatif::ProgressIterator;
use itertools::Itertools;
use std::{ops::Range, fs};
use glam::DVec3 as vec3;

const MAX_VALUE: u8 = 255;

//Image
const ASPECT_RATIO: f64 = 16.0 / 9.0;
const IMAGE_WIDTH: u32 = 400;
const IMAGE_HEIGHT: u32 = (IMAGE_WIDTH as f64 / ASPECT_RATIO) as u32;

//Camera
const FOCAL_LENGTH: f64 = 1.;
const VIEWPORT_HEIGHT: f64 = 2.0;
const VIEWPORT_WIDTH: f64 = VIEWPORT_HEIGHT * (IMAGE_WIDTH as f64/IMAGE_HEIGHT as f64);
const CAMERA_CENTER:vec3 = vec3::ZERO;

// Calculate the vectors across the horizontal and down the vertical viewport edges.
const VIEWPORT_U: vec3 = 
    vec3::new(VIEWPORT_WIDTH, 0., 0.);
const VIEWPORT_V: vec3 = 
    vec3::new(0., -VIEWPORT_HEIGHT, 0.);


fn main() -> Result<(), Box<dyn std::error::Error>>{
    // Calculate the image height, and ensure that it's at least 1.
    if IMAGE_HEIGHT < 1 {
        panic!("image height is at least 1");
    }
    //World
    let mut world = HittableList::new();
    world.add(Sphere{
        center: vec3::new(0.0, 0.0, -1.0),
        radius: 0.5
    });
    world.add(Sphere{
        center: vec3::new(0.0, -100.5, -1.0),
        radius: 100.0
    });
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

    // Render
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
            let pixel_color = ray.color(&world) * 255.0;
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
    fn color<T>(&self, world: &T) -> vec3 
    where
        T: Hittable
    {
        if let Some(rec) = 
            world.hit(&self, (0.)..f64::INFINITY)
        {
            return 0.5 * (rec.normal + vec3::new(1., 1., 1.));
        }
        // let t = hit_sphere(&vec3::new(0., 0., -1.), 0.5, self);
        // if t > 0.0 {
        //     let N = (self.at(t) - vec3::new(0., 0., -1.)).normalize();
        //     return 0.5 * (N + 1.0);
        // }
        let unit_direction: vec3 = 
            self.direction.normalize();
        let a = 0.5 * (unit_direction.y + 1.0);
        // 线性插值
        // blendedValue=(1−a)⋅startValue+a⋅endValue
        return (1.0 - a) * vec3::new(1.0, 1.0, 1.0)
            + a * vec3::new(0.5, 0.7, 1.0);
    }
}

struct Sphere{
    center: vec3,
    radius: f64
}
impl Hittable for Sphere{
    fn hit(
        &self,
        ray: &Ray,
        interval: Range<f64>
    ) -> Option<HitRecord> {
        let oc: vec3 = ray.origin - self.center;  // we use Q - C instead
        let a = ray.direction.length_squared();
        let half_b = ray.direction.dot(oc);
        let c = oc.length_squared() - self.radius * self.radius;
        let discriminant = half_b * half_b - a * c;
        if discriminant < 0. {
            return None;
        }
        let sqrtd = discriminant.sqrt();

        // Find the nearest root that lies in the acceptable range.
        let mut root = (-half_b - sqrtd) / a;
        if !interval.contains(&root) {
            root = (-half_b + sqrtd) / a;
            if !interval.contains(&root) {
                return None;
            }
        }
        let t = root;
        let point = ray.at(t);
        let outward_normal = (point - self.center) / self.radius;
        let rec = HitRecord::with_face_normal(
            point,
            outward_normal,
            t,
            ray
        );
        Some(rec)

    }
}


struct HitRecord{
    _point: vec3,
    normal: vec3,
    t: f64,
    _front_face: bool
}
impl HitRecord{
    fn calc_face_normal(outward_normal: &vec3, ray: &Ray) -> (bool, vec3){
        let front_face = ray.direction.dot(*outward_normal) < 0.;
        let normal = if front_face {
            *outward_normal
        }
        else{
            - *outward_normal
        };
        (front_face, normal)
    }
    fn with_face_normal(
        point: vec3,
        outward_normal: vec3,
        t: f64,
        ray: &Ray
    ) -> Self{
        let (front_face, normal) = HitRecord::calc_face_normal(&outward_normal, ray);
        HitRecord{
            _point: point,
            normal,
            t,
            _front_face: front_face
        }
    }
}

trait Hittable {
    fn hit(
        &self,
        ray: &Ray,
        interval: Range<f64>
    ) -> Option<HitRecord>;
}

struct HittableList{
    objects: Vec<Box<dyn Hittable>>
}
impl HittableList{
    fn new() -> Self{
        HittableList{
            objects: vec![]
        }
    }
    #[allow(dead_code)]
    fn clear(&mut self){
        self.objects = vec![];
    }
    fn add<T>(&mut self, object: T)
    where
        T: Hittable + 'static,
    {
        self.objects.push(Box::new(object));
    }
}

impl Hittable for HittableList{
    fn hit(
        &self,
        ray: &Ray,
        interval: Range<f64>
    ) -> Option<HitRecord> {
        let (_closest, hit_record) = self
            .objects
            .iter()
            .fold((interval.end, None), |acc, item|{
                if let Some(temp_rec) = item.hit(
                    ray,
                    interval.start..acc.0,
                ){
                    (temp_rec.t, Some(temp_rec))
                }
                else {
                    acc
                }
            });
        hit_record
    }
}

