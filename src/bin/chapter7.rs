use indicatif::ProgressIterator;
use itertools::Itertools;
use std::{ops::Range, fs};
use glam::DVec3 as vec3;


//Image
const ASPECT_RATIO: f64 = 16.0 / 9.0;
const IMAGE_WIDTH: u32 = 400;


fn main() -> Result<(), Box<dyn std::error::Error>>{
    let camera = Camera::new(IMAGE_WIDTH, ASPECT_RATIO);
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
    camera.render(world)?;

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

#[derive(Debug)]
struct Camera{
    image_width: u32,
    image_height: u32,
    max_value: u8,
    _aspect_ratio: f64,
    center: vec3,
    pixel_delta_u: vec3,
    pixel_delta_v: vec3,
    pixel00_loc: vec3
}
impl Camera{
    fn new(image_width: u32, aspect_ratio: f64) -> Self{
        let max_value: u8 = 255;
        let image_height: u32 = (image_width as f64 / aspect_ratio) as u32;
        // Calculate the image height, and ensure that it's at least 1.
        if image_height < 1 {
            panic!("image height is at least 1");
        }
        let focal_length: f64 = 1.;
        let viewport_height: f64 = 2.0;
        let viewport_width: f64 = viewport_height * (image_width as f64/image_height as f64);
        let center = vec3::ZERO;

        // Calculate the vectors across the horizontal and down the vertical viewport edges.
        let viewport_u: vec3 = vec3::new(viewport_width, 0., 0.);
        let viewport_v: vec3 = vec3::new(0., -viewport_height, 0.);

        // Calculate the horizontal and vertical delta vectors from pixel to pixel.
        let pixel_delta_u: vec3 = viewport_u / image_width as f64;
        let pixel_delta_v: vec3 = viewport_v / image_height as f64;

        // Calculate the location of the upper left pixel.
        let viewport_upper_left: vec3 = center
            - vec3::new(0., 0., focal_length)
            - viewport_u / 2.
            - viewport_v / 2.;
        let pixel00_loc: vec3 = viewport_upper_left
            + 0.5 * (pixel_delta_u + pixel_delta_v);
        Camera{
            image_width,
            image_height,
            max_value,
            _aspect_ratio: aspect_ratio,
            center,
            pixel_delta_u,
            pixel_delta_v,
            pixel00_loc
        }
    }
    fn render<T>(&self, world: T) -> Result<(), Box<dyn std::error::Error>>
    where
        T: Hittable
        //TODO: + Sync
    {
        let pixels = (0..self.image_height)
        .cartesian_product(0..self.image_width)
        .progress_count(
            self.image_height as u64 * self.image_width as u64
        )
        .map(|(y, x)|{
            let pixel_center = self.pixel00_loc
                + (x as f64 * self.pixel_delta_u)
                + (y as f64 * self.pixel_delta_v);
            let ray_direction = 
                pixel_center - self.center;
            let ray = Ray{
                origin: self.center,
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
{0} {1}
{2}
{pixels}
", self.image_width, self.image_height, self.max_value
        ))?;
        Ok(())
    } 
}