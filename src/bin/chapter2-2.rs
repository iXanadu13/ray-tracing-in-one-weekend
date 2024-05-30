use indicatif::ProgressIterator;
use itertools::Itertools;
use std::fs;
use glam::DVec3 as vec3;

const IMAGE_HEIGHT: u32 = 256;
const IMAGE_WIDTH: u32 = 256;
const MAX_VALUE: u8 = 255;

fn main() -> Result<(), Box<dyn std::error::Error>>{
    let pixels = (0..IMAGE_HEIGHT)
        .cartesian_product(0..IMAGE_WIDTH)
        .progress_count(
            IMAGE_HEIGHT as u64 * IMAGE_WIDTH as u64
        )
        .map(|(y, x)|{
            let r = x as f64 / (IMAGE_WIDTH - 1) as f64;
            let g = y as f64 / (IMAGE_HEIGHT - 1) as f64;
            let b = 0.0;
            format!{
                "{} {} {}",
                r * 255.0,
                g * 255.0,
                b * 255.0
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

#[allow(dead_code)]
struct Ray{
    origin: vec3,
    direction: vec3,
}
#[allow(dead_code)]
impl Ray{
    fn at(&self, t: f64) -> vec3 {
        self.origin + t * self.direction
    }
}