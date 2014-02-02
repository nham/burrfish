use std::fmt::{Default, Formatter};
use std::f64::{sqrt, cos, sin, abs};
use std::f64::consts::PI;

pub fn deg_to_rad(deg: f64) -> f64 {
    deg * PI / 180.0
}

// relative to the first parameter
pub fn rel_err(a: f64, b: f64) -> f64 {
    abs(a - b) / a
}

#[deriving(Clone)]
pub struct Vector2 {
    x: f64,
    y: f64,
}

impl Add<Vector2, Vector2> for Vector2 {
    fn add(&self, v: &Vector2) -> Vector2 {
        Vector2 { x: self.x + v.x,
                  y: self.y + v.y }

    }
}

impl Sub<Vector2, Vector2> for Vector2 {
    fn sub(&self, v: &Vector2) -> Vector2 {
        Vector2 { x: self.x - v.x,
                  y: self.y - v.y }

    }
}

impl Default for Vector2 {
    fn fmt(v: &Vector2, f: &mut Formatter) {
        write!(f.buf, "({}, {})", v.x, v.y);
    }
}

impl Vector2 {
    pub fn scale(&mut self, scalar: f64) {
        self.x *= scalar;
        self.y *= scalar;
    }

    pub fn scale_copy(&self, scalar: f64) -> Vector2 {
        let mut new = self.clone();
        new.x *= scalar;
        new.y *= scalar;
        new
    }

    pub fn add(&mut self, v: Vector2) {
        self.x += v.x;
        self.y += v.y;
    }

    pub fn dot(&self, v: Vector2) -> f64 {
        (self.x * v.x) + (self.y * v.y)
    }

    pub fn normsq(&self) -> f64 {
        self.dot(*self)
    }

    pub fn norm(&self) -> f64 {
        sqrt(self.normsq())
    }

    pub fn rotate(&mut self, ang: f64) {
        let cosa = cos(ang);
        let sina = sin(ang);
        let x = self.x;
        let y = self.y;
        self.x = cosa * x + sina * y;
        self.y = - sina * x + cosa * y;
    }

    pub fn rotate_copy(&self, ang: f64) -> Vector2 {
        let cosa = cos(ang);
        let sina = sin(ang);
        Vector2 { x: cosa * self.x + sina * self.y,
                  y: - sina * self.x + cosa * self.y }
    }

    // calculates the angle of the vector wrt the x-axis [direction (1, 0)]
    pub fn angx(&self) -> f64 {
        self.dot( Vector2{ x:1.0, y: 1.0} ) / self.norm()
    }

    pub fn zero() -> Vector2 {
        Vector2 { x: 0.0, y: 0.0 }
    }

    pub fn rel_err(&self, v: Vector2) -> f64 {
        (self - v).norm() / self.norm()

    }
}
