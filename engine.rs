use std::fmt::{Default, Formatter};

fn main() {
    let mut p1 = Particle { m: 1.0, 
                            pos: Vector2{x:-5.0, y:10.0},
                            vel: Vector2{x:0.0, y:0.0} };
    let mut p2 = Particle { m: 20.0, 
                            pos: Vector2{x:5.0, y:10.0},
                            vel: Vector2{x:0.0, y:0.0} };


    println!("P1: {}\nP2: {}", p1, p2);
    println!("--------");

    let dt = 0.01;
    let g = Vector2{x: 0.0, y: -9.8};

    euler_step(&mut p1, 0.0, dt, |_: f64| -> Vector2 { g });

    println!("P1: {}\nP2: {}", p1, p2);

}

fn euler_step(p: &mut Particle, t: f64, dt: f64, force: |f64| -> Vector2) {
    let acc = force(t).scale_copy(1.0/p.m);
    let delta_v = acc.scale_copy(dt);
    p.vel.add( delta_v );
    p.pos.add( p.vel.scale_copy(dt) );
    
}


struct Body {
    m: f64,
    ps: ~[Particle],
}

struct Particle {
    m: f64,
    pos: Vector2,
    vel: Vector2,
}

impl Default for Particle {
    fn fmt(p: &Particle, f: &mut Formatter) {
        write!(f.buf, "(m: {}, pos: {}, vel: {})", p.m, p.pos, p.vel);
    }
}

impl Particle {
    fn linmom(&self) -> Vector2 {
        self.vel.scale_copy( self.m )
    }

    fn angmom(&self, o: Vector2) -> f64 {
        self.linmom().dot( self.pos.sub_copy( o ) )
        
    }
}


impl Body {
    fn new(ps: ~[Particle]) -> Body {
        let mut sum = 0.0;
        for p in ps.iter() {
           sum += p.m;
        }

        Body { m: sum, ps: ps }
    }

    fn linmom(&self) -> Vector2 {
        self.c_of_m().scale_copy( self.m )
    }

    fn c_of_m(&self) -> Vector2 {
        let mut cm = Vector2{ x: 0.0, y: 0.0};
        for p in self.ps.iter() {
            cm.add( p.linmom() );
        }

        cm.scale_copy( 1.0 / self.m )
    }
}



#[deriving(Clone)]
struct Vector2 {
    x: f64,
    y: f64,
}

impl Default for Vector2 {
    fn fmt(v: &Vector2, f: &mut Formatter) {
        write!(f.buf, "({}, {})", v.x, v.y);
    }
}

impl Vector2 {
    fn scale(&mut self, scalar: f64) {
        self.x *= scalar;
        self.y *= scalar;
    }

    fn scale_copy(&self, scalar: f64) -> Vector2 {
        let mut new = self.clone();
        new.x *= scalar;
        new.y *= scalar;
        new
    }

    fn add(&mut self, v: Vector2) {
        self.x += v.x;
        self.y += v.y;
    }

    fn add_copy(& self, v: Vector2) -> Vector2 {
        let mut new = self.clone();
        new.x += v.x;
        new.y += v.y;
        new
    }

    fn sub_copy(& self, v: Vector2) -> Vector2 {
        self.add_copy( v.scale_copy(-1.0) )
    }

    fn dot(&self, v: Vector2) -> f64 {
        (self.x * v.x) + (self.y * v.y)
    }

    fn norm(&self) -> f64 {
        self.dot(*self)
    }
}
