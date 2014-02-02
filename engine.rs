use std::fmt::{Default, Formatter};
use std::f64::sqrt;

use util::Vector2;

fn main() {
    let mut p1 = Particle { m: 20.0, 
                            pos: Vector2{x:-5.0, y:10.0} };
    let mut p2 = Particle { m: 20.0, 
                            pos: Vector2{x:5.0, y:10.0} };

    let mut bod = Body::new(~[p1, p2], 0.0, Vector2 { x: 0.0, y: 0.0 }, 1.0);

    println!("P1: {}\nP2: {}", p1, p2);
    println!("--------");

    let dt = 0.01;
    let g = Vector2{x: 0.0, y: -9.8};

 //   euler_step(&mut p1, 0.0, dt, |_: f64| -> Vector2 { g });

    println!("P1: {}\nP2: {}", p1, p2);

}

// each step we look at the forces on each point to find the total linear
// acceleration and the total angular acceleration. the former determines
// how the center of mass moves, the latter determines how the object rotates
//
// In order to calculate the angular acceleration, we need to calculate torque,
// which means we need to know the distance from CM to each particle. the fast
// way to do this is store this distance in the body. we shouldn't store the
// complete position of each particle, but merely the relative vector
fn euler_step(body: &mut Body, t: f64, dt: f64, force: ~[|f64| -> Force]) {
    let mut tot_force = Vector2 {x: 0.0, y: 0.0};
    let mut tot_torque = 0.0;
    for (f, p) in force.iter().zip( body.particles.iter() ) {
        tot_force.add( (*f)(t) );
        
        let ang = body.ang - p.init_ang;
        let unitx = Vector2 { x: 1.0, y: 0.0 };
        tot_torque += (*f)(t).dot( unitx.rotate_copy(ang) );

    }

    // total linear acceleration
    let lin_acc = tot_force.scale_copy(1.0/body.m);
    let delta_v = lin_acc.scale_copy(dt);

    body.linv.add( delta_v );
    body.cm.add( body.linv.scale_copy(dt) );

    let ang_acc = tot_torque / body.mom;

    body.angv += ang_acc * dt;
    body.ang += body.angv * dt;
}

type Force = Vector2;

// The only particle information we store is the mass, the distance from CoM to the
// particle (being a rigid body, this never changes) and the initial angle that the
// position vector from CoM to particle makes with the x-axis. This combined with
// the angle of the body stored in ang allows us to calculate torques at each
// moment without having to explicity update the position vectors of each particle.
struct Body {
    particles: ~[RelParticle],
    m: f64, // mass
    mom: f64, // moment of inertia
    ang: f64, // angle
    angv: f64, // angular velocity
    cm: Vector2, // center of mass
    linv: Vector2, // linear velocity
}

struct RelParticle {
    m: f64, // mass
    r: f64, // distance from CoM
    init_ang: f64, // initial angle
}

// Particles used to track their own velocities, but I removed them because generally
// we will be dealing Particles as components of rigid bodies, where we do not need
// to track the particle velocities individually.
//
// To recover the ability to model a single particle, we can create a Body with a
// single particle in it. It's a bit of a hack, but does work.
#[deriving(Clone)]
struct Particle {
    m: f64,
    pos: Vector2,
}

impl Default for Particle {
    fn fmt(p: &Particle, f: &mut Formatter) {
        write!(f.buf, "(m: {}, pos: {})", p.m, p.pos);
    }
}


// returns (Total Mass, Center of Mass vector)
fn c_of_m(particles: ~[Particle]) -> (f64, Vector2) {
    let mut cm = Vector2{ x: 0.0, y: 0.0};
    let mut sum = 0.0;
    for p in particles.iter() {
        cm.add( p.pos.scale_copy(p.m) );
        sum += p.m;
    }

    (sum, cm.scale_copy( 1.0 / sum ))
}

// computes the moment of inertia of a set of particles about some point
fn moment_of_inertia(particles: ~[Particle], o: Vector2) -> f64 {
    let mut sum = 0.0;
    for p in particles.iter() {
        sum += p.m * p.pos.sub_copy(o).normsq();
    }
    sum
}


// two ways we might like to construct a body:
//   1) give a cloud of particles and the initial velocities of the body
//   2) give a center of mass position and, for each particle of the body,
//      its mass and distance/angle from CoM.
impl Body {
    // init_ang here serves to rotate the collection of particles (the given
    // particle positions are taken to be angle = 0)
    fn new(ps: ~[Particle], init_ang: f64, init_linv: Vector2, init_angv: f64) -> Body {
        let (m, cm) = c_of_m(ps.clone());

        let mut mom = 0.0;
        let mut relps: ~[RelParticle] = ~[];

        for p in ps.iter() {
            let rsq = p.pos.sub_copy(cm).normsq();
            relps.push(RelParticle { m: p.m, 
                                     r: sqrt(rsq), 
                                     init_ang: p.pos.angx() + init_ang});
            mom += p.m * rsq;
        }

        Body { 
            m: m,
            mom: mom,
            cm: cm, 
            linv: init_linv,
            ang: 0.0, 
            angv: init_angv,
            particles: relps,
        }
    }


    fn linmom(&self) -> Vector2 {
        self.linv.scale_copy( self.m )
    }

    fn angmom(&self) -> f64 {
        self.mom * self.angv
    }

}



mod util {
    use std::fmt::{Default, Formatter};
    use std::f64::{sqrt, cos, sin};

    #[deriving(Clone)]
    pub struct Vector2 {
        x: f64,
        y: f64,
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

        pub fn add_copy(& self, v: Vector2) -> Vector2 {
            let mut new = self.clone();
            new.x += v.x;
            new.y += v.y;
            new
        }

        pub fn sub_copy(& self, v: Vector2) -> Vector2 {
            self.add_copy( v.scale_copy(-1.0) )
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
    }
}
