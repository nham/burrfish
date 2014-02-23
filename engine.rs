extern crate extra;
extern crate serialize;

use std::fmt::{Default, Formatter};
use std::io;
use std::str;
use std::f64::sqrt;
use std::f64::consts::{SQRT2,PI};
use extra::json;
use serialize::Encodable;

use util::{Vector2, deg_to_rad, rad_to_deg};
use util::{rel_err, negligible_diff};

mod util;

fn main() {
    let mut p1 = Particle { m: 20.0, 
                            pos: Vector2{x:-5.0, y:10.0} };
    let mut p2 = Particle { m: 20.0, 
                            pos: Vector2{x:5.0, y:10.0} };

    // note the angular velocity we start with
    let mut bod = PointBody::new(~[p1, p2], Vector2::zero(), 3.*PI/8.);

    let g = Vector2{x: 0.0, y: -9.8};
    bod.set_force(g);

    simulate(&mut bod, 350, 0.05);

    /////////////

    /*
    let p1 = Particle {m: 5., pos: Vector2{x: -SQRT2/0.2, y: -SQRT2/0.2} };
    let p2 = Particle {m: 5., pos: Vector2{x: 10., y: 0.} };
    let mut bod = PointBody::new(~[p1, p2], Vector2::zero(), 0.);

    let f1 = Vector2{x: 5., y: 0.};
    let f2 = Vector2{x: -5., y: 0.};

    let f1func = |_: f64| -> Vector2 { f1 };
    let f2func = |_: f64| -> Vector2 { f2 };
    let funcs = &[f1func, f2func];

    simulate(&mut bod, 350, 0.05);
    */

    /////////////
    /*
    let s = 30.;
    let p1 = Particle {m: 4., pos: Vector2{x: -s, y: -s} };
    let p2 = Particle {m: 4., pos: Vector2{x: s, y: -s} };
    let p3 = Particle {m: 4., pos: Vector2{x: s, y: s} };
    let p4 = Particle {m: 4., pos: Vector2{x: -s, y: s} };
    let mut bod = PointBody::new(~[p1, p2, p3, p4], Vector2::zero(), 0.);

    let f1 = Vector2{x: 0., y: 4.};
    let f2 = Vector2{x: 0., y: 0.};

    let f1func = |_: f64, ang: f64| -> Vector2 { 
        f1.rotate_copy(ang - 5. * PI / 4.)
    };
    let f2func = |_: f64, _: f64| -> Vector2 { f2 };
    let f3func = |_: f64, _: f64| -> Vector2 { f2 };
    let f4func = |_: f64, _: f64| -> Vector2 { f2 };
    let funcs = &[f1func, f2func, f3func, f4func];

    simulate(&mut bod, 700, 0.05);
    */
}

fn simulate(body: &mut Body, steps: int, dt: f64) {
    // let us simulate
    let mut ds: json::List = ~[];
    ds.push( body.coord_dump_json() );

    for i in range(0, steps) {
        ds.push( euler_step(body, dt)  );
    }
    let report = json::List(ds);

    let mut m = io::MemWriter::new();
    {
        let mut encoder = json::Encoder::new(&mut m);
        report.encode(&mut encoder);
    }

    let z = m.unwrap();
    println!("{:s}", str::from_utf8(z).unwrap() );
}

fn calc_torque(force: Vector2, r: f64, ang: f64) -> f64 {
    let unitx = Vector2 { x: 1.0, y: 0.0 };
    let rot = unitx.rotate_copy(ang + PI/2.);
    r * force.dot( rot )
}

// each step we look at the forces on each point to find the total linear
// acceleration and the total angular acceleration. the former determines
// how the center of mass moves, the latter determines how the object rotates
//
// In order to calculate the angular acceleration, we need to calculate torque,
// which means we need to know the distance from CM to each particle. the fast
// way to do this is store this distance in the body. we shouldn't store the
// complete position of each particle, but merely the relative vector
fn euler_step(body: &mut Body, dt: f64) -> json::Json {
    //let mu = 0.0001; // friction coefficient
    //let friction_force_mag = mu * body.m * 9.8;

    // total linear acceleration

    // this whole section is ewwww. I need some way to do classical inheritance
    // here.
    let lin_acc = body.force()  * body.invmass();
    let new_cmv = body.cmv() + lin_acc * dt;
    body.set_cmv( new_cmv );
    let new_cm = body.cmv() * dt;
    body.set_cm( new_cm );

    let ang_acc = body.torque() * body.invmom();
    let new_angv = body.angv() + ang_acc * dt;
    body.set_angv( new_angv );
    let new_ang = body.ang() + body.angv() * dt;
    body.set_ang( new_ang );

    body.coord_dump_json()
}


struct BaseBody {
    m: f64, // mass
    invm: f64, // inverse of mass. stored here so we don't have to keep recomputing
    mom: f64, // moment of inertia
    invmom: f64, // inverse of moment of inertia

    cm: Vector2, // center of mass
    cmv: Vector2, // linear velocity
    ang: f64, // initially 0, tracks how much the body has rotated since the start
    angv: f64, // angular velocity

    force: Vector2,
    torque: f64,
}

// I really only wanted coord_dump, but theres no way to do classical inheritance
// of member fields, so i have to make the trait be all these getters/setters
trait Body {
    fn coord_dump_json(&self) -> json::Json;
    fn mass(&self) -> f64;
    fn invmass(&self) -> f64;
    fn mom(&self) -> f64;
    fn invmom(&self) -> f64;
    fn cm(&self) -> Vector2;
    fn cmv(&self) -> Vector2;
    fn ang(&self) -> f64;
    fn angv(&self) -> f64;
    fn force(&self) -> Vector2;
    fn torque(&self) -> f64;

    fn set_cm(&mut self, cm: Vector2);
    fn set_cmv(&mut self, cmv: Vector2);
    fn set_ang(&mut self, ang: f64);
    fn set_angv(&mut self, angv: f64);
    fn set_force(&mut self, force: Vector2);
    fn set_torque(&mut self, torque: f64);

    fn linmom(&self) -> Vector2 {
        self.cmv() * self.mass() 
    }

    fn angmom(&self) -> f64 {
        self.mom() * self.angv()
    }
    
}

// The only particle information we store is the mass, the distance from CoM to the
// particle (being a rigid body, this never changes) and the initial angle that the
// position vector from CoM to particle makes with the x-axis. This combined with
// the angle of the body stored in ang allows us to calculate torques at each
// moment without having to explicity update the position vectors of each particle.
struct PointBody {
    body: BaseBody,
    particles: ~[RelParticle],
}

impl PointBody {
    fn coord_dump(&self) -> ~[Vector2] {
        let mut vec: ~[Vector2] = ~[];
        for p in self.particles.iter() {
            let pos = Vector2{x: p.r, y: 0.0}.rotate_copy(self.body.ang + p.init_ang);
            vec.push( self.body.cm + pos );
        }

        vec
    }
}

impl Body for PointBody {
    fn coord_dump_json(&self) -> json::Json {
        let mut report: json::List = ~[];
        for v in self.coord_dump().iter() {
            report.push( json::Number(v.x) );
            report.push( json::Number(v.y) );
        }

        json::List(report)
    }

    fn mass(&self) -> f64 {
        self.body.m
    }
    fn invmass(&self) -> f64 {
        self.body.invm
    }
    fn mom(&self) -> f64 {
        self.body.mom
    }
    fn invmom(&self) -> f64 {
        self.body.invmom
    }
    fn cm(&self) -> Vector2 {
        self.body.cm
    }
    fn cmv(&self) -> Vector2 {
        self.body.cmv
    }
    fn ang(&self) -> f64 {
        self.body.ang
    }
    fn angv(&self) -> f64 {
        self.body.angv
    }
    fn force(&self) -> Vector2 {
        self.body.force
    }
    fn torque(&self) -> f64 {
        self.body.torque
    }

    fn set_cm(&mut self, cm: Vector2) {
        self.body.cm = cm;
    }
    fn set_cmv(&mut self, cmv: Vector2) {
        self.body.cmv = cmv;
    }
    fn set_ang(&mut self, ang: f64) {
        self.body.ang = ang;
    }
    fn set_angv(&mut self, angv: f64) {
        self.body.angv = angv;
    }
    fn set_force(&mut self, force: Vector2) {
        self.body.force = force;
    }
    fn set_torque(&mut self, torque: f64) {
        self.body.torque = torque;
    }
}


struct BoxBody {
    body: BaseBody,
    width: f64,
    height: f64,
}

impl Body for BoxBody {
    fn coord_dump_json(&self) -> json::Json {
        let mut report: json::List = ~[];
        let a = self.width / 2.0;
        let b = self.height / 2.0;
        report.push( json::Number(self.body.cm.x + a));
        report.push( json::Number(self.body.cm.y + b));

        report.push( json::Number(self.body.cm.x - a));
        report.push( json::Number(self.body.cm.y + b));

        report.push( json::Number(self.body.cm.x - a));
        report.push( json::Number(self.body.cm.y - b));

        report.push( json::Number(self.body.cm.x + a));
        report.push( json::Number(self.body.cm.y - b));

        json::List(report)
    }

    fn mass(&self) -> f64 {
        self.body.m
    }
    fn invmass(&self) -> f64 {
        self.body.invm
    }
    fn mom(&self) -> f64 {
        self.body.mom
    }
    fn invmom(&self) -> f64 {
        self.body.invmom
    }
    fn cm(&self) -> Vector2 {
        self.body.cm
    }
    fn cmv(&self) -> Vector2 {
        self.body.cmv
    }
    fn ang(&self) -> f64 {
        self.body.ang
    }
    fn angv(&self) -> f64 {
        self.body.angv
    }
    fn force(&self) -> Vector2 {
        self.body.force
    }
    fn torque(&self) -> f64 {
        self.body.torque
    }

    fn set_cm(&mut self, cm: Vector2) {
        self.body.cm = cm;
    }
    fn set_cmv(&mut self, cmv: Vector2) {
        self.body.cmv = cmv;
    }
    fn set_ang(&mut self, ang: f64) {
        self.body.ang = ang;
    }
    fn set_angv(&mut self, angv: f64) {
        self.body.angv = angv;
    }
    fn set_force(&mut self, force: Vector2) {
        self.body.force = force;
    }
    fn set_torque(&mut self, torque: f64) {
        self.body.torque = torque;
    }
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


// returns (Total Mass, Center of Mass vector)
fn c_of_m(particles: ~[Particle]) -> (f64, Vector2) {
    let mut cm = Vector2::zero();
    let mut sum = 0.0;
    for p in particles.iter() {
        cm.add( p.pos * p.m );
        sum += p.m;
    }

    (sum, cm * ( 1.0 / sum ))
}

// computes the moment of inertia of a set of particles about some point
fn moment_of_inertia(particles: ~[Particle], o: Vector2) -> f64 {
    let mut sum = 0.0;
    for p in particles.iter() {
        sum += p.m * (p.pos - o).normsq();
    }
    sum
}


// two ways we might like to construct a body:
//   1) give a cloud of particles and the initial velocities of the body
//   2) give a center of mass position and, for each particle of the body,
//      its mass and distance/angle from CoM.
impl BaseBody {
    fn new(mass: f64, I: f64, cm: Vector2, cmv: Vector2, ang: f64, angv: f64, force: Vector2, torque: f64) -> BaseBody {
        BaseBody { 
            m: mass,
            invm: 1.0 / mass,
            mom: I,
            invmom: 1.0 / I,
            cm: cm, 
            cmv: cmv,
            ang: ang, 
            angv: angv,
            force: force,
            torque: torque,
        }
    }
}

impl PointBody {
    fn new(ps: ~[Particle], init_cmv: Vector2, init_angv: f64) -> PointBody {
        let (m, cm) = c_of_m(ps.clone());

        let mut mom = 0.0;
        let mut relps: ~[RelParticle] = ~[];

        for p in ps.iter() {
            let radvec = p.pos - cm;
            let rsq = radvec.normsq();
            relps.push(RelParticle { m: p.m,
                                     r: sqrt(rsq),
                                     init_ang: radvec.angx() });
            mom += p.m * rsq;
        }

        let bb = BaseBody { 
                m: m,
                invm: 1.0 / m,
                mom: mom,
                invmom: 1.0 / mom,
                cm: cm, 
                cmv: init_cmv,
                ang: 0.0, 
                angv: init_angv,
                force: Vector2 { x: 0.0, y: 0.0 },
                torque: 0.0,
        };

        PointBody { body: bb, particles: relps }
    }
}

#[test]
fn test_body_new() {
    use std::f64::{sin, cos};

    let p1 = Particle {m: 5., pos: Vector2{x: -SQRT2/0.2, y: -SQRT2/0.2} };
    let p2 = Particle {m: 5., pos: Vector2{x: 10., y: 0.} };
    let b = PointBody::new(~[p1, p2], Vector2::zero(), 0.);
    assert!(b.mass() == 10.);

    let mut v = Vector2{ x: sin(PI/8.), y: -cos(PI/8.) };
    v.scale( 10. * cos(PI/2. - PI/8.) );

    assert!(b.cm().rel_err(v) < 1./100_000_000.);

    let rsq1 = (p1.pos - v).normsq();
    let rsq2 = (p2.pos - v).normsq();
    let real_mom = p1.m * rsq1 + p2.m * rsq2;
    assert!( negligible_diff(b.mom(), real_mom) );

    assert!( b.particles[0].m == p1.m );
    assert!( negligible_diff(b.particles[0].r, (p1.pos - v).norm()) );
    assert!( negligible_diff(b.particles[0].init_ang, PI + PI/8.) );

    assert!( b.particles[1].m == p2.m );
    assert!( negligible_diff(b.particles[1].r, (p2.pos - v).norm()) );
    assert!( negligible_diff(b.particles[1].init_ang, PI/8.) );
}

#[test]
fn test_body_coord_dump() {
    use std::f64::{sin, cos};

    let p1 = Particle {m: 5., pos: Vector2{x: -SQRT2/0.2, y: -SQRT2/0.2} };
    let p2 = Particle {m: 5., pos: Vector2{x: 10., y: 0.} };
    let b = PointBody::new(~[p1, p2], Vector2::zero(), 0.);

    let pos = b.coord_dump();

    assert!( pos[0].rel_err(p1.pos) < 1. / 100_000_000. );
    assert!( pos[1].rel_err(p2.pos) < 1. / 100_000_000. );
}


#[test]
fn test_calc_torque() {
    let f = Vector2 { x: 5., y: 0. };
    let r = 2.;
    let ang = deg_to_rad(225.);

    assert!( negligible_diff(calc_torque(f, r, ang), r * 5. * SQRT2 / 2.) );
}
