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
    /*
    let mut p1 = Particle { m: 20.0, 
                            pos: Vector2{x:-5.0, y:10.0} };
    let mut p2 = Particle { m: 20.0, 
                            pos: Vector2{x:5.0, y:10.0} };

    // note the angular velocity we start with
    let mut bod = Body::new_by_particles(~[p1, p2], Vector2::zero(), 3.*PI/8.);

    let dt = 0.05;
    let steps = 350;
    let g = Vector2{x: 0.0, y: -9.8};

    let gfunc = |_: f64| -> Vector2 { g };
    let gfunc1 = |_: f64| -> Vector2 { g };
    let funcs = &[gfunc, gfunc1];
    */

    /////////////

    /*
    let p1 = Particle {m: 5., pos: Vector2{x: -SQRT2/0.2, y: -SQRT2/0.2} };
    let p2 = Particle {m: 5., pos: Vector2{x: 10., y: 0.} };
    let mut bod = Body::new_by_particles(~[p1, p2], Vector2::zero(), 0.);

    let dt = 0.05;
    let steps = 350;

    let f1 = Vector2{x: 5., y: 0.};
    let f2 = Vector2{x: -5., y: 0.};

    let f1func = |_: f64| -> Vector2 { f1 };
    let f2func = |_: f64| -> Vector2 { f2 };
    let funcs = &[f1func, f2func];
    */

    /////////////
    let s = 30.;
    let p1 = Particle {m: 4., pos: Vector2{x: -s, y: -s} };
    let p2 = Particle {m: 4., pos: Vector2{x: s, y: -s} };
    let p3 = Particle {m: 4., pos: Vector2{x: s, y: s} };
    let p4 = Particle {m: 4., pos: Vector2{x: -s, y: s} };
    let mut bod = Body::new_by_particles(~[p1, p2, p3, p4], Vector2::zero(), 0.);

    let dt = 0.05;
    let steps = 700;

    let f1 = Vector2{x: 0., y: 4.};
    let f2 = Vector2{x: 0., y: 0.};

    let f1func = |_: f64, ang: f64| -> Vector2 { 
        f1.rotate_copy(ang - 5. * PI / 4.)
    };
    let f2func = |_: f64, _: f64| -> Vector2 { f2 };
    let f3func = |_: f64, _: f64| -> Vector2 { f2 };
    let f4func = |_: f64, _: f64| -> Vector2 { f2 };
    let funcs = &[f1func, f2func, f3func, f4func];


    // let us simulate
    let mut ds: json::List = ~[];
    ds.push( bod.pos_dump_json() );
    for i in range(0, steps) {
        ds.push( euler_step(&mut bod, 0.0, dt, funcs) );
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
fn euler_step(body: &mut Body, t: f64, dt: f64, force: &[|f64, f64| -> Force]) -> json::Json {
    let mu = 0.0001; // friction coefficient
    let friction_force_mag = mu * body.m * 9.8;

    let mut tot_force = Vector2::zero();
    let mut tot_torque = 0.0;

    let vels = body.vel_dump();
    let mut vel_iter = vels.iter();

    for (f, p) in force.iter().zip( body.particles.iter() ) {
        let ang = body.ang + p.init_ang;
        let this_force = (*f)(t, ang);

        let this_vel = vel_iter.next().unwrap();
        if this_vel.is_approx_zero() {
            tot_force.add( this_force );
        } else {
            let mut frforce = this_vel.normalize();
            frforce.scale(-friction_force_mag);
            tot_force.add( this_force + frforce );
        }
        
        let this_torque = calc_torque(this_force, p.r, ang);
        debug!("this_torque = {}", this_torque);
        tot_torque += this_torque;

    }

    debug!("total torque = {}", tot_torque);

    // total linear acceleration
    let lin_acc = tot_force.scale_copy(body.invm);
    let delta_v = lin_acc.scale_copy(dt);

    body.linv.add( delta_v );
    body.cm.add( body.linv.scale_copy(dt) );

    let ang_acc = tot_torque * body.invmom;

    body.angv += ang_acc * dt;
    body.ang += body.angv * dt;

    body.pos_dump_json()
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
    invm: f64, // inverse of mass. stored here so we don't have to keep recomputing
    mom: f64, // moment of inertia
    invmom: f64, // inverse of moment of inertia

    cm: Vector2, // center of mass
    linv: Vector2, // linear velocity
    ang: f64, // initially 0, tracks how much the body has rotated since the start
    angv: f64, // angular velocity
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
        cm.add( p.pos.scale_copy(p.m) );
        sum += p.m;
    }

    (sum, cm.scale_copy( 1.0 / sum ))
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
impl Body {
    fn new(ps: ~[RelParticle], mass: f64, I: f64, cm: Vector2, linv: Vector2, angv: f64) -> Body {
        Body { 
            m: mass,
            invm: 1.0 / mass,
            mom: I,
            invmom: 1.0 / I,
            cm: cm, 
            linv: linv,
            ang: 0.0, 
            angv: angv,
            particles: ps,
        }

    }

    fn new_by_particles(ps: ~[Particle], init_linv: Vector2, init_angv: f64) -> Body {
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

        Body { 
            m: m,
            invm: 1.0 / m,
            mom: mom,
            invmom: 1.0 / mom,
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

    // dump the actual position vectors of each particle
    fn pos_dump(&self) -> ~[Vector2] {
        let mut vec: ~[Vector2] = ~[];
        for p in self.particles.iter() {
            let pos = Vector2{x: p.r, y: 0.0}.rotate_copy(self.ang + p.init_ang);
            vec.push( self.cm + pos );
        }

        vec
    }

    fn pos_dump_json(&self) -> json::Json {
        let mut report: json::List = ~[];
        for v in self.pos_dump().iter() {
            report.push( json::Number(v.x) );
            report.push( json::Number(v.y) );
        }

        json::List(report)
    }

    fn vel_dump(&self) -> ~[Vector2] {
        let mut vec: ~[Vector2] = ~[];
        for p in self.particles.iter() {
            let ang = self.ang + p.init_ang;

            let mut rotv = (Vector2 { x: 1., y: 0. });
            rotv.rotate_copy(ang + PI/2.);
            rotv.scale(self.angv);

            vec.push( rotv + self.linv );
        }

        vec
    }

}

#[test]
fn test_body_new() {
    use std::f64::{sin, cos};

    let p1 = Particle {m: 5., pos: Vector2{x: -SQRT2/0.2, y: -SQRT2/0.2} };
    let p2 = Particle {m: 5., pos: Vector2{x: 10., y: 0.} };
    let b = Body::new_by_particles(~[p1, p2], Vector2::zero(), 0.);
    assert!(b.m == 10.);

    let mut v = Vector2{ x: sin(PI/8.), y: -cos(PI/8.) };
    v.scale( 10. * cos(PI/2. - PI/8.) );

    assert!(b.cm.rel_err(v) < 1./100_000_000.);

    let rsq1 = (p1.pos - v).normsq();
    let rsq2 = (p2.pos - v).normsq();
    let real_mom = p1.m * rsq1 + p2.m * rsq2;
    assert!( negligible_diff(b.mom, real_mom) );

    assert!( b.particles[0].m == p1.m );
    assert!( negligible_diff(b.particles[0].r, (p1.pos - v).norm()) );
    assert!( negligible_diff(b.particles[0].init_ang, PI + PI/8.) );

    assert!( b.particles[1].m == p2.m );
    assert!( negligible_diff(b.particles[1].r, (p2.pos - v).norm()) );
    assert!( negligible_diff(b.particles[1].init_ang, PI/8.) );
}

#[test]
fn test_body_pos_dump() {
    use std::f64::{sin, cos};

    let p1 = Particle {m: 5., pos: Vector2{x: -SQRT2/0.2, y: -SQRT2/0.2} };
    let p2 = Particle {m: 5., pos: Vector2{x: 10., y: 0.} };
    let b = Body::new_by_particles(~[p1, p2], Vector2::zero(), 0.);

    let pos = b.pos_dump();

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
