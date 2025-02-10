"""
    N-body simulation with local variables instead of globals.
"""

def create_system():
    """Create the initial solar system with constants and body positions."""
    pi = 3.14159265358979323
    solar_mass = 4 * pi * pi
    days_per_year = 365.24
    
    return {
        'sun': ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], solar_mass),

        'jupiter': ([4.84143144246472090e+00,
                    -1.16032004402742839e+00,
                    -1.03622044471123109e-01],
                   [1.66007664274403694e-03 * days_per_year,
                    7.69901118419740425e-03 * days_per_year,
                    -6.90460016972063023e-05 * days_per_year],
                   9.54791938424326609e-04 * solar_mass),

        'saturn': ([8.34336671824457987e+00,
                   4.12479856412430479e+00,
                   -4.03523417114321381e-01],
                  [-2.76742510726862411e-03 * days_per_year,
                   4.99852801234917238e-03 * days_per_year,
                   2.30417297573763929e-05 * days_per_year],
                  2.85885980666130812e-04 * solar_mass),

        'uranus': ([1.28943695621391310e+01,
                   -1.51111514016986312e+01,
                   -2.23307578892655734e-01],
                  [2.96460137564761618e-03 * days_per_year,
                   2.37847173959480950e-03 * days_per_year,
                   -2.96589568540237556e-05 * days_per_year],
                  4.36624404335156298e-05 * solar_mass),

        'neptune': ([1.53796971148509165e+01,
                    -2.59193146099879641e+01,
                    1.79258772950371181e-01],
                   [2.68067772490389322e-03 * days_per_year,
                    1.62824170038242295e-03 * days_per_year,
                    -9.51592254519715870e-05 * days_per_year],
                   5.15138902046611451e-05 * solar_mass)
    }

def compute_deltas(x1, x2, y1, y2, z1, z2):
    return (x1-x2, y1-y2, z1-z2)
    
def compute_b(m, dt, dx, dy, dz):
    mag = compute_mag(dt, dx, dy, dz)
    return m * mag

def compute_mag(dt, dx, dy, dz):
    return dt * ((dx * dx + dy * dy + dz * dz) ** (-1.5))

def update_vs(v1, v2, dt, dx, dy, dz, m1, m2):
    b1 = compute_b(m2, dt, dx, dy, dz)
    b2 = compute_b(m1, dt, dx, dy, dz)
    
    v1[0] -= dx * b1
    v1[1] -= dy * b1
    v1[2] -= dz * b1
    v2[0] += dx * b2
    v2[1] += dy * b2
    v2[2] += dz * b2

def update_rs(r, dt, vx, vy, vz):
    r[0] += dt * vx
    r[1] += dt * vy
    r[2] += dt * vz

def advance(bodies, dt):
    '''
    Advance the system one timestep
    '''
    seenit = []
    for body1 in bodies:
        for body2 in bodies:
            if (body1 != body2) and not (body2 in seenit):
                ([x1, y1, z1], v1, m1) = bodies[body1]
                ([x2, y2, z2], v2, m2) = bodies[body2]
                (dx, dy, dz) = compute_deltas(x1, x2, y1, y2, z1, z2)
                update_vs(v1, v2, dt, dx, dy, dz, m1, m2)
                seenit.append(body1)
        
    for body in bodies:
        (r, [vx, vy, vz], m) = bodies[body]
        update_rs(r, dt, vx, vy, vz)

def compute_energy(m1, m2, dx, dy, dz):
    return (m1 * m2) / ((dx * dx + dy * dy + dz * dz) ** 0.5)
    
def report_energy(bodies):
    '''
    Compute the energy and return it so that it can be printed
    '''
    e = 0.0
    seenit = []
    
    for body1 in bodies:
        for body2 in bodies:
            if (body1 != body2) and not (body2 in seenit):
                ((x1, y1, z1), v1, m1) = bodies[body1]
                ((x2, y2, z2), v2, m2) = bodies[body2]
                (dx, dy, dz) = compute_deltas(x1, x2, y1, y2, z1, z2)
                e -= compute_energy(m1, m2, dx, dy, dz)
                seenit.append(body1)
        
    for body in bodies:
        (r, [vx, vy, vz], m) = bodies[body]
        e += m * (vx * vx + vy * vy + vz * vz) / 2.
        
    return e

def offset_momentum(bodies, ref):
    '''
    Offset momentum relative to the referenced body
    '''
    px = py = pz = 0.0
    
    for body in bodies:
        (r, [vx, vy, vz], m) = bodies[body]
        px -= vx * m
        py -= vy * m
        pz -= vz * m
        
    (r, v, m) = bodies[ref]
    v[0] = px / m
    v[1] = py / m
    v[2] = pz / m

def nbody(loops, reference, iterations):
    '''
    Run nbody simulation
    loops - number of loops to run
    reference - body at center of system
    iterations - number of timesteps to advance
    '''
    bodies = create_system()
    
    # Set up initial conditions
    offset_momentum(bodies, reference)

    for _ in range(loops):
        energy = report_energy(bodies)
        for _ in range(iterations):
            advance(bodies, 0.01)
        print(report_energy(bodies))

if __name__ == '__main__':
    nbody(100, 'sun', 20000)
