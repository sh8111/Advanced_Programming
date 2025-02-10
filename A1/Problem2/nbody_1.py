"""
    N-body simulation with optimized function calls.
"""

PI = 3.14159265358979323
SOLAR_MASS = 4 * PI * PI
DAYS_PER_YEAR = 365.24

BODIES = {
    'sun': ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], SOLAR_MASS),

    'jupiter': ([4.84143144246472090e+00,
                 -1.16032004402742839e+00,
                 -1.03622044471123109e-01],
                [1.66007664274403694e-03 * DAYS_PER_YEAR,
                 7.69901118419740425e-03 * DAYS_PER_YEAR,
                 -6.90460016972063023e-05 * DAYS_PER_YEAR],
                9.54791938424326609e-04 * SOLAR_MASS),

    'saturn': ([8.34336671824457987e+00,
                4.12479856412430479e+00,
                -4.03523417114321381e-01],
               [-2.76742510726862411e-03 * DAYS_PER_YEAR,
                4.99852801234917238e-03 * DAYS_PER_YEAR,
                2.30417297573763929e-05 * DAYS_PER_YEAR],
               2.85885980666130812e-04 * SOLAR_MASS),

    'uranus': ([1.28943695621391310e+01,
                -1.51111514016986312e+01,
                -2.23307578892655734e-01],
               [2.96460137564761618e-03 * DAYS_PER_YEAR,
                2.37847173959480950e-03 * DAYS_PER_YEAR,
                -2.96589568540237556e-05 * DAYS_PER_YEAR],
               4.36624404335156298e-05 * SOLAR_MASS),

    'neptune': ([1.53796971148509165e+01,
                 -2.59193146099879641e+01,
                 1.79258772950371181e-01],
                [2.68067772490389322e-03 * DAYS_PER_YEAR,
                 1.62824170038242295e-03 * DAYS_PER_YEAR,
                 -9.51592254519715870e-05 * DAYS_PER_YEAR],
                5.15138902046611451e-05 * SOLAR_MASS)}

def advance(dt):
    '''
    Advance the system one timestep - optimized version with reduced function calls
    '''
    seenit = []
    for body1 in BODIES.keys():
        for body2 in BODIES.keys():
            if (body1 != body2) and not (body2 in seenit):
                ([x1, y1, z1], v1, m1) = BODIES[body1]
                ([x2, y2, z2], v2, m2) = BODIES[body2]
                
                # Combined compute_deltas inline
                dx = x1 - x2
                dy = y1 - y2
                dz = z1 - z2
                
                # Combined compute_mag and compute_b calculations inline
                dist_squared = dx * dx + dy * dy + dz * dz
                mag = dt * (dist_squared ** (-1.5))
                
                b1 = m2 * mag
                b2 = m1 * mag
                
                # Combined update_vs calculations inline
                v1[0] -= dx * b1
                v1[1] -= dy * b1
                v1[2] -= dz * b1
                v2[0] += dx * b2
                v2[1] += dy * b2
                v2[2] += dz * b2
                
                seenit.append(body1)
        
    # Combined update_rs calculations inline
    for body in BODIES.keys():
        (r, [vx, vy, vz], m) = BODIES[body]
        r[0] += dt * vx
        r[1] += dt * vy
        r[2] += dt * vz

def report_energy():
    '''
    Compute the energy and return it so that it can be printed - optimized version
    '''
    e = 0.0
    seenit = []
    
    for body1 in BODIES.keys():
        for body2 in BODIES.keys():
            if (body1 != body2) and not (body2 in seenit):
                ((x1, y1, z1), v1, m1) = BODIES[body1]
                ((x2, y2, z2), v2, m2) = BODIES[body2]
                
                # Combined compute_deltas and compute_energy calculations inline
                dx = x1 - x2
                dy = y1 - y2
                dz = z1 - z2
                
                # Simplified energy calculation
                e -= (m1 * m2) / ((dx * dx + dy * dy + dz * dz) ** 0.5)
                seenit.append(body1)
        
    # Combined kinetic energy calculations
    for body in BODIES.keys():
        (r, [vx, vy, vz], m) = BODIES[body]
        e += m * (vx * vx + vy * vy + vz * vz) / 2
        
    return e

def offset_momentum(ref):
    '''
    Offset momentum relative to the referenced body - optimized version
    '''
    px = py = pz = 0.0
    
    for body in BODIES.keys():
        (r, [vx, vy, vz], m) = BODIES[body]
        px -= vx * m
        py -= vy * m
        pz -= vz * m
        
    (r, v, m) = ref
    v[0] = px / m
    v[1] = py / m
    v[2] = pz / m

def nbody(loops, reference, iterations):
    '''
    nbody simulation - optimized version
    '''
    # Set up global state
    offset_momentum(BODIES[reference])

    for _ in range(loops):
        report_energy()
        for _ in range(iterations):
            advance(0.01)
        print(report_energy())

if __name__ == '__main__':
    nbody(100, 'sun', 20000)
