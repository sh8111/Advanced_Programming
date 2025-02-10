"""
    Optimized N-body simulation with reduced function calls, better membership testing,
    local variables, and reduced loop overhead through data aggregation.
"""

def main():
    # Local constants
    PI = 3.14159265358979323
    SOLAR_MASS = 4 * PI * PI
    DAYS_PER_YEAR = 365.24

    # Use tuple for immutable data and local dictionary
    BODIES = {
        'sun': ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], SOLAR_MASS),
        'jupiter': ([4.84143144246472090e+00, -1.16032004402742839e+00, -1.03622044471123109e-01],
                   [1.66007664274403694e-03 * DAYS_PER_YEAR, 7.69901118419740425e-03 * DAYS_PER_YEAR, -6.90460016972063023e-05 * DAYS_PER_YEAR],
                   9.54791938424326609e-04 * SOLAR_MASS),
        'saturn': ([8.34336671824457987e+00, 4.12479856412430479e+00, -4.03523417114321381e-01],
                  [-2.76742510726862411e-03 * DAYS_PER_YEAR, 4.99852801234917238e-03 * DAYS_PER_YEAR, 2.30417297573763929e-05 * DAYS_PER_YEAR],
                  2.85885980666130812e-04 * SOLAR_MASS),
        'uranus': ([1.28943695621391310e+01, -1.51111514016986312e+01, -2.23307578892655734e-01],
                  [2.96460137564761618e-03 * DAYS_PER_YEAR, 2.37847173959480950e-03 * DAYS_PER_YEAR, -2.96589568540237556e-05 * DAYS_PER_YEAR],
                  4.36624404335156298e-05 * SOLAR_MASS),
        'neptune': ([1.53796971148509165e+01, -2.59193146099879641e+01, 1.79258772950371181e-01],
                   [2.68067772490389322e-03 * DAYS_PER_YEAR, 1.62824170038242295e-03 * DAYS_PER_YEAR, -9.51592254519715870e-05 * DAYS_PER_YEAR],
                   5.15138902046611451e-05 * SOLAR_MASS)
    }

    # Pre-compute body pairs to reduce loop overhead
    BODY_PAIRS = [(b1, b2) for i, b1 in enumerate(BODIES.keys()) 
                          for b2 in list(BODIES.keys())[i + 1:]]

    def advance(dt, bodies=BODIES, pairs=BODY_PAIRS):
        """Advance the system one timestep with reduced function calls and loop overhead"""
        for body1, body2 in pairs:
            ([x1, y1, z1], v1, m1) = bodies[body1]
            ([x2, y2, z2], v2, m2) = bodies[body2]
            
            # Inline compute_deltas
            dx = x1 - x2
            dy = y1 - y2
            dz = z1 - z2
            
            # Inline compute_mag and compute_b
            mag = dt * ((dx * dx + dy * dy + dz * dz) ** (-1.5))
            b1 = m2 * mag
            b2 = m1 * mag

            # Inline update_vs
            v1[0] -= dx * b1
            v1[1] -= dy * b1
            v1[2] -= dz * b1
            v2[0] += dx * b2
            v2[1] += dy * b2
            v2[2] += dz * b2

        # Inline update_rs
        for body in bodies.values():
            r, [vx, vy, vz], m = body
            r[0] += dt * vx
            r[1] += dt * vy
            r[2] += dt * vz

    def report_energy(bodies=BODIES, pairs=BODY_PAIRS):
        """Report energy of the system with reduced function calls and loop overhead"""
        e = 0.0
        
        # Compute potential energy
        for body1, body2 in pairs:
            ((x1, y1, z1), v1, m1) = bodies[body1]
            ((x2, y2, z2), v2, m2) = bodies[body2]
            
            # Inline delta computation
            dx = x1 - x2
            dy = y1 - y2
            dz = z1 - z2
            
            # Inline energy computation
            e -= (m1 * m2) / ((dx * dx + dy * dy + dz * dz) ** 0.5)

        # Compute kinetic energy
        for (r, [vx, vy, vz], m) in bodies.values():
            e += m * (vx * vx + vy * vy + vz * vz) / 2.
            
        return e

    def offset_momentum(ref, bodies=BODIES):
        """Offset momentum of the system with reduced function calls"""
        px = py = pz = 0.0
        
        for (r, [vx, vy, vz], m) in bodies.values():
            px -= vx * m
            py -= vy * m
            pz -= vz * m
            
        (r, v, m) = ref
        v[0] = px / m
        v[1] = py / m
        v[2] = pz / m

    def nbody(loops, reference, iterations):
        """N-body simulation with optimized performance"""
        # Set up initial conditions
        offset_momentum(BODIES[reference])
        
        # Main simulation loop
        for _ in range(loops):
            report_energy()
            for _ in range(iterations):
                advance(0.01)
            print(report_energy())

    # Run simulation
    nbody(100, 'sun', 20000)

if __name__ == '__main__':
    main()
