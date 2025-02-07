import time

# Record the start time
start_time = time.time()

"""
    N-body simulation - Final Optimized Version
    Combines all optimizations:
    - Reduced function call overhead
    - Optimized membership testing
    - Used local variables
    - Data aggregation and loop optimization
    - Pre-computed pairs
    - Improved tuple unpacking
"""

from itertools import combinations

def initialize_bodies():
    """Initialize the bodies with their positions, velocities, and masses."""
    PI = 3.14159265358979323
    SOLAR_MASS = 4 * PI * PI
    DAYS_PER_YEAR = 365.24
    
    return {
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

def advance(bodies, pairs, dt):
    """Advance the system one timestep."""
    for (body1, body2) in pairs:
        # Local references for better performance
        ([x1, y1, z1], v1, m1) = bodies[body1]
        ([x2, y2, z2], v2, m2) = bodies[body2]
        
        # Position deltas
        dx = x1-x2
        dy = y1-y2
        dz = z1-z2
        
        # Update velocities
        mag = dt * ((dx * dx + dy * dy + dz * dz) ** (-1.5))
        b1 = m2 * mag
        b2 = m1 * mag
        
        v1[0] -= dx * b1
        v1[1] -= dy * b1
        v1[2] -= dz * b1
        v2[0] += dx * b2
        v2[1] += dy * b2
        v2[2] += dz * b2
    
    # Update positions
    for body, (r, [vx, vy, vz], m) in bodies.items():
        r[0] += dt * vx
        r[1] += dt * vy
        r[2] += dt * vz

def report_energy(bodies, pairs):
    """Calculate and return the system's energy."""
    e = 0.0
    
    # Calculate potential energy
    for (body1, body2) in pairs:
        ((x1, y1, z1), v1, m1) = bodies[body1]
        ((x2, y2, z2), v2, m2) = bodies[body2]
        
        dx = x1-x2
        dy = y1-y2
        dz = z1-z2
        
        e -= (m1 * m2) / ((dx * dx + dy * dy + dz * dz) ** 0.5)
    
    # Add kinetic energy
    for body, (r, [vx, vy, vz], m) in bodies.items():
        e += m * (vx * vx + vy * vy + vz * vz) / 2.
    
    return e

def offset_momentum(bodies, ref):
    """Offset the system's momentum relative to the referenced body."""
    px = py = pz = 0.0
    
    for body, (r, [vx, vy, vz], m) in bodies.items():
        px -= vx * m
        py -= vy * m
        pz -= vz * m
    
    # Update reference body's velocity
    (r, v, m) = bodies[ref]
    v[0] = px / m
    v[1] = py / m
    v[2] = pz / m

def nbody(loops, reference='sun', iterations=20000):
    """
    Run nbody simulation.
    
    Args:
        loops: Number of loops to run
        reference: Body at center of system (default: 'sun')
        iterations: Number of timesteps to advance (default: 20000)
    """
    # Initialize system
    bodies = initialize_bodies()
    
    # Pre-compute body pairs for optimization
    body_pairs = list(combinations(bodies.keys(), 2))
    
    # Offset momentum relative to sun
    offset_momentum(bodies, reference)
    
    # Run simulation
    for _ in range(loops):
        report_energy(bodies, body_pairs)
        for _ in range(iterations):
            advance(bodies, body_pairs, 0.01)
        print(report_energy(bodies, body_pairs))

if __name__ == '__main__':
    nbody(100)  # Run simulation with default parameters

# Record the end time
end_time = time.time()

# Calculate the elapsed time
finaltime = end_time - start_time

# Print the elapsed time
print(f"Elapsed time: {finaltime} seconds")
