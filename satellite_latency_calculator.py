"""
satellite_latency_sim.py
------------------------------------------------
Calculation of communication latency between a
ground station and a satellite.
------------------------------------------------
"""

import numpy as np
import time
from datetime import datetime, timezone


#  CONSTANTS

c = 3e8                         # speed of light (m/s)
Re = 6.371e6                    # Earth radius (m)
mu = 3.986004418e14             # Earth's gravitational parameter (m^3/s^2)
Earth_rotations_per_day = 1.0027379093   # Number of rotations the earth makes in a day relative to the stars

def mean_motion(a: float) -> float:
    """Mean motion n = sqrt(mu / a^3) """
    return np.sqrt(mu / a**3)

def solve_kepler(M: float, e: float, tol: float = 1e-8) -> float:
    """
    Newton–Raphson solution for Kepler’s Equation
    M = E - e*sin(E)
    """
    E = M   #initial guess
    count = 0
    while True:
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E -= dE
        count += 1
        if abs(dE) < tol or count > 10000:  #10000 iterations
            return E
            

def true_anomaly(E: float, e: float) -> float:
    """Compute true anomaly θ from eccentric anomaly E """
    cos_theta = (e - np.cos(E)) / (e*np.cos(E) - 1)
    return np.arccos(cos_theta)

def rotation_matrix(r_asc: float, i: float, ω: float) -> np.ndarray:
    """Combined 3-axis rotation R3(Ω)*R1(i)*R3(ω) """
    R3Ω = np.array([[np.cos(r_asc), -np.sin(r_asc), 0],
                    [np.sin(r_asc),  np.cos(r_asc), 0],
                    [0, 0, 1]])
    R1i = np.array([[1, 0, 0],
                    [0, np.cos(i), -np.sin(i)],
                    [0, np.sin(i),  np.cos(i)]])
    R3ω = np.array([[np.cos(ω), -np.sin(ω), 0],
                    [np.sin(ω),  np.cos(ω), 0],
                    [0, 0, 1]])
    return R3Ω @ R1i @ R3ω

def satellite_position(a, e, θ, r_asc, i, w):
    """Satellite position vector in ECI frame """
    r_mag = a * (1 - e**2) / (1 + e*np.cos(θ))
    r_ijk = np.array([r_mag*np.cos(θ), r_mag*np.sin(θ), 0])
    R = rotation_matrix(r_asc, i, w)
    return R @ r_ijk

def ground_station(phi: float, lam: float) -> np.ndarray:
    """
    Ground-station position vector
    """
    agi = 120.07  # From Astronomers Almanac, recorded on Jan 21, 2019 8:00:16.7646 UTC
    start = datetime(2019, 1, 21, tzinfo=timezone.utc)
    days_passed = (datetime.now(timezone.utc) - start).days
    print(days_passed)
    # (Eq. 4.1) Greenwich sidereal angle progression
    ag = agi + np.deg2rad(360 * (Earth_rotations_per_day* days_passed))
    # Reduce to within 0–2π
    ag = np.mod(ag, 2 * np.pi)

    # (Eq. 4.2) Local angle α = αg + λ
    a = ag + lam

    # (Eq. 4.6) Ground-station position vector
    x = Re * np.cos(phi) * np.cos(a)
    y = Re * np.cos(phi) * np.sin(a)
    z = Re * np.sin(phi)

    return np.array([x, y, z])


def latency(r_s, r_site):
    """Round-trip latency t = 2|r_s − r_site| / c """
    ρ = np.linalg.norm(r_s - r_site)
    return (2 * ρ) / c

# Main function

def jupiter3_example():
    """
    Calculates latency in communication with a satellite.
    """
    ### For changing satellite data, modify these values:

    # Orbital elements (from NORAD data in the paper)
    a  = 42164e3                    # semi-major axis (m)
    e  = 0.0001971                  # eccentricity
    i  = np.deg2rad(0.0132)         # inclination
    r_asc  = np.deg2rad(63.4038)    # Right Ascension
    w  = np.deg2rad(255.8347)       # Argument of Perigee
    M0 = np.deg2rad(40.7219)        # Mean Anomaly
    # Ground-station coordinates (Germantown, MD)
    phi = np.deg2rad(39.1732)       # Latitude
    lam = np.deg2rad(-77.2717)      # Longitude

    
    # Current Epoch Time
    delta_t = time.time()

    # Step 1: Mean motion & Mean anomaly
    n = mean_motion(a)
    M = M0 + n * delta_t

    # Step 2: Solve Kepler’s Eqn
    E = solve_kepler(M, e)  #eccentric anomaly

    # Step 3: True anomaly
    θ = true_anomaly(E, e)

    # Step 4: Satellite & ground-station vectors
    r_s = satellite_position(a, e, θ, r_asc, i, w)
    r_site = ground_station(phi, lam)

    # Step 5: Latency
    t_latency = latency(r_s, r_site)

    # Display results
    print("\n================= SATELLITE LATENCY SIMULATION =================")
    print(f"Simulation timestamp (UTC): {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Current Unix Epoch time (s): {delta_t:.2f}\n")

    print("Satellite Orbital Parameters:")
    print(f"  Semi-major axis (a)        = {a/1e3:.2f} km")
    print(f"  Eccentricity (e)           = {e}")
    print(f"  Inclination (i)            = {np.rad2deg(i):.6f}°")
    print(f"  Right Ascension (Ω)        = {np.rad2deg(r_asc):.6f}°")
    print(f"  Argument of Perigee (ω)    = {np.rad2deg(w):.6f}°")
    print(f"  Mean Anomaly (M₀)          = {np.rad2deg(M0):.6f}°")
    print("\nGround Station:")
    print(f"  Latitude (φ)               = {np.rad2deg(phi):.4f}°")
    print(f"  Longitude (λ)              = {np.rad2deg(lam):.4f}°")
    print("-----------------------------------------------------------------")

    print("Calculated Orbital Results:")
    print(f"  Mean motion (n)            = {n:.8e} rad/s")
    print(f"  Mean anomaly (M)           = {np.rad2deg(M % (2*np.pi)):.6f}°")
    print(f"  Eccentric anomaly (E)      = {np.rad2deg(E):.6f}°")
    print(f"  True anomaly (θ)           = {np.rad2deg(θ):.6f}°")
    print(f"  Satellite–Ground distance  = {np.linalg.norm(r_s - r_site)/1e3:.2f} km")
    print(f"  Round-trip latency         = {t_latency*1e3:.2f} ms")
    print("=================================================================\n")



# Run directly
if __name__ == "__main__":
    jupiter3_example()
