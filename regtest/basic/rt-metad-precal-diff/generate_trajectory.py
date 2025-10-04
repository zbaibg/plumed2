#!/usr/bin/env python3
"""
Generate a controlled trajectory with specified autocorrelation time (OU process)
for testing PRECAL_DIFF under ADAPTIVE=DIFF in metadynamics.

The distance between atoms 1 and 2 follows an Ornstein-Uhlenbeck (OU) process
around a mean value, with known stationary variance and exponential
autocorrelation: C(t) = sigma_stat^2 * exp(-t / tau_corr).

This allows testing combinations of UPDATE_FROM relative to tau_corr and the
effect of enabling/disabling PRECAL_DIFF on the first hill sigma.
"""

import numpy as np
import argparse
from scipy import signal
from scipy.optimize import curve_fit

def compute_correlation_time(x, dt, max_lag=None):
    """
    Compute autocorrelation time using scipy and exponential fit.
    
    Uses statsmodels-style acf or scipy correlate for ACF computation,
    then fits exp(-t/tau) using scipy.optimize.curve_fit for robustness.
    
    Parameters:
      - x: time series data
      - dt: time step
      - max_lag: maximum lag to compute (default: min(len(x)//4, 200))
    
    Returns:
      - tau_fitted: fitted correlation time in same units as dt
    """
    x = np.array(x)
    n = len(x)
    if max_lag is None:
        max_lag = min(n // 4, 200)
    
    # Center the data
    x_centered = x - np.mean(x)
    
    # Compute ACF using scipy correlate (unbiased normalization)
    acf_values = signal.correlate(x_centered, x_centered, mode='full', method='fft')
    acf_values = acf_values[n-1:]  # take positive lags only
    
    # Unbiased normalization: divide by (n - lag)
    lags = np.arange(min(max_lag + 1, len(acf_values)))
    counts = n - lags
    acf_normalized = acf_values[:len(lags)] / counts
    acf_normalized = acf_normalized / acf_normalized[0]  # normalize to 1 at lag 0
    
    # Select fit region: lag > 0, acf > 0, and before first negative or very small value
    fit_mask = (lags > 0) & (acf_normalized > 0.05) & (acf_normalized < 0.99)
    if not np.any(fit_mask):
        # Fallback: use AR(1) estimate from lag-1 autocorrelation
        rho1 = acf_normalized[1] if len(acf_normalized) > 1 else 0.5
        if rho1 > 0 and rho1 < 1:
            return -dt / np.log(rho1)
        return np.nan
    
    lags_fit = lags[fit_mask]
    acf_fit = acf_normalized[fit_mask]
    times_fit = lags_fit * dt
    
    # Fit exponential: C(t) = A * exp(-t / tau)
    def exp_decay(t, A, tau):
        return A * np.exp(-t / tau)
    
    try:
        # Initial guess: A ~ 1, tau from first e-folding
        tau_guess = times_fit[len(times_fit)//2] / np.log(2) if len(times_fit) > 1 else 0.2
        popt, _ = curve_fit(exp_decay, times_fit, acf_fit, 
                           p0=[1.0, tau_guess],
                           bounds=([0.5, 0.01], [1.5, 10.0]),
                           maxfev=5000)
        tau_fitted = float(popt[1])
        return tau_fitted
    except:
        # Fallback to AR(1) if fit fails
        rho1 = acf_normalized[1] if len(acf_normalized) > 1 else 0.5
        if rho1 > 0 and rho1 < 1:
            return -dt / np.log(rho1)
        return np.nan

def generate_ou_xyz_trajectory(nsteps=1000, natoms=2, timestep=0.005,
                               d0=1.5, sigma_stat=0.15, tau_corr_ps=0.2,
                               seed=42, outfile="trajectory_ou.xyz"):
    """
    Generate XYZ trajectory where the interatomic distance follows an
    Ornstein-Uhlenbeck process with given stationary std dev and correlation time.

    OU discretization (Euler-Maruyama):
      x_{k+1} = x_k + (- (x_k - d0) / tau) * dt + sqrt(2 * sigma_stat^2 / tau) * sqrt(dt) * N(0,1)

    Parameters:
      - nsteps: number of steps
      - natoms: 2 atoms
      - timestep: ps
      - d0: mean distance (nm)
      - sigma_stat: stationary std dev (nm)
      - tau_corr_ps: correlation time (ps)
    """

    rng = np.random.default_rng(seed)
    dt = timestep
    tau = tau_corr_ps
    kappa = 1.0 / tau  # mean-reversion rate (1/ps)
    noise_coeff = np.sqrt(2.0 * sigma_stat**2 * kappa)

    x = d0
    distances = []
    for _ in range(nsteps):
        xi = rng.normal(0.0, 1.0)
        x = x + (-kappa * (x - d0)) * dt + noise_coeff * np.sqrt(dt) * xi
        distances.append(x)

    filename = outfile
    with open(filename, 'w') as f:
        for step, distance in enumerate(distances):
            f.write(f"{natoms}\n")
            # PLUMED driver expects a box line with three numbers
            f.write(f"10.0 10.0 10.0\n")
            f.write(f"AR  0.000000  0.000000  0.000000\n")
            f.write(f"AR  {distance:.6f}  0.000000  0.000000\n")

    distances = np.array(distances)
    actual_mean = float(np.mean(distances))
    actual_std = float(np.std(distances, ddof=1))

    # Measure correlation time using scipy-based method
    tau_measured = compute_correlation_time(distances, dt)

    print(f"Generated OU trajectory: {filename}")
    print(f"  Seed: {seed}")
    print(f"  Output file: {filename}")
    print(f"  Steps: {nsteps}, Total time: {nsteps * timestep:.3f} ps")
    print(f"  Time step: {timestep:.6f} ps")
    print(f"  Mean distance (target): {d0:.6f} nm, measured: {actual_mean:.6f} nm")
    print(f"  Std deviation (target): {sigma_stat:.6f} nm, measured: {actual_std:.6f} nm")
    print(f"  Correlation time (target): {tau_corr_ps:.3f} ps, measured: {tau_measured:.3f} ps")

    return filename, actual_mean, actual_std, tau_measured

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate OU trajectory for MetaD tests")
    parser.add_argument('--nsteps', type=int, default=600)
    parser.add_argument('--dt', type=float, default=0.005, help='timestep in ps')
    parser.add_argument('--tau', type=float, default=0.2, help='correlation time in ps')
    parser.add_argument('--sigma', type=float, default=0.15, help='stationary std dev in nm')
    parser.add_argument('--d0', type=float, default=1.5, help='mean distance in nm')
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--out', type=str, default='trajectory_ou.xyz')
    args = parser.parse_args()

    filename, mean_measured, sigma_measured, tau_measured = generate_ou_xyz_trajectory(
        nsteps=args.nsteps, natoms=2, timestep=args.dt, d0=args.d0,
        sigma_stat=args.sigma, tau_corr_ps=args.tau, seed=args.seed,
        outfile=args.out
    )
