#!/usr/bin/env python3
"""
Example: Computing H-values for Riemann Zeros
Derivative amplitude H(ρ) = |ζ'(ρ)| computation using mpmath

Author: Gongshan Liu
License: MIT
"""

from mpmath import mp, zeta as mp_zeta

# Set precision to 50 decimal digits
mp.dps = 50


def compute_H(t, h=1e-8):
    """
    Compute derivative amplitude H(t) = |ζ'(1/2 + it)| using Richardson extrapolation.
    
    Parameters:
    -----------
    t : float
        Imaginary part of zero location
    h : float
        Step size for finite differences (default: 1e-8)
    
    Returns:
    --------
    float : H(t) = |ζ'(1/2 + it)|
    """
    s = mp.mpc(0.5, t)
    h1, h2 = mp.mpf(h), mp.mpf(h/2)
    
    # Five-point Richardson extrapolation
    f_p_h1 = mp_zeta(s + mp.mpc(0, h1))
    f_m_h1 = mp_zeta(s - mp.mpc(0, h1))
    d1 = (f_p_h1 - f_m_h1) / (2 * h1)
    
    f_p_h2 = mp_zeta(s + mp.mpc(0, h2))
    f_m_h2 = mp_zeta(s - mp.mpc(0, h2))
    d2 = (f_p_h2 - f_m_h2) / (2 * h2)
    
    derivative = (4 * d2 - d1) / 3
    return float(abs(derivative))


# Example usage
if __name__ == "__main__":
    # First few Riemann zeros
    zeros = [14.134725, 21.022040, 25.010858, 30.424876, 32.935062]
    
    print("Computing H-values for first 5 Riemann zeros:")
    print("-" * 50)
    print(f"{'t':<12} {'H(t) = |ζ\'(1/2+it)|':<20}")
    print("-" * 50)
    
    for t in zeros:
        H = compute_H(t)
        print(f"{t:<12.6f} {H:<20.6f}")
    
    print("-" * 50)
    print("\nNote: Computation may take a few seconds per zero")
    print("due to high-precision arithmetic (50 digits).")