import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Symbolic variables setup
r, ε, k, A = sp.symbols('r ε k A', positive=True)
π = sp.pi

# Core pressure functions
def p(r_val):
    return 1/(56 * np.pi * r_val**2)

def m(r_val, ε_val):
    return (3/(14)) * (r_val - ε_val)

# Metric function g(r)
def g(r_val, ε_val):
    numerator = m(r_val, ε_val) + 4 * np.pi * p(r_val) * r_val**3
    denominator = r_val**2 * (1 - (2 * m(r_val, ε_val))/r_val)
    return numerator / denominator

# Pressure perturbation term (dp)
def dp(r_val, ε_val, k_val):
    numerator = 49 * k_val * ε_val * np.sqrt(r_val * (r_val - 2 * m(r_val, ε_val)))
    denominator = (4*r_val + 3*ε_val)**2 * (
        1 + 4*k_val*np.pi*quad(lambda t: (
            49 * t**2 * ε_val / (
                (4*t + 3*ε_val)**2 * np.sqrt(1 - 2*m(t, ε_val)/t)
            )
        ), ε_val, r_val)[0]
    )
    return numerator / denominator

# User interaction and computation
def main():
    # Get user inputs
    ε_input = float(input("Enter ε value (e.g., 0.001): "))
    A_input = float(input("Enter A value (pressure constant): "))
    r_max = float(input("Enter maximum r value for plot: "))
    
    # Compute perturbation parameter k
    k_value = A_input - 1/(56 * np.pi * ε_input**2)
    
    # Generate r values
    r_values = np.linspace(ε_input*1.1, r_max, 100)
    
    # Compute pressure profiles
    base_pressure = [p(r) for r in r_values]
    perturbed_pressure = [p(r) + dp(r, ε_input, k_value) for r in r_values]

    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(r_values, base_pressure, label='Base Pressure p(r)')
    plt.plot(r_values, perturbed_pressure, '--', label='Perturbed Pressure p(r)+dp(r)')
    plt.xlabel('Radius (r)')
    plt.ylabel('Pressure')
    plt.title('Conformal Fluid Pressure Profile')
    plt.legend()
    plt.grid(True)
    plt.show()
    
    # Print key parameters
    print(f"\nComputed parameters:")
    print(f"Perturbation parameter k = {k_value:.4e}")
    print(f"Central pressure p(ε) = {p(ε_input):.4e}")
    print(f"Metric function g(ε) = {g(ε_input, ε_input):.4f}")

if __name__ == "__main__":
    main()
