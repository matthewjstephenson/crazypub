import numpy as np
from scipy.linalg import solve
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def spectral_setup(Ncut, Lmax):
    """Create Chebyshev-Gauss-Lobatto grid and differentiation matrix"""
    # Chebyshev grid points
    j = np.arange(Ncut + 1)
    z = Lmax/2 * (np.cos(np.pi * j/Ncut) + 1)
    
    # Differentiation matrix
    D = np.zeros((Ncut+1, Ncut+1))
    for i in range(Ncut+1):
        for j in range(Ncut+1):
            if i != j:
                c_i = 1 if i == 0 or i == Ncut else 2
                c_j = 1 if j == 0 or j == Ncut else 2
                D[i,j] = (c_i/c_j) * (-1)**(i+j) / (np.sin((i+j)*np.pi/(2*Ncut)) * np.sin((i-j)*np.pi/(2*Ncut)))
    
    D -= np.diag(np.sum(D, axis=1))  # Diagonal entries
    return z, D

def solve_harmonic_oscillator(Ncut, Lmax, x0, v0, force_func=None, t_span=(0, 2)):
    """Spectral method solver for harmonic oscillator"""
    z, D = spectral_setup(Ncut, Lmax)
    N = len(z)
    
    # Create system matrix (y' = A@y + b)
    A = np.block([[D[1:,1:], np.eye(N-1)],
                 [-np.eye(N-1), D[1:,1:]]])
    
    # Create boundary condition vector
    b = np.concatenate([v0 * D[1:,0], x0 * D[1:,0]])
    
    if force_func:
        # Add forcing term
        force = force_func(z[1:])
        b -= np.concatenate([force, np.zeros(N-1)])
    
    # Solve linear system
    y = solve(A, -b)
    
    # Extract position solution
    x_sol = np.concatenate([[x0], y[N-1:]])
    return z, x_sol

def analytical_solution(t, x0, v0):
    """Analytical solution for harmonic oscillator"""
    return x0 * np.cos(t) + v0 * np.sin(t)

# User input section
if __name__ == "__main__":
    # Get parameters from user
    Ncut = int(input("Enter number of grid points (recommend 16-32): "))
    Lmax = float(input("Enter domain size (e.g., 1.5): "))
    x0 = float(input("Enter initial position: "))
    v0 = float(input("Enter initial velocity: "))
    use_force = input("Add forcing function? (y/n): ").lower() == 'y'
    
    # Define forcing function if needed
    if use_force:
        force_func = lambda t: np.sin(10*t)
    else:
        force_func = None
    
    # Spectral method solution
    z, x_spectral = solve_harmonic_oscillator(Ncut, Lmax, x0, v0, force_func)
    
    # Numerical reference solution
    def oscillator(t, y):
        x, v = y
        dxdt = v
        dvdt = -x + (force_func(t) if force_func else 0)
        return [dxdt, dvdt]
    
    sol = solve_ivp(oscillator, [0, Lmax], [x0, v0], t_eval=z, rtol=1e-8)
    
    # Plot results
    plt.figure(figsize=(10, 6))
    plt.plot(z, x_spectral, 'bo-', label='Spectral Method')
    plt.plot(sol.t, sol.y[0], 'r--', label='Numerical Solution')
    plt.xlabel('Time')
    plt.ylabel('Position')
    plt.title('Harmonic Oscillator Solution Comparison')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Print key metrics
    print("\nSolution Summary:")
    print(f"{'Time':<10}{'Spectral':<15}{'Numerical':<15}")
    for t, xs, xn in zip(z, x_spectral, sol.y[0]):
        print(f"{t:<10.3f}{xs:<15.6f}{xn:<15.6f}")
