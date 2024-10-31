# src/main.py
import matplotlib.pyplot as plt
from config import Config, generate_initial_condition, generate_exact_solution
from methods.numerical_methods import UpwindSolver, LaxWendroffSolver, BeamWarmingSolver


def run_simulation(smooth=True):
    """Run the complete simulation with all methods."""
    # Initialize configuration
    config = Config()

    # Generate initial conditions and exact solution
    U0 = generate_initial_condition(config, smooth)
    U_exact = generate_exact_solution(config, config.t_final, smooth)

    # Initialize solvers
    upwind = UpwindSolver(config)
    lax = LaxWendroffSolver(config)
    beam = BeamWarmingSolver(config)

    # Solve using each method
    U_upwind = upwind.solve(U0)
    U_lax = lax.solve(U0)
    U_beam = beam.solve(U0)

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(range(config.nodes_num), U0, label="Initial function, t = 0")
    plt.plot(range(config.nodes_num), U_exact, label=f"Exact solution at t = {config.t_final}")
    plt.plot(range(config.nodes_num), U_upwind, label=f"Upwind scheme, t = {config.t_final}")
    plt.plot(range(config.nodes_num), U_lax, label=f"Lax-Wendroff, t = {config.t_final}")
    plt.plot(range(config.nodes_num), U_beam, label=f"Beam-Warming, t = {config.t_final}")

    plt.title("Advection Equation Solutions" + (" (Smooth)" if smooth else " (Non-smooth)"))
    plt.xlabel("Node number")
    plt.ylabel("u(x,t)")
    plt.legend()
    plt.grid(True)
    plt.savefig(f'results/solution_{"smooth" if smooth else "nonsmooth"}.png')
    plt.show()


if __name__ == "__main__":
    # Run both smooth and non-smooth cases
    run_simulation(smooth=True)
    run_simulation(smooth=False)