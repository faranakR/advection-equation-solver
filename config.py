# src/config.py
import math


class Config:
    """Configuration class for advection equation solver."""

    def __init__(self):
        self.a = 1  # advection coefficient
        self.omega = 1
        self.nodes_num = 200
        self.dx = self.omega / (self.nodes_num - 1)
        self.dt = 0.5 * self.dx
        self.t_final = 1
        self.final_time_step = int(self.t_final / self.dt)


def generate_initial_condition(config, smooth=True):
    """Generate initial conditions for the simulation.

    Args:
        config: Configuration object
        smooth (bool): If True, use smooth initial condition, else use non-smooth

    Returns:
        list: Initial conditions
    """
    if smooth:
        return [math.sin(2 * math.pi * i * config.dx) for i in range(config.nodes_num)]
    else:
        U0 = []
        n = int((1 / 3) / config.dx)
        return [1 if i <= n else 0 for i in range(config.nodes_num)]


def generate_exact_solution(config, t, smooth=True):
    """Generate exact solution for comparison.

    Args:
        config: Configuration object
        t (float): Time point
        smooth (bool): If True, use smooth solution, else use non-smooth

    Returns:
        list: Exact solution
    """
    if smooth:
        return [math.sin(2 * math.pi * (i * config.dx - t)) for i in range(config.nodes_num)]
    else:
        return [1 if 2 / 3 <= i * config.dx <= 1 else 0 for i in range(config.nodes_num)]