# src/methods/numerical_methods.py

class AdvectionSolver:
    """Base class for advection equation numerical solvers."""

    def __init__(self, config):
        self.config = config

    def solve(self, U0):
        """Solve the advection equation using the specified method.

        Args:
            U0 (list): Initial conditions

        Returns:
            list: Solution at final time step
        """
        U_list = [U0]
        for _ in range(self.config.final_time_step):
            U_new = self.update_u(U_list[-1])
            U_list.append(U_new)
        return U_list[-1]


class UpwindSolver(AdvectionSolver):
    """Implementation of Upwind method."""

    def update_u(self, U_old):
        CFL = -self.config.a * self.config.dx * self.config.dt
        U_new = []

        for i in range(self.config.nodes_num):
            if i == 0:
                U_new_i = -CFL * U_old[self.config.nodes_num - 2] + (1 + CFL) * U_old[i]
            else:
                U_new_i = -CFL * U_old[i - 1] + (1 + CFL) * U_old[i]
            U_new.append(U_new_i)

        return U_new


class LaxWendroffSolver(AdvectionSolver):
    """Implementation of Lax-Wendroff method."""

    def update_u(self, U_old):
        alpha = -self.config.dt / (2 * self.config.dx)
        beta = self.config.dt ** 2 / (2 * self.config.dx ** 2)
        U_new = []

        for i in range(self.config.nodes_num):
            if i == 0:
                U_new_i = alpha * (U_old[i + 1] - U_old[self.config.nodes_num - 2]) + \
                          beta * (U_old[i + 1] + U_old[self.config.nodes_num - 2]) + \
                          (1 - 2 * beta) * U_old[i]
            elif i == self.config.nodes_num - 1:
                U_new_i = alpha * (U_old[1] - U_old[i - 1]) + \
                          beta * (U_old[1] + U_old[self.config.nodes_num - 2]) + \
                          (1 - 2 * beta) * U_old[i]
            else:
                U_new_i = alpha * (U_old[i + 1] - U_old[i - 1]) + \
                          beta * (U_old[i + 1] + U_old[i - 1]) + \
                          (1 - 2 * beta) * U_old[i]
            U_new.append(U_new_i)

        return U_new


class BeamWarmingSolver(AdvectionSolver):
    """Implementation of Beam-Warming method."""

    def update_u(self, U_old):
        alpha = -self.config.dt / (2 * self.config.dx)
        beta = self.config.dt ** 2 / (2 * self.config.dx ** 2)
        U_new = []

        for i in range(self.config.nodes_num):
            if i == 0:
                U_new_i = (3 * alpha + beta + 1) * U_old[i] + \
                          (-4 * alpha - 2 * beta) * U_old[i - 1] + \
                          (alpha + beta) * U_old[self.config.nodes_num - 3]
            else:
                U_new_i = (3 * alpha + beta + 1) * U_old[i] + \
                          (-4 * alpha - 2 * beta) * U_old[i - 1] + \
                          (alpha + beta) * U_old[i - 2]
            U_new.append(U_new_i)

        return U_new