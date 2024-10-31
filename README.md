# Advection Equation Solver
Implementation of the advection equation solver using Upwind, Lax-Wendroff, and Beam-Warming methods, with comparisons to exact solutions.

## Methods Implemented

Upwind method (First-order HJ-ENO)

Lax-Wendroff scheme

Beam-Warming method

## Key Features

Analytical solutions using method of characteristics
Support for both smooth and non-smooth initial conditions
Comparison plots between numerical and exact solutions
Configurable parameters (nodes, time steps, etc.)
 
## Usage 

1. Configure parameters in config.py
2. Run the simulation: python src/main.py

## Results

The solver compares different numerical methods with the exact solution, showing:

1. Phase lag between numerical and exact solutions
2. Behavior at singularities
3. Oscillation characteristics in non-locally discontinuous regions