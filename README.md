# Shock System Solver

This repository contains MATLAB scripts for solving oblique shock problems in supersonic flow.

## Files

1. `Shock_System_Solver.m`: Main script for solving a specific shock system problem.
2. `MP.m`: Function to calculate Mach number and Pressure after an oblique shock.
3. `obliquerelations.m`: Function to solve oblique shock relations for various input combinations.

## Description

This set of MATLAB scripts is designed to analyze and solve problems involving oblique shocks in supersonic flow. The system can handle various input combinations and solve for unknown parameters using numerical methods.

### Shock_System_Solver.m

The main script that sets up and solves a specific shock system problem. It uses an iterative approach to find the solution for a given set of initial conditions[1].

Key features:
- Defines initial conditions (Mach number, angles, pressure)
- Iteratively solves for shock angles and pressures
- Outputs final pressures and shock angle

### MP.m

A function that calculates the Mach number and Pressure after an oblique shock, given the initial conditions[2].

### obliquerelations.m

A comprehensive function for solving oblique shock relations. It can handle various input combinations[3]:
- Mach number and shock angle (beta) to find deflection angle (theta)
- Deflection angle (theta) and shock angle (beta) to find Mach number
- Mach number and deflection angle (theta) to find shock angle (beta)

The function uses numerical methods (Newton's method) for cases that cannot be solved directly.

## Usage

To use this shock system solver:

1. Ensure all three files are in the same MATLAB directory.
2. Open and run `Shock_System_Solver.m` in MATLAB.
3. Modify the initial conditions in `Shock_System_Solver.m` as needed for different problems.

## Notes

- The solver is currently set up for a specific problem, but can be adapted for other shock system configurations.
- The `obliquerelations.m` function is versatile and can be used independently for various oblique shock calculations.
- Ensure you have a basic understanding of supersonic flow and oblique shock theory to interpret the results correctly.

## Dependencies

These scripts require MATLAB to run. No additional toolboxes are necessary."# Shock-System-Solver" 
