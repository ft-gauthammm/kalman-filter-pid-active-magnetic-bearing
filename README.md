# Kalman Filter-Assisted PID Control for an Active Magnetic Bearing

## Overview

This project investigates PID control and Kalman-filter-assisted
PID control for a simulated single-axis Active Magnetic Bearing (AMB).

The objective is to compare controller behaviour under measurement
noise, plant-parameter variations, and rotating-unbalance disturbances.

## Mathematical Model

The linearised plant is represented by:

m * x_ddot = ks * x + ki * i + Fd

where:
- m = 46 kg
- ki = 129.7 N/A
- ks = 2.4 × 10^6 N/m
- Reference displacement = 100 µm
- Controller sampling frequency = 5 kHz

The model represents an open-loop unstable system.

## Controllers

- Baseline PID
- Optimised PID
- Two-state Kalman-filter-assisted PID
- Three-state Kalman-filter-assisted PID

## Implementation

The project uses MATLAB for plant modelling, controller
optimization, simulation, and performance evaluation.

The main script also exports CSV tables and a MAT-file
containing simulation results.

## Repository Structure

- `matlab/` — MATLAB source files
- `results/` — CSV results and MATLAB data
- `figures/` — system diagrams
- `docs/` — mathematical model and verification notes

## Limitations

This is a single-axis linearised simulation, not a physical
AMB implementation. Amplifier dynamics, cross-axis coupling,
and hardware validation are outside the current scope.

## Execution Status

Document the MATLAB version, execution date, and verification
status after running the scripts. Do not describe results as
reproduced until the simulation has been executed and checked.
