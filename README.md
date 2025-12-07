# TurtleBot Hybrid Obstacle Avoidance

A MATLAB implementation of Lyapunov-based hybrid obstacle avoidance for TurtleBot robots, based on the theoretical framework from Sanfelice et al. (ACC 2006).

## Overview

This project implements a hybrid controller that enables two TurtleBot robots to navigate to a common goal while avoiding a circular obstacle. The robots take different paths (homotopy classes) around the obstacle using mode-dependent Lyapunov functions with logarithmic barrier potentials.

### Key Features

- **Hybrid Control**: Flow/jump logic with hysteresis for robust mode switching
- **Barrier Functions**: Logarithmic barrier potentials prevent collision with obstacles
- **Dual-Path Navigation**: Two robots navigate via different homotopy classes (above/below obstacle)
- **Real-Time Control**: Integration with Vicon motion capture for real TurtleBot experiments
- **Simulation Mode**: Virtual simulation for testing and parameter tuning

## Theoretical Background

The controller implements the synergistic Lyapunov functions from Sanfelice et al.  (ACC 2006):

```
V_q(x) = 0.5 * ||x - a||² + B(d_q(x))
```

Where:
- `a` is the target position
- `d_q(x)` is the signed distance to the mode-dependent free set boundary
- `B(·)` is a logarithmic barrier function that penalizes proximity to obstacles

### Barrier Function

The barrier function (Eq. 11 from the paper):
```
B(d) = (d-1)² * log(1/d),  for 0 < d < 1
B(d) = 0,                   for d ≥ 1
B(d) = +∞,                  for d ≤ 0
```

## Files

| File | Description |
|------|-------------|
| `simulate_hybrid_two_turtlebots.m` | Virtual simulation of two robots with the hybrid controller |
| `client_hybrid_obstacle_avoidance.m` | Real-time control client for physical TurtleBots via Vicon |
| `hybrid_controller_step.m` | Core hybrid controller implementing flow/jump logic |
| `hybrid_potentials.m` | Mode-dependent Lyapunov functions and gradients |
| `barrier_B.m` | Logarithmic barrier potential function |
| `barrier_derivative_numeric.m` | Analytical derivative of the barrier function |

## Usage

### Simulation

Run the virtual simulation to visualize the hybrid controller behavior:

```matlab
simulate_hybrid_two_turtlebots
```

This will:
1. Initialize two virtual robots at positions `(0, 0. 5)` and `(0, -0.5)`
2.  Navigate both to the goal at `(3, 0)` while avoiding the obstacle at `(1, 0)`
3. Display trajectory plots, mode evolution, and convergence metrics

### Real Hardware

For physical TurtleBot experiments with Vicon motion capture:

```matlab
client_hybrid_obstacle_avoidance
```

**Prerequisites:**
- Vicon motion capture system configured
- Two TurtleBots tracked as `Object1` and `Object2`
- Network connectivity to robot ports (50804, 50805)

## Parameters

### Controller Parameters (Example 5.2 from paper)

| Parameter | Value | Description |
|-----------|-------|-------------|
| Target | `(3, 0)` | Goal position [m] |
| Obstacle center | `(1, 0)` | Obstacle position [m] |
| Obstacle radius | `1/(20√2) ≈ 0.0354` | Obstacle size [m] |
| Wedge angle | `π/4` (45°) | Mode region boundary angle |

### Hysteresis Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `μ` | 1.2 | Flow set multiplier (>1) |
| `λ` | 0.05 | Jump set margin for noise robustness |
| `γ` | 0.01 | Relaxation near goal |

## How It Works

1. **Mode Selection**: Each robot is assigned an initial mode (`q=1` or `q=2`) determining its path around the obstacle
2.  **Potential Evaluation**: At each timestep, compute `V_q(y)` using the measured position
3. **Flow/Jump Logic**: 
   - Flow (continue) if `V_q ≤ μ·min(V)` or `V_q ≤ γ`
   - Jump (switch mode) if `V_q ≥ (μ-λ)·min(V)` and another mode has lower potential
4. **Control**: Apply `u = -∇V_q(y)` to drive toward the goal while avoiding obstacles
5. **Unicycle Mapping**: Convert Cartesian velocity to wheel commands via heading alignment

## Results

The simulation produces:
- **Trajectory plot**: Shows both robots navigating around the obstacle
- **Mode evolution**: Displays discrete mode switches over time
- **Potential plots**: Visualizes `V_1` and `V_2` for each robot
- **Convergence plot**: Distance to goal over time

## References

- R. G. Sanfelice, M.  J. Messina, S. E. Tuna, and A. R. Teel, "Robust hybrid controllers for continuous-time systems with applications to obstacle avoidance and regulation to disconnected set of points," *American Control Conference (ACC)*, 2006.

## License

This project is provided for educational and research purposes. 