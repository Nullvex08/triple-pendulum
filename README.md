# Triple Pendulum Simulation

A 2D animation of a triple pendulum system with asymmetric masses.

## Physics Background

Three rods are each attached at their midpoint to the end of the previous rod, with a mass at both ends of each rod. The equations of motion are derived from the Lagrangian and form a 3×3 linear system solved at every timestep. The triple pendulum behaves chaotically, so tiny changes in initial conditions produce dramatically different trajectories.

## Features

- Real-time 2D animation of all three rods and their end masses
- Fourth-order Runge-Kutta (RK4) integrator for accurate numerical integration
- Asymmetric rod masses: each rod has an independently configurable upper and lower mass
- Total energy plot after the animation to monitor numerical drift
- Configurable initial angles, angular velocities, and masses for each pendulum
