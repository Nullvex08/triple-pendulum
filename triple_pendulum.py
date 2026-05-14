#!/usr/bin/env python

# Name: energy_triple.py
# Description: This Code visualizes a triple pendulum, 
# where all rods are attached at its middle point, and a
# Version: 1.0.0

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation

# Const
g = 9.81        # m/s^2
dt = 0.0008      # s
L = 1           # length of rods in m
STEPS_PER_FRAME = 60
ANIMATION_AXIS_LIMIT = L * 3.5

# Variables
energy_history = []
time_history = []

class pendulum:
    length = L
    mass_r = 0
    all_pendulums = []
    def __init__(self, mass_up, mass_down, init_angle, init_angle_velocity):
        self.mass_d = mass_down                                        # kg
        self.mass_u = mass_up                                          # kg
        self.init_angle = np.deg2rad(init_angle)                       # theta_ pendulum.all_pendulums.index(self)  in rad
        self.init_angle_velocity = init_angle_velocity/360*2*np.pi     # rad/s
        pendulum.all_pendulums.append(self)

def energy(state):
    theta_1, omega_1, theta_2, omega_2, theta_3, omega_3 = state
    l = pendulum.length
    rod = pendulum.mass_r
    m11 = P_1.mass_d
    m12 = P_1.mass_u
    m21 = P_2.mass_d
    m22 = P_2.mass_u
    m31 = P_3.mass_d
    m32 = P_3.mass_u
    tot_energy = (1/2 * (m11+m12+m21+m22+m31+m32) *
                  omega_1**2 * l**2 + 1/2 * (m21+m22+m31+m32)
                  * omega_2**2 * l**2 + (m21+m31+m32-m22) * l**2
                  * omega_1 * omega_2 * np.cos(theta_1 - theta_2) +1/2 * (m31+m32) * l**2 *
                  omega_3 **2 + (m31-m32) * l**2 * (omega_1 * omega_3 * np.cos(theta_1 - theta_3) + omega_2 * omega_3 * np.cos(theta_2 - theta_3))
                  + 1/2 * rod * (2* omega_1**2 * l**2 + omega_2**2 * l**2 + 2 * l**2 * omega_1 * omega_2 * np.cos(theta_1 - theta_2)) + 1/6 * (omega_1**2 + omega_2**2 + omega_3**2) * l**2 * rod + (g * l * np.cos(theta_1) * (m12-m11-m21-m22-m31-m32) + g * l * np.cos(theta_2) * (m22 - m21-m31-m32) + g * l * np.cos(theta_3) * (m32-m31) - g * rod * l * (2 * np.cos(theta_1) + np.cos(theta_2))))
    return tot_energy

def calc_derivative(state, t):
    theta_1, omega_1, theta_2, omega_2, theta_3, omega_3 = state
    dtheta_1 = omega_1
    dtheta_2 = omega_2
    dtheta_3 = omega_3
    l = pendulum.length
    rod = pendulum.mass_r
    m11 = P_1.mass_d
    m12 = P_1.mass_u
    m21 = P_2.mass_d
    m22 = P_2.mass_u
    m31 = P_3.mass_d
    m32 = P_3.mass_u
    M = np.array([
        [l**2 * (m11+m12+m21+m22+m31+m32 + 7/3*rod), l**2 * np.cos(theta_1-theta_2) * (m21+m31+m32-m22 + rod) , l**2 * np.cos(theta_1-theta_3) *(m31-m32)],
        [l**2 * np.cos(theta_1 - theta_2) * (m21+m31+m32-m22+rod) , l**2 * (m21+m22+m31+m32+4/3*rod) , l**2 * np.cos(theta_2-theta_3) * (m31-m32)],
        [l**2 * np.cos(theta_1 - theta_3) * (m31-m32) , l**2 * np.cos(theta_2 - theta_3) * (m31-m32) , l**2 * (m31+m32+1/3 * rod)]
    ])
    answer = np.array([
        ((m21+m31+m32-m22) * l**2 * omega_2 * np.sin(theta_1 - theta_2) * (omega_1 - omega_2) + (m31-m32) * l**2 * omega_3 * np.sin(theta_1 - theta_3) * (omega_1 - omega_3) + l**2 * rod * omega_2 * np.sin(theta_1 - theta_2) * (omega_1 - omega_2) - (m21+m31+m32-m22-rod) * l**2 * omega_1 * omega_2 * np.sin(theta_1 - theta_2) - (m31-m32) * l**2 * omega_1 * omega_3 * np.sin(theta_1 - theta_3) - g * l * np.sin(theta_1) * (m11+m21+m22+m31+m32-m12) - np.sin(theta_1) * 2 * g * rod * l),
        ((m21+m31+m32-m22) * l**2 * omega_1 * np.sin(theta_1 - theta_2) * (omega_1 - omega_2) + (m31-m32) * l**2 * omega_3 * np.sin(theta_2 - theta_3) * (omega_2 - omega_3) + l**2 * rod * omega_1 * np.sin(theta_1 - theta_2) * (omega_1 - omega_2) + (m21+m31+m32-m22+rod) * l**2 * omega_1 * omega_2 * np.sin(theta_1 - theta_2) - (m31-m32) * l**2 * omega_2 * omega_3 * np.sin(theta_2 - theta_3) - g * l * np.sin(theta_2) * (m21+m31+m32-m22) - g * rod * l * np.sin(theta_2)),
        ((m31-m32) * l**2 * omega_1 * np.sin(theta_1 - theta_3) * (omega_1 - omega_3) + (m31-m32) * l**2 * omega_2 * np.sin(theta_2 - theta_3) * (omega_2 - omega_3) + (m31-m32) * l**2 * (omega_1 * omega_3 * np.sin(theta_1 - theta_3) + omega_2 * omega_3 * np.sin(theta_2 - theta_3)) - g * l * np.sin(theta_3) * (m31-m32))
    ])

    domega_1, domega_2, domega_3 = np.linalg.solve(M, answer)
    derivative_state = np.array([dtheta_1, domega_1, dtheta_2, domega_2, dtheta_3, domega_3])
    return derivative_state

def rk4_step(y, t, h):
    k1 = calc_derivative(y, t)
    k2 = calc_derivative(y + 0.5*h*k1, t + 0.5*h)
    k3 = calc_derivative(y + 0.5*h*k2, t + 0.5*h)
    k4 = calc_derivative(y + h*k3, t + h)
    return y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

def rod_endpoints(cx, cy, theta, rod_length):
    dx = rod_length * np.sin(theta)
    dy = -rod_length * np.cos(theta)
    xA, yA = cx - dx, cy - dy
    xB, yB = cx + dx, cy + dy
    return xA, yA, xB, yB

if __name__ == "__main__":
    t = 0.0
    fig, ax = plt.subplots()
    ax.plot(0, 0, 'ko', markersize=8)
    ax.set_xlim(-ANIMATION_AXIS_LIMIT, ANIMATION_AXIS_LIMIT)
    ax.set_ylim(-ANIMATION_AXIS_LIMIT, ANIMATION_AXIS_LIMIT)
    ax.set_aspect('equal')

    rod1, = ax.plot([], [], '-', lw=2)
    rod2, = ax.plot([], [], '-', lw=2)
    rod3, = ax.plot([], [], '-', lw=2)

    mass_points, = ax.plot([], [], 'o', markersize=8)
    P_1 = pendulum(0.5, 0.4, 180, 0)
    P_2 = pendulum(0.5, 0.2, -90, 0)
    P_3 = pendulum(0.5, 0.2, -150, 0)
    if len(pendulum.all_pendulums) !=3:
        raise ValueError("It is not a triple pendulum")

    state_start = np.array([P_1.init_angle, P_1.init_angle_velocity, P_2.init_angle, P_2.init_angle_velocity, P_3.init_angle, P_3.init_angle_velocity])
    state = state_start.copy()

    def update(frame):
        global state, t
        for _ in range(STEPS_PER_FRAME):
            state = rk4_step(state, t, dt)
            t += dt
            energy_history.append(energy(state))
            time_history.append(t)
        theta1, omega1, theta2, omega2, theta3, omega3 = state
        # Pivot points (chained)
        x0, y0 = 0.0, 0.0                          # fixed pivot of rod 1
        x1 = x0 + L * np.sin(theta1)               # pivot of rod 2
        y1 = y0 - L * np.cos(theta1)
        x2 = x1 + L * np.sin(theta2)               # pivot of rod 3
        y2 = y1 - L * np.cos(theta2)

        # Rod endpoints (each rod extends L/2 in both directions from pivot)
        xA1, yA1, xB1, yB1 = rod_endpoints(x0, y0, theta1, L)
        xA2, yA2, xB2, yB2 = rod_endpoints(x1, y1, theta2, L)
        xA3, yA3, xB3, yB3 = rod_endpoints(x2, y2, theta3, L)
        rod1.set_data([xA1, xB1], [yA1, yB1])
        rod2.set_data([xA2, xB2], [yA2, yB2])
        rod3.set_data([xA3, xB3], [yA3, yB3])

        # Masses at both ends of each rod + pivot points
        mass_points.set_data(
            [xA1, xB1, xA2, xB2, xA3, xB3],
            [yA1, yB1, yA2, yB2, yA3, yB3]
        )
        return rod1, rod2, rod3, mass_points

    ani = FuncAnimation(fig, update, interval=20, blit=True)

    plt.show()
    plt.figure()
    plt.plot(time_history, energy_history)
    plt.title("Total Energy over Time")
    plt.xlabel("Time (s)")
    plt.ylabel("Energy (J)")
    plt.grid(True)
    plt.show()
