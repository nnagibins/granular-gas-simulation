# 2D Granular Gas Simulation

Interactive simulation of a two-dimensional granular gas with inelastic collisions, real-time visualization, and kinetic energy analysis.

## Overview

This project implements a numerical simulation of a granular gas consisting of identical particles moving inside a square container. Unlike an ideal gas, particle collisions are inelastic, meaning that a part of kinetic energy is dissipated during each collision.

The simulation allows studying:

* Kinetic energy dissipation;
* Influence of the restitution coefficient;
* Formation of particle clusters (inelastic collapse);
* Computational performance of spatial partitioning algorithms.

The project was developed as a high school research project combining classical mechanics, numerical modelling, and scientific programming in Python.

---

## Features

* 2D hard-sphere particle model;
* Adjustable coefficient of restitution (`0.50 ≤ e ≤ 1.00`);
* Real-time particle animation;
* Real-time normalized kinetic energy graph;
* External energy injection ("Shake" function);
* Saving and comparing multiple energy curves;
* Collision detection accelerated using the Cell List (Spatial Hashing) algorithm;
* Built-in self-tests for physical validation.

---

## Physical Model

Assumptions:

* All particles have identical mass and radius;
* Gravity is neglected;
* Container walls are perfectly elastic;
* Particle interactions occur only during direct contact;
* Inelastic collisions are described by the coefficient of restitution `e`.

For `e = 1.0`, the total kinetic energy remains constant.

For `e < 1.0`, the system continuously loses kinetic energy and may exhibit clustering and inelastic collapse.

---

## Algorithms

### Collision Detection

A naive implementation requires checking every particle pair:

O(N²)

To improve performance, the simulation space is divided into square cells (Cell List / Spatial Hashing approach). Each particle only checks collisions inside its own cell and neighboring cells, reducing the average complexity to approximately:

O(N)

for systems with bounded particle density.

---

## Technologies

* Python 3
* NumPy
* Matplotlib

---

## Installation

Clone the repository:

```bash
git clone https://github.com/your-username/granular-gas-simulation.git
cd granular-gas-simulation
```

Install dependencies:

```bash
pip install numpy matplotlib
```

---

## Running the Simulation

Start the graphical application:

```bash
python granular_gas_simulation.py
```

Run physical self-tests:

```bash
python granular_gas_simulation.py --self-test
```

Launch with custom parameters:

```bash
python granular_gas_simulation.py --particles 500 --radius 0.007 --restitution 0.9
```

---

## Controls

| Control            | Description                                         |
| ------------------ | --------------------------------------------------- |
| Restitution Slider | Changes the coefficient of restitution in real time |
| Steps per Frame    | Changes simulation speed                            |
| Shake              | Adds random thermal energy to particles             |
| Save Curve         | Saves the current energy curve                      |
| Clear              | Removes saved curves                                |

---

## Research Applications

The simulation can be used for educational and research purposes, including:

* Granular material dynamics;
* Dissipative systems;
* Statistical mechanics;
* Computational physics;
* Scientific visualization;
* Algorithm optimization.

---

## Author

Developed as a high school research project on numerical modelling of granular gases and kinetic energy dissipation in dissipative particle systems.
