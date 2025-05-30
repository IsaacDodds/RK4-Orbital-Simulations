# RK4-Orbital-Simulations

This repository contains a set of orbital simulations implemented using the 4th-order Runge-Kutta (RK4) method. Each script numerically solves Newtonian equations of motion under gravitational interactions to produce orbital trajectories.

---

## 📂 Directory Structure

```
RK4-Orbital-Simulations/
├── Q1_HalleysComet/
│   └── id408_Q1a.py
├── Q1_TwoPlanets/
│   └── id408_Q1b.py
├── Q2_SimplifiedModel/
│   └── id408_Q2a.py
├── Q2_CenterOfMass/
│   └── id408_Q2b.py
├── README.md
├── requirements.txt
```

---

## 🚀 Simulations

### Q1a - Halley's Comet

* Simulates Halley's comet orbiting the Sun in 2D.
* Uses RK4 integration on a two-variable system.
* Visualizes the elliptical orbit over 76 years.

### Q1b - Two Planets

* Simulates two bodies (e.g., Jupiter and Saturn) orbiting the Sun.
* Includes mutual gravitational interactions.
* Visualizes planetary orbits over time.

### Q2a - Simplified Model

* Generalization of Q1b for different planetary masses and radii.
* Includes parameterized solar system setup.
* Repeats simulations over short and long terms.

### Q2b - Center-of-Mass Frame

* Adds the Sun as a moving body.
* Computes and shifts to the center-of-mass frame.
* More accurate representation of system dynamics.

---

## 📦 Requirements

Install dependencies via pip:

```bash
pip install numpy matplotlib
```

---

## 📈 Output

Each script produces a 2D plot of the orbits in billion meters. The plots illustrate the dynamics of bodies over realistic astronomical timescales.

---

## 🛠️ Usage

Navigate to any subfolder and run:

```bash
python3 <filename>.py
```

---

## 📚 License

MIT License (add `LICENSE` file if needed)

---

## 📝 Acknowledgements

Gravitational parameters from:

* NASA Planetary Fact Sheets
* IAU 2015 Resolution B3 constants

Numerical method based on standard RK4 ODE integration.

