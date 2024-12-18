# GERG-2008 Equation of State

This repository provides an implementation of the GERG-2008 equation of state in Python. It calculates thermodynamic properties of gas mixtures based on given components and conditions.

## Features

- Calculate thermodynamic properties (P, T, rho, h, s, u, etc.) from specified inputs.
- Supports a variety of components defined in the `src/components` directory.
- Uses `numpy` and `scipy` for mathematical computations and solving nonlinear equations.

## Getting Started

### Prerequisites

- Python 3.7 or higher
- `numpy`
- `scipy`

### Running the Code
Run an example:
```
python examples/main.py
```

## References
Kunz, O.; Wagner, W. The GERG-2008 Wide-Range Equation of State for Natural Gases and Other Mixtures. J. Chem. Eng. Data, 2012, 57, 3032–3091.

Zhao Q, Wang Y, Chen C. Numerical simulation of the impact of different cushion gases on underground hydrogen storage in aquifers based on an experimentally-benchmarked equation-of-state. International Journal of Hydrogen Energy. 2024 Jan 2;50:495-511.