# LR-bound-for-approximation

This repository contains the code used in https://arxiv.org/abs/2311.12732 in order to compute a lower bound on the achievable approximation ratio of MaxCut over cubic graphs using quantum annealing.


Installation
=============

Clone the repository:

```bash
git clone https://github.com/Arts-Braido/LR-bound-for-approximation.git
```

Move to the cloned repo:


```bash
cd LR-bound-for-approximation
```

Create a virtual env:

```bash
virtualenv venv
```

Activate it:

```bash
source venv/bin/activate
```

Install the repo:

```bash
pip install .
```


Enumeration and simulation
==========================

A first script `scripts/enum_ball.py` will enumerate all balls of radius up to 3 and store them on disk.
The enumeration process is quite time-consuming (up to 20h on a laptop).

A second script `scripts/balls_energy.py` will load each ball, simulate the corresponding quantum annealing and compare the result against the analytical LR bound and the target approximation ratio (see the paper for more details). Store printouts in a text file named energy_log.txt. If any ball induces a ratio lower than the target ratio, the script will crash with an assertion error. This should never happen. This script tooks several weeks to complete on a HPC hardware. 


A jupyter notebook processes the information stored in the energy_log.txt file and retrieves the values needed to calculate the approximation ratio. It is also possible to obtain the energy plots of the central edge X of the "worst" balls against alpha.