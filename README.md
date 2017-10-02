# SpicePy
**SpicePy** is a name coming from the merge of **SPICE** (*Simulation Program with Integrated Circuit Emphasis*) and **Python**, hence, it goes without saying that it is a _Circuit simulator written in python_

**SpicePy** borns as a teaching project. It is shared with students of *basic circuit theory* with two aims:

* to allow them to check the results of exercises solved analytically
* to show them how a numerical code to solve circuit is made

This document provides information about features and installation procedure. For the user's guide please refer to the [Wiki section](https://github.com/giaccone/SpicePy/wiki).

# Installation

1. Clone the repository:
`git clone https://github.com/giaccone/SpicePy.git`

2. add SpicePy folder to the python path
3. done!

## Requirements
At the time of writing (October 4, 2017), the project is based on:

* `Python 3.6.2`
* `numpy 1.13.1`
* `scipy 0.19.1`
* `matplotlib 2.0.2`

This is not a frozen configuration. I'm use to update these packages as soon as an update is available.

Finally, SpicePy should work also with earlier version of Python (but **not** Python 2), numpy, scipy and matplotlib.

# Feature (and changelog)

## October 2, 2017

* new methods to compute/print/plot generic voltages and current
* added examples in the folder 'demo'.

(Check commit history for more details).

## September 1, 2017

* solution of dynamic networks. So far, tested first and second order circuits.

(Check commit history for more details).

## August 29, 2017

* solution of DC network
* computation of the operating point (including dynamic components)
* solution of AC network (only single frequency right now)




