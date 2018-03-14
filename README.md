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
At the time of writing (January 6, 2018), the project is based on:

* `Python 3.6.4`
* `numpy 1.14.2`
* `scipy 1.0.0`
* `matplotlib 2.2.0`

This is not a frozen configuration. I'm use to update these packages as soon as an update is available.

Finally, SpicePy should work also with earlier version of Python (but **not** Python 2), numpy, scipy and matplotlib.

Last but not least, as many other Python tools a convenient method to handle SpicePy is `iPython`, a powerful interactive Python shell where Python can be used interactively.

# Work in progress

Here you can read about aspects that I'm thinking to include in this project. The list does not certifies that I will integrate all it is described. It's a simple list of the topics that potentially will be developed.

* **Controlled sources** In a local branch I've already added VCVSs. I will push this work on GitHub only when all controlled sources will be ready an tested.
* **non-linearity** I'm planning to implement a non-linear solver based on Newthon-Rapshon method. A rough implementation already exists but it is far from be harmonized to the entire project. Since I prefer *stability* over *new features* I will include non-linearity when I will have some more (free) time that currently I do not have.
* **GUI** one of my students is implementing a very basic GUI for this code. Basically, it is text-editor for netlists that includes also the tools for the simulation. Right now it is not ready for the final use because it is not stable.

# Feature (and changelog)

## March 14, 2018

* support for unit-prefix in netlists: (`'f = 1e-15'`, `'p = 1e-12'`, `'n  = 1e-9'`, `'u = 1e-6'`, `'m = 1e-3'`, `'k = 1e3'`, `'meg = 1e6'`, `'g = 1e9'`, `'t = 1e-12'`) 

## January 6, 2018

* support for ac-multi-frequency analysis
* examples of ac-multi-frequency analysis

(Check commit history for more details).

## November 16, 2017

* support for transient-sources

(Check commit history for more details).

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




