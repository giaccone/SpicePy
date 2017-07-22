# SpicePy
**SpicePy** is a name coming from the merge of **SPICE** (*Simulation Program with Integrated Circuit Emphasis*) and **Python**, hence, it goes without saying that it is a _Circuit simulator written in python_

**SpicePy** borns as a teaching project. It is shared with students of *basic circuit theory* with two aims:

* to allow them to check the results of exercises solved analytically
* to show them how a numerical code to solve circuit is made

# Example
Let us consider the network described by the file *neltlist.net* including the following lines:
```
V1 3 0 5
V2 2 0 15
R1 1 3 10
R2 1 2 5
R3 1 0 2
I1 2 1 2
```

one can solve it executing the following script:
```python
import netlist as ntl
import scipy.sparse.linalg as ssl

net = ntl.Network('network.net')
net.dc_solve()
net.branch_voltage()
```

On the terminal/shell the following result is obtained:
```
=====================
   branch voltages   
=====================
v(V1) = 5.0000 V
v(V2) = 15.0000 V
v(R1) = 1.8750 V
v(R2) = -8.1250 V
v(R3) = 6.8750 V
v(I1) = 8.1250 V
```
