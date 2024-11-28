# Deterministic-Network-Calculus
DNC code for the paper 'Efficient Exploration on Worst-Case Delay Performance of Networked Industrial Control Systems via Network Calculus and Deep Learning,' where I am the second author. 
Accepted by the Symposium on Reliable Distributed Systems (SRDS 2024).
# Deterministic-Network-Calculus

## Input File

* `filename`：file name
* `pred_filename`：Predicted topology filename for GNN

### Network Data's Format

The input data uses a `.csv` file where:

* The first line is the network contained (including sink) nodes

```
0,1,2,3,4,5...
```

* The second row is for all nodes in the network except the sink node:

```
0,1,2,3,4...
```

* The third line is the arrival curve constraints, a $$ \rho $$ followed by a $$ \delta $$. In order, the 0th, 1st, 2nd... Arrival curve constraints for bar flow

```
5000000,512,10000000,512,10000000,512...
```

* The fourth line is the service curve constraints, a $$R$$ followed by a $$T$$, in order to represent the service curve constraints of node 0,1,2... The fourth line is service curve constraints, a $$R$$ followed by a $$T$$ in order of node 0,1,2.

```
100000000,0.000005,100000000,0.000005,100000000,0.000005,100000000,0.000005,100000000,0.000005....

```

* A number in the fifth line indicating the amount of flow

```
3
```

* The next n lines, each with x numbers, indicate the nodes this flow passes through from the beginning to the end.

```
-1, 0, 1, 7, 8
-1, 2, 3, 7, 9
-1, 4, 5, 6, 10
```

The first node is specified as -1 for the convenience of traversal when introducing constraints, node number -1 does not need to be represented in the network, and also the last node must be a sink node.

### CTO's Format

Use `csv` file input containing n lines, each line representing a CTO

Each line is formatted as follows:

```
t1, t2, t3, t4, t5, t6, t7, t8, t9
```

## run

`python fcfs.py`

