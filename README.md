# Deterministic-Network-Calculus
The code associated with the paper "Efficient Exploration on Worst-Case Delay Performance of Networked Industrial Control Systems via Network Calculus and Deep Learning."

# Deterministic-Network-Calculus

## 输入数据

* `filename`：输入网络数据的路径
* `pred_filename`：输入GNN预测CTO的路径

### 输入网络格式

输入数据使用`.csv`文件，其中：

* 第一行是网络包含的（包含sink）节点

```
0,1,2,3,4,5...
```

* 第二行是网络中除了sink节点以外的所有节点：

```
0,1,2,3,4...
```

* 第三行是arrival curve constraints，一个$$ \rho $$跟着一个$$\delta$$。按顺序表示第0,1,2...条流的arrival curve constraints

```
5000000,512,10000000,512,10000000,512...
```

* 第四行是service curve constraints，一个$$R$$跟着一个$$T$$，按顺序表示第0,1,2...号节点的service curve constraints

```
100000000,0.000005,100000000,0.000005,100000000,0.000005,100000000,0.000005,100000000,0.000005....

```

* 第五行一个数字，表示flow的数量

```
3
```

* 接下来n行，每行x个数字，从开始到结束表示这一flow经过的节点。

```
-1, 0, 1, 7, 8
-1, 2, 3, 7, 9
-1, 4, 5, 6, 10
```

第一个节点规定为-1，是为了引入约束时遍历的方便，-1号节点不需要在网络中表示出来，同时，最后一个节点一定是sink节点。

### CTO格式

使用`csv`文件输入，包含n行，每一行表示一个CTO

每一行的格式如下：

```
t1, t2, t3, t4, t5, t6, t7, t8, t9
```

## 运行

`python fcfs.py`

