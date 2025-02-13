# Implementation of Learning Automata-enabled Ant Colony Optimization for Group Influence Maximization in Social Networks: A Two-Stage Approach"

## Overview

This project implements an Ant Colony Optimization algorithm with an additional layer of adaptive parameter setting using Learning Automata for selecting optimal seed sets in a graph. It is particularly designed for maximizing the group influence spread in social network graphs.

## Requirements

- Python 3.x
- Required Libraries:
	- networkx
	- matplotlib
	- numpy
	- csv
	- random
	- math

## Input Format

- dataset.csv: The dataset containing the edges of the graph. Each row represents an edge with the format:

```
node1,node2,weight
```
- group.csv: The group information with each row representing a group of nodes. Each column corresponds to a node in the group.

```
node1,node2,...,nodeN
```

### Parameters

- `k`: Length of the seed set.
- `num_ants`: Number of ants in the colony.
- `num_iterations`: Number of iterations to run the optimization.
- `evaporation_rate`: Rate at which pheromones evaporate.
- `alpha, beta`: ACO control parameters.
- `delta`: Step value used in LA operations.
- `div`: Diversity parameter
- `theta`: Activation threshold