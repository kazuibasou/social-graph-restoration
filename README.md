# Social Graph Restoration via Random Walk Sampling

This repository provides C++ and Python codes of the method proposed in 

Social Graph Restoration via Random Walk Sampling. Kazuki Nakajima and Kazuyuki Shudo. 38th IEEE International Conference on Data Engineering (ICDE 2022), 2022. (to appear) [<a href="http://arxiv.org/abs/2111.11966">arXiv</a>]

Social graph restoration is a problem in which we attempt to generate a graph whose structural properties are as close as possible to the corresponding properties of the original graph from the small sample obtained using a crawling method. 
The generated graph enables us to estimate the local and global structural properties and predict the visual representation of the original graph.
In this paper, we proposed a method for restoring the original social graph from the small sample obtained by a random walk.

The Python code is easy to use and should help you to reproduce our method.
The C++ code is much faster than the Python code.

# Python

## Requirements
Require Python 3 series.

We have confirmed that our code works in the following environments.

- macOS 11.4 and 11.5
- Ubuntu 16.04.3 LTS

## Usage

First, clone this repository:

	git clone git@github.com:kazuibasou/social-graph-restoration.git

### Input file

We need a file, *graph*.txt, where *graph* indicates the name of the network and is arbitrary. 
The file should be placed in `social-graph-restoration/data/`.
In the file, each line contains two integers separated by a half-width space, where each integer represents the node's index.
The file format follows a common format for network data.

### Restoring the original graph

Go to `social-graph-restoration/py` and run the following command:

	python main.py <graph> <sample_size>

The two arguments are as follows.

#### `<graph>`
The name of the graph.

#### `<sample_size>`
The number of nodes to be queried (see [1] for the definition of querying a node).
The number must be not less than 1 and not more than N, where N denotes the number of nodes in the original graph.

#### Example
To restore the graph named `syn10000` using 1000 queried nodes via a simple random walk, go to `social-graph-restoration/py` and run the following command:

	python main.py syn10000 1000

### Output file
The generated graph, i.e., *graph*\_.txt, will be created in the folder `social-graph-restoration/gen_graph/`.
The format of the output file is the same as that of the input file.

## Note
- The generation process for *syn10000* graph will finish in a few minutes.
However, the scalability of our method to the graph size (i.e., the number of edges) is not high at present.
For example, our method implemented in Python took approximately 20 minutes to restore the Anybeat graph with 12645 nodes and 49132 edges using 10\% queried nodes.

# C++

## Requirements
Require gcc version 4.2.1 or later.

We have confirmed that our code works in the following environments.

- macOS 11.4, 11.5, 11.6
- Ubuntu 20.04.2 LTS

## Build
(i) Clone this repository:

	git clone git@github.com:kazuibasou/social-graph-restoration.git

(ii) Go to `social-graph-restoration/cpp/src`:

	cd social-graph-restoration/cpp/src

(iii) Run the make command:

	make

This generates the following structure of the directory.

	social-graph-restoration/
	??? cpp/
	   ??? bin/
	   ??? src/
	??? data/
	??? gen_graph/

If you find a file `main` in the folder `social-graph-restoration/cpp/bin/`, the build has been successfully completed.

## Usage

### Input file

We need a file, *graph*.txt, where *graph* indicates the name of the network and is arbitrary. 
The file should be placed in `social-graph-restoration/data/`.
In the file, each line contains two integers separated by a half-width space, where each integer represents the node's index.
The file format follows a common format for network data.

### Restoring the original graph

Go to `social-graph-restoration/cpp/bin/` and run the following command:

	./main <graph> <sample_size>

The two arguments are as follows.

#### `<graph>`
The name of the graph.

#### `<sample_size>`
The number of nodes to be queried (see [1] for the definition of querying a node).
The number must be not less than 1 and not more than N, where N denotes the number of nodes in the original graph.

#### Example
To restore the graph named `syn10000` using 1000 queried nodes via a simple random walk, go to `social-graph-restoration/cpp/bin/` and run the following command:

	./main syn10000 1000

### Output file
The generated graph, i.e., *graph*\_.txt, will be created in the folder `social-graph-restoration/gen_graph/`.
The format of the output file is the same as that of the input file.

### Note
- The generation process for *syn10000* graph will finish in a few tens of seconds.
However, the scalability of our method to the graph size (i.e., the number of edges) is not high at present.
For example, our method implemented in C++ took approximately 12 hours to restore the YouTube graph with three million edges using 1\% queried nodes (see [1] for details).
- All the simulation code to reproduce our experimental results shown in the paper is available at <a href="https://www.dropbox.com/sh/qrtxb1p7ifhd58f/AADBseKsUzVqPge2ZEvtvDKNa?dl=0">here</a>.

# Reference

[1] Social Graph Restoration via Random Walk Sampling. Kazuki Nakajima and Kazuyuki Shudo. 38th IEEE International Conference on Data Engineering (ICDE 2022), 2022. (to appear) [<a href="http://arxiv.org/abs/2111.11966">arXiv</a>]

# License

This source code is released under the MIT License, see LICENSE.txt.

# Contact
- Kazuki Nakajima (https://kazuibasou.github.io/index_en.html)
- nakajima.k.an [at] m.titech.ac.jp

(Last update: 2021/12/12)