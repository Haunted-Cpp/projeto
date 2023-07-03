# Hypermotif
## An efficient tool for motif discovery in hypergraphs



##### Building
In order to build our application simply execute:

```sh
make
```
An executable filed called "hypermotif" will be created.

To execute follow the given structure:
```sh
./hypermotif -s <integer> -i <input_file> [... optinal arguments...]
```

##### Input Format

An edge list format is used.
A Hypergraph is encoded as a list of edges in the following form:

```
e11 e12 ... e1n         
e21 e22 ... e2m         
...
en1 en2 ... enk 
```

A valid input is for the following hypergraph:

<img src="./images/hypergraph.svg" width=500>

is:

```
3 5
1 2 3
1 4 6 7
```
 
###### Input restrictions

* Each edge should have size between 2 and 4. No duplicate vertices should be included.
* However, the code will try to detect such cases and recover from a badly formatted input file.


List of arguments valid arguments:

Mandatory:

-s <size>           [Motif size - Currently only accepts 3 or 4]
-i <input>          [Input filename]

Optional:


-d                  [Specify if the output should be detailed or not]
-o <output>         [Output filename to where save results - Default is "report.txt"]
-n <k>              [
                      Specify the number of nodes in the hypergraph (In order to allow for isolated vertices) - 
                      Should be at least as large as the number of nodes used
                      By default the minimum amount of nodes needed will be used. This means all isolated vertices will be removed
                    ]
-a <k>              [
                      Select which algorithm should be used - Currently a number between 1 and 4
                      For more information about each one of them, please read the Project Report.

                      Currently each one corresponds to:
                        ESU Baseline - 1
                        ESU Modified - 2
                        Triangle - 3
                        FaSE - 4
                    ]

<--- The following are options only used for Motif discovery ---> 

-m                  [Perform motif discovery- By default a network-census is executed in the input hypergraph]
-z                  [Display the z_score instead of the SP - By default, SP is shown]
-r <k>              [Number of random networks to sample from - Integer greater than 0]
-e <k>              [Number of random shuffles between edges - Integer greater than 0]


Notice that algorithm 3 - Triangle, cannot be used with size 4. 

Sample valid commands:

./hypermotif -s 4 -i datasets/hs.edges -r 10 -m -d -a 4
./hypermotif -s 3 -i datasets/geology.edges -r 10 -m -d -a 1
./hypermotif -s 3 -i datasets/geology.edges
./hypermotif -s 4 -i datasets/geology.edges -r 10 -m -d -a 4
./hypermotif -s 3 -i datasets/geology.edges -a 1 -o geology/baseline && ./hypermotif -s 3 -i datasets/geology.edges -a 2 -o geology/mod && ./hypermotif -s 3 -i datasets/geology.edges -a 3 -o geology/triangle && ./hypermotif -s 3 -i datasets/geology.edges -a 4 -o geology/fase
./hypermotif -s 3 -i datasets/geology.edges -a 1 -o hs1/baseline -m -r 3 -e 3
./hypermotif -s 4 -i datasets/hs.edges -a 1 -o k4/m11.out -r 2 -e 2 -m && ./hypermotif -s 4 -i datasets/hs.edges -a 2 -o k4/m22.out -r 2 -e 2 -m && ./hypermotif -s 4 -i datasets/hs.edges -a 3 -o k4/m33.out -r 2 -e 2 -m&& ./hypermotif -s 4 -i datasets/hs.edges -a 4 -o k4/m44.out -r 2 -e 2 -m


A list of sample datasets is provided in the folder "datasets".

Random Dataset:
https://drive.google.com/file/d/1T3cov10rd3RAxpXrbYXrQqFXy5eBhR8D/view?usp=sharing


All the steps have been tested on:
- Ubuntu 20.04.5 LTS
- g++ (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0


###### License

A 

However, this software uses nauty by Brendan McKay. Therefore, nauty's license restrictions apply to the usage of hypermotif.