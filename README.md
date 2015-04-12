<script type="text/javascript" src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS_HTML"></script>

# Table of contents

[TOC]

# Introduction

ACE is a software package designed to quickly and accurately infer [Ising](https://en.wikipedia.org/wiki/Ising_model) or [Potts](https://en.wikipedia.org/wiki/Potts_model) models based on correlation data from a variety of biological and artificial systems. This software makes use of the <b>A</b>daptive <b>C</b>luster <b>E</b>xpansion (ACE) algorithm.

Given a set of correlation data or sequence input in [FASTA](http://en.wikipedia.org/wiki/FASTA_format) format, ACE will produce a Ising or Potts model that reproduces the input correlations to within the expected error due to finite sampling.


# Installation

Download and unzip the package, then run the following commands in the terminal from the new directory:

```bash
$ ./configure
$ make
```

If you'd like to be able to run the program from any directory, you can then enter:

```bash
$ make install
```


# Required Input

Running the algorithm requires a set of correlations as input, to be computed from your data.

As an example, let's consider a system of $N$ variables described by the configuration $\underline{x}=\{x_1, x_2,\ldots,x_N\}$, with each variable $x_i$ taking one of $q_i$ possible values, $x_i\in\{1,2,\ldots,q_i\}$. From a set of $B$ observations of the system, we can compute the frequency of each variable as well as the pairwise correlations,

$$
\begin{aligned}
\begin{align}
p_i(a) &= \frac{1}{B}\sum_{k=1}^{B}\delta(x_i,a)\,,\\
p_{ij}(a,b) &= \frac{1}{B}\sum_{k=1}^{B}\delta(x_i,a)\delta(x_j,b)\,.
\end{align}
\end{aligned}
$$

Here $\delta$ represents the [Kronecker delta function](http://en.wikipedia.org/wiki/Kronecker_delta). These correlations should be saved in a file ending with the extension `.p`, in the following format:

> $p_1(1)$ $p_1(2)$ ... $p_1(q_1-1)$ 
> $p_2(1)$ $p_2(2)$ ... $p_2(q_2-1)$
> ...
> $p_N(1)$ $p_N(2)$ ... $p_N(q_N-1)$
> $p_{1,2}(1,1)$ $p_{1,2}(1,2)$ ... $p_{1,2}(1,q_2-1)$ $p_{1,2}(2,1)$ $p_{1,2}(2,2)$ ... $p_{1,2}(q_1-1,q_2-1)$
> $p_{1,3}(1,1)$ ...

In other words, the first $N$ lines of the file record the frequency that each state is observed at each site, and the next $N(N-1)/2$ lines record the pairwise correlations. Note that, because $\sum_{a=1}^{q_i} p_i(a)=1$, the frequency (and corresponding pair correlations) for one state at each site need not be specified explicitly.

These values should be given in floating point or scientific format, with whitespace (e.g. `'\t'`) between successive values and a newline character (`'\n'`) at the end of each line. In order for the correlations to be read in properly, there should be **no** whitespace between the final correlation value and the newline character on each line. 

 For examples, see the `examples/` directory. Instructions on how to automatically generate a correlations file from a sequence alignment in FASTA format using Matlab can be found [here](#generating-correlations-from-a-sequence-alignment).


# Running the program

Here we show a simple example of how to run the program and interpret the output, using a set of sample data for the HIV protein p6. Full explanations for the possible options are given [here](#command-line-options).

### Running ACE

We begin running the ACE algorithm on the example p6 dataset with the command:

```bash
$ ./bin/sce -d examples -i p6 -o p6-out -b 4064.0 -g2 2.5e-4 -kmax 6
```

This creates two new files, `examples/p6-out.sce` and `examples/p6-out.j`, which record general output on the inference procedure and the current inferred Potts parameters, respectively.

Output from the first file, `examples/p6-out.sce`, should appear something like the following:

>1.000000e+00 &nbsp;&nbsp;3.189915e-01 &nbsp;&nbsp; 8.262700e-01 &nbsp;&nbsp; 5.258494e+01 &nbsp;&nbsp; 1.604867e+01 &nbsp;&nbsp; 2 &nbsp;&nbsp; 1378 &nbsp;&nbsp; 52
>9.523810e-01 &nbsp;&nbsp; 3.189915e-01 &nbsp;&nbsp; 8.262700e-01 &nbsp;&nbsp; 5.258494e+01 &nbsp;&nbsp; 1.604867e+01 &nbsp;&nbsp; 2 &nbsp;&nbsp; 1378 &nbsp;&nbsp; 52
>9.070295e-01 &nbsp;&nbsp; 3.189915e-01 &nbsp;&nbsp; 8.262700e-01 &nbsp;&nbsp; 5.258494e+01 &nbsp;&nbsp; 1.604867e+01 &nbsp;&nbsp; 2 &nbsp;&nbsp; 1378 &nbsp;&nbsp; 52
>8.638376e-01 &nbsp;&nbsp; 3.189915e-01 &nbsp;&nbsp; 8.262700e-01 &nbsp;&nbsp; 5.258494e+01 &nbsp;&nbsp; 1.604867e+01 &nbsp;&nbsp; 2 &nbsp;&nbsp; 1378 &nbsp;&nbsp; 52
>8.227025e-01 &nbsp;&nbsp; 3.189915e-01 &nbsp;&nbsp; 8.262700e-01 &nbsp;&nbsp; 5.258494e+01 &nbsp;&nbsp; 1.604867e+01 &nbsp;&nbsp; 2 &nbsp;&nbsp; 1378 &nbsp;&nbsp; 52
>7.835262e-01 &nbsp;&nbsp; 3.189915e-01 &nbsp;&nbsp; 8.262700e-01 &nbsp;&nbsp; 5.258494e+01 &nbsp;&nbsp; 1.604867e+01 &nbsp;&nbsp; 2 &nbsp;&nbsp; 1378 &nbsp;&nbsp; 52
>7.462154e-01 &nbsp;&nbsp; 3.189915e-01 &nbsp;&nbsp; 8.262700e-01 &nbsp;&nbsp; 5.258494e+01 &nbsp;&nbsp; 1.604867e+01 &nbsp;&nbsp; 2 &nbsp;&nbsp; 1378 &nbsp;&nbsp; 52
>...

These columns represent, respectively: the current value of the threshold $\theta$, error on the one-point correlations $\epsilon_{p1}$, error on the pairwise correlations $\epsilon_{p2}$, normalized maximum error $\epsilon_{\rm max}$, current estimate of the entropy $S$, maximum cluster size, total number of clusters in the expansion, and the number of selected clusters (i.e. those for which $| \Delta S |>\theta$).

The inferred Potts parameters in the second file, `examples/p6-out.j`, are output in the same format as the input correlations, as shown [above](#required-input). In this case, the first $N$ lines record the Potts fields $h_i(a)$, and the following $N(N-1)/2$ lines record the couplings $J_{ij}(a,b)$.

With the command line options entered above, the program will terminate when the maximum cluster size reaches 6. The final line of `examples/p6-out.sce` should then appear something like:

> 1.274301e-05 &nbsp;&nbsp; 5.750033e+02 &nbsp;&nbsp; 4.429895e+00	 &nbsp;&nbsp; 4.301941e+06 &nbsp;&nbsp; 1.450155e+01 &nbsp;&nbsp; 6 &nbsp;&nbsp; 29478 &nbsp;&nbsp; 6194

As you can see, the current Potts parameters do not accurately reproduce the input correlations (the error terms $\epsilon_{p1}, \epsilon_{p2}, \epsilon_{\rm max}>1$). However, the entropy has nearly converged (see column 6 in `examples/p6-out.sce`), suggesting that the inferred Potts parameters may be very close to ones that do recover the correlations accurately. Our next step is to run the <b>M</b>onte <b>C</b>arlo (MC) learning algorithm to refine the inferred Potts parameters so that they reproduce the input correlations.

### Running the MC learning algorithm QLS

We now run the MC algorithm on the output we previously obtained from ACE, using the command:

```bash
$ ./bin/qls -d examples -i p6-out -o p6-out-learn -b 4064.0 -g2 2.5e-4 -c p6
```

This creates two additional output files, `examples/p6-out-learn.fit` and `examples/p6-out-learn.j`, which record progress on the MC learning procedure and the current refined Potts parameters, respectively.

Output from the first file, `examples/p6-out-learn.fit`, should appear something like the following:

>1 &nbsp;&nbsp; 5.776764e+02 &nbsp;&nbsp; 8.605859e+02 &nbsp;&nbsp; 4.342585e+06 &nbsp;&nbsp; 1.900000e+00

>2 &nbsp;&nbsp; 5.768460e+02 &nbsp;&nbsp; 8.592596e+02 &nbsp;&nbsp; 4.330094e+06 &nbsp;&nbsp; 3.610000e+00
>3 &nbsp;&nbsp; 5.773766e+02 &nbsp;&nbsp; 8.598271e+02 &nbsp;&nbsp; 4.338444e+06 &nbsp;&nbsp; 6.859000e+00
>4 &nbsp;&nbsp; 5.772833e+02 &nbsp;&nbsp; 8.593048e+02 &nbsp;&nbsp; 4.337158e+06 &nbsp;&nbsp; 1.303210e+01
>5 &nbsp;&nbsp; 5.764251e+02 &nbsp;&nbsp; 8.572838e+02 &nbsp;&nbsp; 4.324644e+06 &nbsp;&nbsp; 2.476099e+01
>6	 &nbsp;&nbsp; 5.775672e+02 &nbsp;&nbsp; 8.575524e+02 &nbsp;&nbsp; 4.342923e+06 &nbsp;&nbsp; 4.704588e+01
>7 &nbsp;&nbsp; 5.765323e+02 &nbsp;&nbsp; 8.534015e+02 &nbsp;&nbsp; 4.329062e+06 &nbsp;&nbsp; 8.938717e+01
>...

These columns represent, respectively: the current iteration, error on the one-point correlations $\epsilon_{p1}$, error on the pairwise correlations $\epsilon_{p2}$, normalized maximum error $\epsilon_{\rm max}$, and the maximum size of the weight parameter used in the MC learning update step.

After about **XX** iterations, the MC learning algorithm should converge and the program will terminate. The Potts parameters recorded in the second file, `examples/p6-out-learn.j`, now specify a model that accurately recovers the input correlations to within fluctuations expected due to finite sampling.

### Verifying the output with QGT

To be added.

# Command line options

### Options for both ACE and QLS

- `-d` gives the path to the directory where data files are located, and where output will be written
- `-i` specifies the name of the input file (excluding the extension)
- `-o` specifies the name of the output file (excluding the extension)
- `-v` enables verbose output
- `-b` tells the program how many samples were used to generate the input correlations, so that the expected error in the correlations due to finite sampling can be estimated
- `-mcb` gives the number of Monte Carlo steps used to estimate the inference error
- `-mcr` gives the number of independent Monte Carlo trajectories to use when estimating the inference error
- `-g2` sets the $L_2$-norm regularization strength (note that a natural value for this parameter is $1/B$, where $B$ is the number of samples used to generate the input correlations)
- `-ag` automatically sets the $L_2$-norm regularization strength equal to $1/B$, using the number of samples $B$ passed with the `-b` option

### Additional ACE options

- `-kmin` sets the minimum cluster size required before the program will terminate
- `-kmax` sets the maximum cluster size; the program terminates automatically after a cluster of this size is created
- `-t` specifies a single value of the threshold $\theta$ at which the algorithm will run, then exit
- `-tmax` specifies the maximum (starting) value of the threshold
- `-tmin` specifies the minimum allowed value of the threshold; the program terminates automatically after $\theta$ falls below this minimum value
- `-ts` specifies the logarithmic step size to for between successive values $\theta$, through $\theta_{i+1} = \theta_i / \theta_{\rm step}$
- `-r` enables the expansion of the entropy $S$ around a mean-field reference entropy $S_{0}$, which may be helpful in particular for inferring models described by dense networks of weak interactions (note: works only if all variables are binary)
- `-g0` sets the $L_0$-norm regularization strength, and turns on $L_0$-norm regularization, enforcing sparsity for Potts couplings
- `-l0` turns on $L_0$-norm regularization, but **without** setting the regularization strength
- `-ss` specifies an input "secondary structure" file used to specify the initial set of clusters to consider in the expansion
- `-lax` enables a laxer cluster construction rule, increasing the number of clusters included in the cluster expansion routine

### Additional QLS options

- `-c` specifies the set of (true) correlations to compare with for the MC learning routine
- `-e` sets the maximum tolerable error threshold; the program will run until all of the error terms $\epsilon_{p1}, \epsilon_{p2}, \epsilon_{\rm max}<e$


# Generating correlations from a sequence alignment

To be added.

# References

1. ACE reference here.
2. [Cocco, S. and Monasson, R. (2011). Adaptive Cluster Expansion for Inferring Boltzmann Machines with Noisy Data. <i>Physical Review Letters</i>, <b>106</b>, 090601][2].
[2]: http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.106.090601

3. [Cocco, S. and Monasson, R. (2012). Adaptive Cluster Expansion for the Inverse Ising Problem: Convergence, Algorithm and Tests.<i> Journal of Statistical Physics</i>, <b>147</b>(2), 252–314][3].
[3]: http://link.springer.com/article/10.1007/s10955-012-0463-4#page-1

4. [Barton, J. and Cocco, S. (2013). Ising models for neural activity inferred via selective cluster expansion: structural and coding properties. <i>Journal of Statistical Mechanics: Theory and Experiment</i>, <b>2013</b>(03), P03002][4].
[4]: http://iopscience.iop.org/1742-5468/2013/03/P03002


> Written with [StackEdit](https://stackedit.io/).