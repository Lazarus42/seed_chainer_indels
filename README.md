# Seed-chain-extend simulations on a random string passed through indel and substitution channels.

Experiments on seed-chain-extend on random strings by implementing a basic seed-chain-extend algorithm in rust. The program generates a random string `S` with letters chosen i.i.d. from {A,C,T,G} and obtains `S'` by passing a substring of `S` through a mutation channel where independently at each position the following may occur: (1) a substitution, (2) a deletion of that position, (3) an insertion of a uniformly random string of geometric length inserted to the left. The strings `S` and `S'` are aligned with seed-chain-extend. For each command-line run, the average recoverability and empirical runtime (as defined in the paper) are recorded to the `recoverability.csv` file.

The value of the `k`-mers is increasing as `k = C log n` where `n` is the sequence length, and `C` is defined to be `3/(1 - 2 * alpha)` where `alpha = -log_4 (1 - theta)`. We simulate alignments for `k = 20, 21, ...` up to a user specified value. Other simulation parameters can be specified and are outlined below. For a quick overview of the algorithm:

**Seeding**: using open syncmer or minimizer seeds. NOTE: We don't use any sort of bitwise algorithms for representing k-mers. We have not optimized for k-mer matching, seeding, etc so it will be slow.  

**Chaining**: using a linear gap cost with a MRQ data structure as described in the minigraph paper.

**Extension**: using rust-bio's simple dynamic programming extension algorithm. 

### Requirements 

1. [rust](https://www.rust-lang.org/tools/install) and associated tools such as cargo are required and assumed to be in PATH.
### Install
```
git clone https://github.com/Lazarus42/seed_chainer_indels.git
cd seed_chainer_indels
cargo build --release
./target/release/basic_seed_chainer 100 0.10 -k 45 --substring
#100 iterations at theta 0.10 for k = 20,...,45
```

### Parameters

One can specify a few parameters for the simulation described by the `basic_seed_chainer -h` menu. The following are supported: 
1. `-k INT`: maximum `k`-mer size to iterate up to.
2. `num_iters`: you must specify the number of iterations for each value of `k`. For example: if the command run is `./target/release/basic_seed_chainer 100 0.10 -k 45 --substring`, then 100 iterations occur for each value of `k`.
3. `theta_T`: you must specify the total mutation rate `theta_T`, which is the sum of the insertion, substitution, and deletion rates (`theta_T = theta_i + theta_s + theta_d`). Example for `theta_T = 0.10`: `./target/release/basic_seed_chainer 100 0.10 -k 45 --substring`
4. To replicate the experiments and results used in the paper, be sure to include the `--substring` parameter.

### Results
For each value of `k` between 20 and the maximum `k` specified by the user, the tuple:

`(|S|, average length of |S'|, k, average recoverability, average empirical runtime of chaining + extension) `

is written to the `recoverability.csv` file. Note: the average is computed over the number of iterations, `num_iters`, specified by the user.

Additionally, the values of `k, n = |S|, m' = length of substring passed through mutation channel, lower bound for average recoverability, average empirical runtime of extension, and  average empirical runtime of chaining` are printed in the following way:

```
Theta = 0.1, alpha = 0.07600154672252499, C = 3.5377487545181285
n = 2542,  m' = 415
Value for k 20
Mean lower bound for recoverability for k 20 is 0.808363145421188
Mean extend time 0.00033765155
Mean chain time 0.000063527346
n = 3758,  m' = 560
Value for k 21
Mean lower bound for recoverability for k 21 is 0.817117074973988
Mean extend time 0.0006030821
Mean chain time 0.00007562292
...
```


## Plotting results

The resulting `recoverability.csv` file can be analyzed to obtain graphs appearing in the paper with the `sce_simulation_analysis.ipynb` file found in the `analysis_code` folder.
