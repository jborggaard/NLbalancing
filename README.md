# NLbalancing
Software that uses polynomials to approximately solve the nonlinear (NL) balancing problem for systems with quadratic nonlinearities.  The description of the NL balancing algorithms are provided in the papers

- *Balanced Truncation Model Reduction for Large-Scale Polynomial Systems*

by Boris Kramer, Jeff Borggaard, and Serkan Gugercin.  The Kronecker product
solvers are based on those in the QQR software described in

- *Approximating Polynomial-Quadratic Regulator Problems, Arxiv*

by Jeff Borggaard and Lizette Zietsman (full references included below).

## Installation Notes
Clone this repository: 
```
  git clone https://www.github.com/jborggaard/NLbalancing.git
```

The installation can be tested in Matlab (we used R2020b) by typing
```
>> examplesForPaper
```

The details of some of our functions and test examples are provided below.  


## How to use this software

Let A be n-by-n, B be n-by-m, C be p-by-n, and N be n-by-n^2 , with [A,B] a controllable pair and [A,C] a detectable pair.  The parameter eta is used we can compute the coefficients of the solution to future and past energy functions in Matlab as
```
>>  [v] = solveFutureEnergy(A,N,B,C,eta,degree);

and

>>  [w] = solvePastEnergy(A,N,B,C,eta,degree);
```
The variable _v_ is a cell array with _v{2}_ being n-by-n^2 , up to _v{degree+1}_ which is n-by-n^(degree+1) .  These are coefficients of the polynomial approximation to the value function.  From an initial _x0_, we can compute the approximation to the energy function as
```
>>  E = v{2}*kron(x0,x0) + ... + v{degree+1}*kron(kron(... ,x0),x0);
```
(to do: create an efficient function that evaluates E)

For details on how to run **HJBbalance**, type
```
>>  help HJBbalance
```

for examples how to run **HJBbalance** see those in
```
>> examplesForPaper
```
and found in the examples directory


## Description of Files
#### setKroneckerSumPath

Defines the path to the functions for working with Kronecker product expressions.

#### CT2Kron and Kron2CT

These compute mappings between coefficients of a multidimensional polynomial in compact Taylor series format and those in a Kronecker product format.  As a simple example, if p(x) = c1 x1^2 + c2 x1 x2 + c3 x2^2 , then n=2, degree=2.  We have

p(x) = [c1 c2 c3] * [x1^2 x1x2 x2^2 ].' written as

p(x) = ( CT2Kron(n,degree)*[c1 c2 c3].' ).' * kron([x1;x2],[x1;x2])

or

p(x) = [c1 c2/2 c2/2 c3] * kron([x1;x2],[x1;x2]) written as

p(x) = ( Kron2CT(n,degree) * [c1 c2/2 c2/2 c3].' ).' * [x1^2 x1x2 x2^2 ].'

There mappings are used to balance coefficients of the feedback and value functions.  (e.g., in the Kronecker form, we seek the same coefficient for x1 x2 and x2 x1).

#### LyapProduct

Efficiently computes the product of a special Kronecker sum matrix (aka an N-Way Lyapunov matrix) with a vector.  This is done by reshaping the vector, performing matrix-matrix products, then reshaping the answer.  We could also utilize the matrization of the associated tensor.

## Examples

### Example01.m

Approximates the future and past energy functions for a one-dimensional model problem motivated by the literature.  This appears as example 1 in Kramer, Borggaard, and Gugercin.



### References
```
  @misc{kramer2021balancedtruncation,
    title={Balanced Truncation Model Reduction for Large-Scale Polynomial Systems},
    author={Boris Kramer, Jeff Borggaard, and Serkan Gugercin},
    year={2021},
    eprint={pending},
    archivePrefix={arXiv},
    primaryClass={math.OC}
  }
```

```
  @misc{borggaard2019quadraticquadratic,
    title={The Quadratic-Quadratic Regulator Problem: 
     Approximating feedback controls for quadratic-in-state nonlinear systems},
    author={Jeff Borggaard and Lizette Zietsman}, 
    note={to appear in the Proceedings of the 2020 American Conference on Control},
    year={2019},
    eprint={1910.03396},
    archivePrefix={arXiv},
    primaryClass={math.OC}
  }
```

```
  @misc{borggaard2020polynomialquadratic,
    title={On Approximating Polynomial-Quadratic Regulator Problems},
    author={Jeff Borggaard and Lizette Zietsman},
    year={2020},
    note={submitted}
  }
```


