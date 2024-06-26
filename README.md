# NLbalancing
Software that uses polynomials to approximately solve the nonlinear (NL) balancing problem for systems with quadratic nonlinearities.  The description of the NL balancing algorithms are provided in the papers

- *Nonlinear Balanced Truncation: Part 1-computing energy functions*

- *Nonlinear Balanced Truncation: Part 2-nonlinear manifold model reduction*

by Boris Kramer, Serkan Gugercin, and Jeff Borggaard.  The Kronecker product
solvers are based on those in the KroneckerTools repository and also used for
the QQR software described in

- *Approximating Polynomial-Quadratic Regulator Problems, Arxiv*

by Jeff Borggaard and Lizette Zietsman (full references included below).

## Installation Notes
Clone this repository: 
```
  git clone https://www.github.com/jborggaard/KroneckerTools.git
  git clone https://www.github.com/jborggaard/NLbalancing.git
```
then modify the path in **setKroneckerToolsPath.m**

The installation can be tested in Matlab (we used R2020b) by typing
```
>> examplesForPaper1
```
and
```
>> examplesForPaper2
```
that provide the numerical results for our nonlinear balanced truncation papers.

The details of some of our functions and test examples are provided below.  


## How to use this software

Let A be n-by-n, B be n-by-m, C be p-by-n, and N be n-by-n^2 , with [A,B] a controllable pair and [A,C] a detectable pair.  The parameter eta is used we can compute the coefficients of the solution to future and past energy functions in Matlab as
```
>>  [w] = approxFutureEnergy(A,N,B,C,eta,degree);

and

>>  [v] = approxPastEnergy(A,N,B,C,eta,degree);
```
The variable _v_ is a cell array with _v{2}_ being n-by-n^2 , up to _v{degree+1}_ which is n-by-n^(degree+1) .  These are coefficients of the polynomial approximation to the value function.  From an initial _x0_, we can compute the approximation to the energy function as
```
>>  E = (1/2)*( v{2}*kron(x0,x0) + ... + v{degree+1}*kron(kron(... ,x0),x0) );
```
or, using the utility function,
```
>>  E = (1/2)*kronPolyEval(v(1:degree),x0,degree);
```

For details on how to compute input-normal balancing with **inBalance**, type
```
>>  help inBalance
```
Examples of input-normal balancing are found in the examples folder (inExample1 and inExample2).

The **inBalance** function uses **inputNormalTransformation** and **approximateSingularValueFunctions** with a provided tolerance to build a balanced reduced-order model.

For details on how to run **HJBbalance**, type
```
>>  help HJBbalance
```

for examples how to run **HJBbalance** see those in
```
>> examplesForPaper3
```
and the files in the examples directory.


## Description of Files
#### setKroneckerToolsPath

Defines the path to the KroneckerTools directory containing functions for working with Kronecker product expressions.  KroneckerTools can be downloaded from github.com/jborggaard/KroneckerTools  The default assumes that NLbalancing and KroneckerTools lie in the same directory and uses relative pathnames.  This should be changed if you use different locations.  (setKroneckerToolsPath also lies in the examples and tests directories, so should be changed there as well if you plan to run functions from those directories.

#### CT2Kron and Kron2CT

These compute mappings between coefficients of a multidimensional polynomial in compact Taylor series format and those in a Kronecker product format.  As a simple example, if p(x) = c1 x1^2 + c2 x1 x2 + c3 x2^2 , then n=2, degree=2.  We have

p(x) = [c1 c2 c3] * [x1^2 x1x2 x2^2 ].' written as

p(x) = ( CT2Kron(n,degree)*[c1 c2 c3].' ).' * kron([x1;x2],[x1;x2])

or

p(x) = [c1 c2/2 c2/2 c3] * kron([x1;x2],[x1;x2]) written as

p(x) = ( Kron2CT(n,degree) * [c1 c2/2 c2/2 c3].' ).' * [x1^2 x1x2 x2^2 ].'

There mappings are used to balance coefficients of the feedback and value functions.  (e.g., in the Kronecker form, we seek the same coefficient for x1 x2 and x2 x1).  This is done automatically using the provided function, kronPolySymmetrize.

#### LyapProduct

Efficiently computes the product of a special Kronecker sum matrix (aka an N-Way Lyapunov matrix) with a vector.  This is done by reshaping the vector, performing matrix-matrix products, then reshaping the answer.  We could also utilize the matrization of the associated tensor.

## Examples

### Example01.m

Approximates the future and past energy functions for a one-dimensional model problem motivated by the literature.  This appears as example 1 in Kramer, Borggaard, and Gugercin.

### Example02.m

Approximates the future and past energy functions, then computes an approximation to the (input-normal) balancing transformation and computes a reduced model.  The example is based on a two-dimensional problem found in Kawano and Scherpen, IEEE Transactions on Automatic Control, 2016 (we ignore their bilinear term in this example).


## Algorithms from Kramer, Gugercin, and Borggaard, Part 1:

### Algorithm 1 is implemented in _approxFutureEnergy.m_ and _approxPastEnergy.m_

## Algorithms from Kramer, Gugercin, and Borggaard, Part 2:

### Algorithm 1 is implemented in _inputNormalTransformation.m_

### Algorithm 2 is implemented in _approximateSingularValueFunctions.m_

## References
```
  @misc{kramer2022balancedtruncation1,
    title={Nonlinear Balanced Truncation: Part 1-computing energy functions},
    author={Boris Kramer, Jeff Borggaard, and Serkan Gugercin},
    year={2022},
    eprint={2209.07645},
    archivePrefix={arXiv},
    primaryClass={math.OC}
  }
```

```
  @misc{kramer2022balancedtruncation2,
    title={Nonlinear Balanced Truncation: Part 2-nonlinear manifold model reduction},
    author={Boris Kramer, Jeff Borggaard, and Serkan Gugercin},
    year={2022},
    eprint={pending},
    archivePrefix={arXiv},
    primaryClass={math.OC}
  }
```

```
  @inproceedings{borggaard2019quadraticquadratic,
    title={The Quadratic-Quadratic Regulator Problem: 
     Approximating feedback controls for quadratic-in-state nonlinear systems},
    author={Jeff Borggaard and Lizette Zietsman}, 
    booktitle={Proceedings of the 2020 American Conference on Control},
    year={2020},
    eprint={1910.03396},
    archivePrefix={arXiv},
    primaryClass={math.OC}
  }
```

```
  @article{borggaard2021polynomialquadratic,
    title={On Approximating Polynomial-Quadratic Regulator Problems},
    author={Jeff Borggaard and Lizette Zietsman},
    journal={IFAC-PaersOnLine},
    volume=54,
    number=9,
    pages={329--334},
    year={2021}
  }
```


