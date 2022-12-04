function [T,c,sigma,v,w] = inBalance(A,N,B,C,eta,maxDegree)
%inBalance Creates input-normal balancing transformation for quadratic systems.
%   Given a system of the form
%       \dot{x} = A*x + N*kron(x,x) + B*u
%             y = C*x,
%
%   the parameters eta and maxDegree, produces the polynomial
%   transformation coefficients T of the transformation
%       x = T(z) == T{1}*z + T{2}*kron(z,z) + ... + T{maxDegree}*kron(x,...),x)
%
%   as well as the coefficients of the singular value functions c 
%       c(z) = sigma + c{1}*z + c{2}*kron(z,z) + ... + c{maxDegree-1}*kron(z,...)
%
%   such that the past energy function 
%       v(x) = 0.5*(v{1}*x + v{2}*kron(x,x) + ... + v{maxDegree+1}*kron(x,...),x))
%
%   and the future energy function
%       w(x) = 0.5*(w{1}*x + w{2}*kron(x,x) + ... + w{maxDegree+1}*kron(x,...),x))
%
%   have the property
%       v(T(z)) = 0.5*z.'*z     and   w(T(z)) = 0.5*z.'*sigma.^2*z.
%
%   Usage:
%       [T,c,sigma,v,w] = inBalance(A,N,B,C,eta,maxDegree)
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  License: MIT
%
%  References: Nonlinear balanced truncation: Part 1--Computing energy functions,
%              by Kramer, Gugercin, and Borggaard, arXiv:2209.07645.
%
%              See Algorithm 1.
%
%              Nonlinear balanced truncation: Part 2--Nonlinear manifold model 
%              reduction, by Kramer, Gugercin, and Borggaard, arXiv.
%
%              See Algorithms 1 and 2.
%
%  Part of the NLbalancing repository.
%%


%% approximate the energy functions 
  [w] = approxFutureEnergy(A,N,B,C,eta,maxDegree+1);
  [v] = approxPastEnergy(A,N,B,C,eta,maxDegree+1);

  %% Approximate the Transformations using Algorithm 1.
  [sigma,T] = inputNormalTransformation(v,w,maxDegree);

  %% Approximate the singular value functions using Algorithm 2.
  [c] = approximateSingularValueFunctions(T,w,sigma,maxDegree-1);

end