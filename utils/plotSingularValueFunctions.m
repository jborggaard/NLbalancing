function [] = plotSingularValueFunctions(sigma,c,zRange,n,maxDegree)
%  Plots polynomial approximations to singular value functions.
%
%   plotSingularValueFunctions(sigma,c,zRange,n,l)
%
%   Given an approximation to the singular value functions with the form
%
%       xi_i(z) = sigma(i) + sum_{k=1}^l diag(c{k})*z.^k
%
%   this function plots the above approximation for i=1,...,n at the points
%   specified in zRange. 
% 
%   The default value of n (the number of singular value functions plotted)
%   is determined from the length of sigma.
% 
%   The default value of the maximum degree of the approximation, l, is 
%   determined from the length of the cell array c.
%
%   Note that we explicitly assume the same range of z in all of the
%   coordinate directions.
%
%   Author: Jeff Borggaard, Virginia Tech
%
%   License: MIT
%
%   Reference: Nonlinear balanced trunction: Part 2--Model reduction on
%              manifolds, by Kramer, Gugercin, and Borggaard, arXiv.
%
%   Part of the NLbalancing repository.
%%
  if (nargin<4)
    n = length(sigma);
  end
  if (nargin<5)
    maxDegree = length(c);
  end

  legendText = cell(1,n);
  for i=1:n
    xi = sigma(i)*ones(size(zRange));
    for k=1:maxDegree
      xi = xi + c{k}(i)*zRange.^k;
    end
    legendText{i} = sprintf('xi_{%d}',i);
    plot(zRange,xi); hold on
  end
  legend(legendText)


end