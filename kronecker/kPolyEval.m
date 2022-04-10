function [x] = kPolyEval(c,z,degree)
%kPolyEval evaluates a kronecker polynomial out to degree=degree
%
%  Usage:
%     x = kPolyEval(c,z,degree);
%
%  Input Variables:
%     c - coefficients of the polynomial (cell array)
%     z - value to calculate the polynomial at (vector)
%     degree - polynomial will be evaluated out to degree=degree
%              (default is the length of c)
%
%  Returns:  x = c{1}*z + c{2}*kron(z,z) + ... c{degree}*kron(...)
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Part of the NLbalancing repository.
%%
  
  d = length(c);
  if (nargin<3)
    degree = d;
  else
    if (d<degree)
      error('kPolyEval: not enough entries in the coefficient cell array')
    end
  end
  
  if isempty(c{1})
    n = size(c{2},1);
    x = zeros(n,1);
  else
    x = c{1}*z;
  end
  
  zk = z;
  for k=2:degree
    zk = kron(zk,z);
    x  = x + c{k}*zk;
  end
  
end % function kPolyEval
  


