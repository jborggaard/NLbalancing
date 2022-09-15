function [A,B,C,N] = getSystem2()
%getSystem2  Generates a quadratic system for testing energy functions.
%
%        \dot(x1) = -x1 + x2 - x2^2 + u1
%        \dot(x2) =     - x2        + u1
%              y1 =  x1 + x2
%
%%

  A = [-1 1;0 -1];
  N = [0 0 0 -1;0 0 0 0];
  B = [1;1];
  C = [1 1];

end