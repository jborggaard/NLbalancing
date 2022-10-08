function [T1,Sigma] = linearBalancing(A,B,C)
%linearBalancing solves the linear balancing problem for comparison.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  License: MIT
%
%  Reference:  Nonlinear balanced truncation: Part 2--Model reduction on
%              manifolds, by Kramer, Gugercin, and Borggaard, arXiv.
%            
%  Part of the NLbalancing library.

  %  Compute the infinite reachability Gramian.  This can be
  %  used to compute the minimal energy required to drive the
  %  state from 0 to x as x'inv(P)x where P solves
  %    AP + PA' + BB' = 0.
  P = lyap(A,B*B');
  
  %  Compute the infinite observability Gramian.  This can be used
  %  to compute the maximal energy produced by observing the 
  %  output of a system y=C*x from the initial state x as
  %  x'*Q*w where Q solves
  %    A'Q + QA + C'C = 0
  Q = lyap(A',C'*C);
  
  L = chol(P,'lower');      % P = LL'
  R = chol(Q,'lower');      % Q = RR'
  
  [U,Sigma,V] = svd(L'*R);  % L'QL = U Sigma^2 U'
  
  T1 = sqrt(Sigma)*U'/L;    % then T1*P*T1' = inv(T1)'*Q*inv(T1) = Sigma
  
end

