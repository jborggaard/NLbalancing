function [K,PI,E] = hamiltonian(A,B,Q,R,destabilizing)
%-------------------------------------------------------------------------------
%  HAMILTONIAN  Solve the Riccati equation using eigenvectors of Hamiltonian
%
%    This version optionally finds the destabilizing solution.
%
%    [K,PI,E] = HAMILTONIAN(A,B,Q,R,flag) solves the algebraic Riccati equation,
%    A'*PI + PI*A + PI*B*inv(R)*B'*PI + Q = 0, by forming the Hamiltonian matrix
%    M = [A  -B*inv(R)*B' ;  -Q  -A'], then splitting the eigenvectors into 
%    those associated with positive and negative eigenvalues.  The solutions to
%    the Riccati equation are then given by either PI = V21*inv(V11) or
%    PI = V22*inv(V12).  
%
%    if   flag = false (default), the stabilizing solution is returned
%
%    else flag = true             returns the destabilizing solution
%
%
%    To match the output of Matlab's LQR, we also return the values 
%       K = inv(R)*B'*PI 
%    and 
%       E = eig(A-B*K)  (optionally).
% 
%    Author: Jeff Borggaard, Virginia Tech
%            part of the QQGbalancing library.
%
%    License: MIT public license
%    
%    See also <a href="matlab:help lqr">lqr</a>
%-------------------------------------------------------------------------------

  if ( nargin<5 )
    destabilizing = false;
  end

  [n,m] = size(B);

  RinvBT = R\B';

  % Form the Hamiltonian
  H = [A -B*RinvBT; -Q -A'];

  [U,T] = schur(H,'complex');

  if ( destabilizing )
    [U,T] = ordschur(U,T,'rhp');
  else
    [U,T] = ordschur(U,T,'lhp');
  end

  [V,e] = eig(T);
  U = U*V;
  
  PI = U( n+1:2*n, 1:n )/( U( 1:n, 1:n ) );
  K  = RinvBT*PI;

  if ( nargout==3 )
    E  = eig(A-B*K);
  end
  
  %  test the residual?
  % norm(A'*PI+PI*A-PI*B*RinvBT*PI+Q)
end % function hamiltonian