function [w] = approxFutureEnergy(A,N,B,C,eta,d,verbose)
%  Calculates a polynomial approximation to the future energy function
%  for a quadratic system.
%
%  w = approxFutureEnergy(A,N,B,C,eta,d) 
%
%  Computes a degree d polynomial approximation to the future energy function 
%
%          E^+(x) = 1/2 ( w{2}'*kron(x,x) + ... + w{d}'*kron(.. x) )
%
%  for the polynomial system
%
%    \dot{x} = Ax + Bu + N*kron(x,x)
%          y = Cx
%
%  where eta = 1-gamma^(-2), gamma is the parameter in the algebraic Riccati
%  equation for W2,
%
%    A'*W2 + W2*A - eta*W2*B*B'*W2 + C'*C = 0.
%
%  Note that w2 = vec(W2).  Details are in Section III.B of the reference.
%
%  Requires the following functions from the KroneckerTools repository:
%      KroneckerSumSolver
%      kronMonomialSymmetrize
%      LyapProduct
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  License: MIT
%
%  Reference: Nonlinear balanced truncation: Part 1--Computing energy functions,
%             by Kramer, Gugercin, and Borggaard, arXiv:2209.07645.
%
%             See Algorithm 1.
%
%  Part of the NLbalancing repository.
%%

  if (nargin<7)
    verbose = false;
  end

  n = size(A,1);   % A should be n-by-n
  m = size(B,2);   % B should be n-by-m
  p = size(C,1);   % C should be p-by-n

  Acell = cell(1,d);
  for i=1:d
    Acell{i} = A.';
  end

  % Create a vec function for readability
  vec = @(X) X(:);


  %% k=2 case
  R = eye(m)/eta;

  if ( eta>0 )
    [W2] = icare(A,B,(C.'*C),R);
    if ( isempty(W2) && verbose )
      warning('approxFutureEnergy: icare couldn''t find stabilizing solution')
    end
    
  elseif ( eta<0 )
    [W2] = icare(A,B,(C.'*C),R,'anti');
    
    if ( isempty(W2) && verbose )
      warning('approxFutureEnergy: icare couldn''t find stabilizing solution')
      warning('approxFutureEnergy: using the hamiltonian')
      [~,W2,~] = hamiltonian(A,B,C.'*C,R,true);
    end
    
  else % eta==0
    [W2] = lyap(A.',(C.'*C));
    
  end

  if (isempty(W2))
    error('approxFutureEnergy: Can''t find a stabilizing solution')
  end
  
  %  Check the residual of the Riccati/Lyapunov equation
  if (verbose)
    RES = A'*W2 + W2*A - eta*(W2*B)*(B'*W2) + C'*C ;
    fprintf('The residual of the Riccati equation is %g\n',norm(RES,'inf'));
  end

  %  Reshape the resulting quadratic coefficients
  w2 = W2(:);
  w{2} = w2;

  %% k=3 case
  if ( d>2 )
    b = -LyapProduct(N.',w2,2);
    W2BB = W2*(B*B.');
    [w{3}] = KroneckerSumSolver(Acell(1:3),b,2,-3*eta*W2BB);
    
    [w{3}] = kronMonomialSymmetrize(w{3},n,3);
  end
  
  %% k>3 cases (up to d)
  BWk = cell(1,d-1);
  if ( d>3 )

    for k=4:d
      BWk{k-1} = B.'*reshape(w{k-1},n,n^(k-2));

      b = -LyapProduct(N.',w{k-1},k-1);
      
      for i=3:(k+1)/2
        j   = k+2-i;
        tmp = BWk{i}.'*BWk{j};
        b   = b + 0.25*eta*i*j*(vec(tmp) + vec(tmp.'));
      end

      if (mod(k,2)==0) % k is even
        i   = (k+2)/2; 
        j   = i;
        tmp = BWk{i}.'*BWk{j};
        b   = b + 0.25*eta*i*j*vec(tmp);
      end

      [w{k}] = KroneckerSumSolver(Acell(1:k),b,k-1,-k*eta*W2*(B*B.'));

      [w{k}] = kronMonomialSymmetrize(w{k},n,k);
    end
  end

end
