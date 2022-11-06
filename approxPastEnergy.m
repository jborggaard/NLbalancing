function [v] = approxPastEnergy(A,N,B,C,eta,d,verbose)
%  Calculates a polynomial approximation to the past energy function
%  for a quadratic system.
%
%  v = approxPastEnergy(A,N,B,C,eta,d) 
%
%  Computes a degree d polynomial approximation to the past energy function 
%
%          E^-(x) = 1/2 ( v{2}'*kron(x,x) + ... + v{d}'*kron(.. x) )
%
%  for the polynomial system
%
%    \dot{x} = Ax + N*kron(x,x) + Bu
%          y = Cx
%
%  where eta = 1-gamma^(-2), gamma is the parameter in the algebraic Riccati
%  equation
%
%    A'*V2 + V2*A + V2*B*B'*V2 - eta*C'*C = 0.
%
%  Note that v{2} = vec(V2) = V2(:).  Details are in Section III.C of the paper.
%
%  Requires functions from the KroneckerTools repository
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
  R = eye(m);
  
  if ( eta~=0 )
    %  We multiply the ARE by -1 to put it into the standard form in icare,
    %  and know (-A,-B) is a controllable pair if (A,B) is.
    [V2] = icare(-A,-B,eta*(C.'*C),R);

    if ( isempty(V2) && verbose )
      warning('approxPastEnergy: icare couldn''t find stabilizing solution')
    end
    
    if ( isempty(V2) )
      if (verbose) 
        warning('approxPastEnergy: using matrix inverse')
      end
      [Yinf] = icare(A.',C.',B*B.',eye(p)/eta);
      if ( isempty(Yinf) )
        V2 = [];
      else
        V2 = inv(Yinf);  % yikes!!!
      end
    end

    if ( isempty(V2) )
      if (verbose), fprintf("trying the anti-solution\n"), end
      [V2] = icare(-A,B,eta*(C.'*C),R,'anti');
    end
    
    if ( isempty(V2) )
      if (verbose)
        warning('approxPastEnergy: icare couldn''t find stabilizing solution')
        fprintf('approxPastEnergy: using the hamiltonian\n')
      end
      [~,V2,~] = hamiltonian(-A,B,eta*(C.'*C),R,true);
      V2 = real(V2);
    end
     
    if ( isempty(V2) )
      error('Could not find a solution to the ARE, adjust the eta parameter')
    end
  
  else % eta==0
    %  This case is described in Section II.B of the paper and requires a
    %  matrix inverse to calculate E_c.
    [V2] = lyap(A,(B*B.'));
    V2 = inv(V2); % yikes!!!!!!!!

    %  To do: look at approximating this by [V2] = icare(-A,-B,eta*(C.'*C),R)
    %  with a small value of eta (and perhaps other choices for C.'*C)
    
  end
  
  %  Reshape the resulting quadratic coefficients
  v2 = V2(:);
  v{2} = v2;

  %% k=3 case
  if ( d>2 )
    V2BB = V2*(B*B.');
  
    b = -LyapProduct(N.',v2,2);
    [v{3}] = KroneckerSumSolver(Acell(1:3),b,2,3*V2BB);

    [v{3}] = kronMonomialSymmetrize(v{3},n,3);
  end
  
  %% k>3 case (up to d)
  BVk = cell(1,d-1);
  if ( d>3 )

    for k=4:d
      BVk{k-1} = B.'*reshape(v{k-1},n,n^(k-2));

      b = -LyapProduct(N.',v{k-1},k-1);

      for i=3:(k+1)/2
        j   = k+2-i;
        tmp = BVk{i}.'*BVk{j};
        b   = b - 0.25*i*j*(vec(tmp) + vec(tmp.'));
      end

      if (mod(k,2)==0) % k is even
        i   = (k+2)/2;
        j   = i;
        tmp = BVk{i}.'*BVk{j};
        b   = b - 0.25*i*j*vec(tmp);
      end

      [v{k}] = KroneckerSumSolver(Acell(1:k),b,k-1,k*V2BB);

      [v{k}] = kronMonomialSymmetrize(v{k},n,k);
    end
  end
  
end
