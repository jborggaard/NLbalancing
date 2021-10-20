function [w] = approxFutureEnergy(A,N,B,C,eta,d)
%  Calculates a polynomial approximation to the future energy function
%  for a quadratic system.
%
%  w = approxFutureEnergy(A,N,B,C,eta,d) 
%
%  Computes a degree d polynomial approximation to the future energy function 
%
%          E^+(x) = w{2}'*kron(x,x) + ... + w{d}'*kron(.. x)
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
%  Note that w2 = vec(W2).  Details are in Section 3.1 of the reference.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Reference:  Nonlinear balanced truncation model reduction for large-scale
%              polynomial systems, arXiv
%
%              See Algorithm 1.
%
%  Part of the NLbalancing repository.
%%
  n = size(A,1);
  m = size(B,2);

  R = eye(m)/eta;

  Acell = cell(1,d);
  for i=1:d
    Acell{i} = A.';
  end

  if ( eta>0 )
    [W2] = icare(A,B,(C.'*C),R);
    if ( isempty(W2) )
      warning('approxFutureEnergy: icare couldn''t find stabilizing solution')
    end
    
  elseif ( eta<0 )
    [W2] = icare(A,B,(C.'*C),R,'anti');
    
    if ( isempty(W2) )
      warning('approxFutureEnergy: icare couldn''t find stabilizing solution')
      warning('approxFutureEnergy: using the hamiltonian')
      [~,W2,~] = hamiltonian(A,B,C.'*C,R,true);
    end
    
  else % eta==0
    [W2] = lyap(A.',(C.'*C));
  end

  %  Reshape the resulting second-order term 
  w2 = W2(:);
  w{2} = w2;

  if ( d>2 )
    b = -LyapProduct(N.',w2,2);
    [w3] = KroneckerSumSolver(Acell(1:3),b,2,-3*eta*W2*(B*B.'));
    
    S  = Kron2CT(n,3);
    C  = CT2Kron(n,3);
    w3 = C*S*w3;

    w{3} = w3;
  end
  
  if ( d>3 )
    W3 = reshape(w3,n,n^2);
    W3BBW3 = (W3.'*B)*(B.'*W3);  % just for now...    new variables BW3, etc.
    b = -LyapProduct(N.',w3,3) + 9*eta*W3BBW3(:)/4;
    [w4] = KroneckerSumSolver(Acell(1:4),b,3,-4*eta*W2*(B*B.'));

    S  = Kron2CT(n,4);
    C  = CT2Kron(n,4);
    w4 = C*S*w4;
    
    w{4} = w4;
  end

  if ( d>4 )
    W4 = reshape(w4,n,n^3);
    W3BBW4 = (W3.'*B)*(B.'*W4);
    W4BBW3 = W3BBW4.';
    b = -LyapProduct(N.',w4,4) + 12*eta*W3BBW4(:)/4 + ...
                                 12*eta*W4BBW3(:)/4;
    [w5] = KroneckerSumSolver(Acell(1:5),b,4,-5*eta*W2*(B*B.'));

    S  = Kron2CT(n,5);
    C  = CT2Kron(n,5);
    w5 = C*S*w5;
    
    w{5} = w5;
  end

  if ( d>5 )
    W5 = reshape(w5,n,n^4);
    W3BBW5 = (W3.'*B)*(B.'*W5);
    W4BBW4 = (W4.'*B)*(B.'*W4);
    W5BBW3 = W3BBW5.';
    b = -LyapProduct(N.',w5,5) + 15*eta*W3BBW5(:)/4 + ...
                                 16*eta*W4BBW4(:)/4 + ...
                                 15*eta*W5BBW3(:)/4;
    [w6] = KroneckerSumSolver(Acell(1:6),b,5,-6*eta*W2*(B*B.'));

    S  = Kron2CT(n,6);
    C  = CT2Kron(n,6);
    w6 = C*S*w6;
    
    w{6} = w6;
  end

  if ( d>6 )
    W6 = reshape(w6,n,n^5);
    W3BBW6 = (W3.'*B)*(B.'*W6);
    W4BBW5 = (W4.'*B)*(B.'*W5);
    W5BBW4 = W4BBW5.';
    W6BBW3 = W3BBW6.';
    b = -LyapProduct(N.',w6,6) + 18*eta*W3BBW6(:)/4 + ...
                                 20*eta*W4BBW5(:)/4 + ...
                                 20*eta*W5BBW4(:)/4 + ...
                                 18*eta*W6BBW3(:)/4;
    [w7] = KroneckerSumSolver(Acell(1:7),b,6,-7*eta*W2*(B*B.'));

    S  = Kron2CT(n,7);
    C  = CT2Kron(n,7);
    w7 = C*S*w7;
    
    w{7} = w7;
  end
  
  if ( d>7 )
    W7 = reshape(w7,n,n^6);
    W3BBW7 = (W3.'*B)*(B.'*W7);
    W4BBW6 = (W4.'*B)*(B.'*W6);
    W5BBW5 = (W5.'*B)*(B.'*W5);
    W6BBW4 = W4BBW6.';
    W7BBW3 = W3BBW7.';
    b = -LyapProduct(N.',w7,7) + 21*eta*W3BBW7(:)/4 + ...
                                 24*eta*W4BBW6(:)/4 + ...
                                 25*eta*W5BBW5(:)/4 + ...
                                 24*eta*W6BBW4(:)/4 + ...
                                 21*eta*W7BBW3(:)/4;
    [w8] = KroneckerSumSolver(Acell(1:8),b,7,-8*eta*W2*(B*B.'));
    w{8} = w8;
  end
  
  if ( d>8 )
    warning('approxFutureEnergy: monomial terms past degree 8 are not computed')
  end
end
