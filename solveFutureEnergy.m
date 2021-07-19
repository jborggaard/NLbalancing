function [v] = solveFutureEnergy(A,N,B,C,eta,d)
%  Calculates a polynomial approximation to the future energy function
%  for a quadratic system.
%
%  V = solveFutureEnergy(A,N,B,C,eta,d) computes a degree d polynomial approximation
%  to the future energy function  E^+(x) = V{2}*kron(x,x) + ... + V{d}*kron(.. x)
%  for the polynomial system
%
%    \dot{x} = Ax + Bu + N*kron(x,x)
%
%    eta = 1-gamma^(-2)
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Reference:  Nonlinear balanced truncation model reduction for large-scale
%              polynomial systems, arXiv
%
%
%  Part of the QQGbalancing repository.
%%
  n = size(A,1);
  m = size(B,2);

  R = eye(m)/eta;

  Acell = cell(1,d);
  for i=1:d
    Acell{i} = A.';
  end

  if ( eta>0 )
    [V2] = icare(A,B,(C.'*C),R);
    if ( isempty(V2) )
      warning('solveFutureEnergy: icare couldn''t find stabilizing solution')
    end
    
  elseif ( eta<0 )
    [V2] = icare(A,B,(C.'*C),R,'anti');
    
    if ( isempty(V2) )
      warning('solveFutureEnergy: icare couldn''t find stabilizing solution')
      warning('solveFutureEnergy: using the hamiltonian')
      [~,V2,~] = hamiltonian(A,B,C.'*C,R,true);
    end
    
  else % eta==0
    [V2] = lyap(A.',(C.'*C));
  end

  %  Reshape the resulting second-order term 
  v2 = V2(:);
  v{2} = v2;

  if ( d>2 )
    b = -LyapProduct(N.',v2,2);
    [v3] = KroneckerSumSolver(Acell(1:3),b,2,-3*eta*V2*(B*B.'));
    
    S  = Kron2CT(n,3);
    C  = CT2Kron(n,3);
    v3 = C*S*v3;

    v{3} = v3;
  end
  
  if ( d>3 )
    V3 = reshape(v3,n,n^2);
    V3BBV3 = (V3.'*B)*(B.'*V3);  % just for now...    new variables BV3, etc.
    b = -LyapProduct(N.',v3,3) + 9*eta*V3BBV3(:)/4;
    [v4] = KroneckerSumSolver(Acell(1:4),b,3,-4*eta*V2*(B*B.'));

    S  = Kron2CT(n,4);
    C  = CT2Kron(n,4);
    v4 = C*S*v4;
    
    v{4} = v4;
  end

  if ( d>4 )
    V4 = reshape(v4,n,n^3);
    V3BBV4 = (V3.'*B)*(B.'*V4);
    V4BBV3 = V3BBV4.';
    b = -LyapProduct(N.',v4,4) + 12*eta*V3BBV4(:)/4 + ...
                                 12*eta*V4BBV3(:)/4;
    [v5] = KroneckerSumSolver(Acell(1:5),b,4,-5*eta*V2*(B*B.'));

    S  = Kron2CT(n,5);
    C  = CT2Kron(n,5);
    v5 = C*S*v5;
    
    v{5} = v5;
  end

  if ( d>5 )
    V5 = reshape(v5,n,n^4);
    V3BBV5 = (V3.'*B)*(B.'*V5);
    V4BBV4 = (V4.'*B)*(B.'*V4);
    V5BBV3 = V3BBV5.';
    b = -LyapProduct(N.',v5,5) + 15*eta*V3BBV5(:)/4 + ...
                                 16*eta*V4BBV4(:)/4 + ...
                                 15*eta*V5BBV3(:)/4;
    [v6] = KroneckerSumSolver(Acell(1:6),b,5,-6*eta*V2*(B*B.'));

    S  = Kron2CT(n,6);
    C  = CT2Kron(n,6);
    v6 = C*S*v6;
    
    v{6} = v6;
  end

  if ( d>6 )
    V6 = reshape(v6,n,n^5);
    V3BBV6 = (V3.'*B)*(B.'*V6);
    V4BBV5 = (V4.'*B)*(B.'*V5);
    V5BBV4 = V4BBV5.';
    V6BBV3 = V3BBV6.';
    b = -LyapProduct(N.',v6,6) + 18*eta*V3BBV6(:)/4 + ...
                                 20*eta*V4BBV5(:)/4 + ...
                                 20*eta*V5BBV4(:)/4 + ...
                                 18*eta*V6BBV3(:)/4;
    [v7] = KroneckerSumSolver(Acell(1:7),b,6,-7*eta*V2*(B*B.'));

    S  = Kron2CT(n,7);
    C  = CT2Kron(n,7);
    v7 = C*S*v7;
    
    v{7} = v7;
  end
  
  if ( d>7 )
    V7 = reshape(v7,n,n^6);
    V3BBV7 = (V3.'*B)*(B.'*V7);
    V4BBV6 = (V4.'*B)*(B.'*V6);
    V5BBV5 = (V5.'*B)*(B.'*V5);
    V6BBV4 = V4BBV6.';
    V7BBV3 = V3BBV7.';
    b = -LyapProduct(N.',v7,7) + 21*eta*V3BBV7(:)/4 + ...
                                 24*eta*V4BBV6(:)/4 + ...
                                 25*eta*V5BBV5(:)/4 + ...
                                 24*eta*V6BBV4(:)/4 + ...
                                 21*eta*V7BBV3(:)/4;
    [v8] = KroneckerSumSolver(Acell(1:8),b,7,-8*eta*V2*(B*B.'));
    v{8} = v8;
  end
  
  if ( d>8 )
    warning('solveFutureEnergy: monomial terms past degree 8 are not computed')
  end
end
