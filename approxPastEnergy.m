function [v] = approxPastEnergy(A,N,B,C,eta,d)
%  Calculates a polynomial approximation to the past energy function
%  for a quadratic system.
%
%  v = approxPastEnergy(A,N,B,C,eta,d) 
%
%  Computes a degree d polynomial approximation to the past energy function 
%
%          E^-(x) = v{2}*kron(x,x) + ... + v{d}*kron(.. x)
%
%  for the polynomial system
%
%    \dot{x} = Ax + Bu + N*kron(x,x)
%          y = Cx
%
%  where eta = 1-gamma^(-2), gamma is the parameter in the algebraic Riccati
%  equation
%
%    A'*V2 + V2*A + V2*B*B'*V2 - eta*C'*C = 0.
%
%  Note that v{2} = vec(V2) = V2(:).  Details are in Section 3.2 of the paper.
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

  R = eye(m);
  
  Acell = cell(1,d);
  for i=1:d
    Acell{i} = A.';
  end

  if ( eta~=0 )
    %  We multiply the ARE by -1 to put it into the standard form,
    %  and know (-A,B) is a controllable pair if (A,B) is.
    [V2] = icare(-A,B,eta*(C.'*C),R);
    if ( isempty(V2) )
      warning('approxPastEnergy: icare couldn''t find stabilizing solution')
    end
    
    if ( isempty(V2) )
      fprintf("trying the anti-solution")
      [V2] = icare(-A,B,eta*(C.'*C),R,'anti');
    end
    
    if ( isempty(V2) )
      warning('approxPastEnergy: icare couldn''t find stabilizing solution')
      fprintf('approxPastEnergy: using the hamiltonian')
      [~,V2,E] = hamiltonian(-A,B,eta*(C.'*C),R,true);
    end
     
    if ( isempty(V2) )
      error('Could not find a solution to the ARE, adjust the eta parameter')
    end
    %[~,V2,~] = hamiltonian(-A,B,eta*(C.'*C),R);
    %V2 = real(V2);
  
  else % eta==0
    warning('read the paper')
    [V2] = lyap(A,(B*B.'));
    V2 = inv(V2); % yikes!!!!!!!!
    
  end
  
  v2 = V2(:);
  v{2} = v2;

  if ( d>2 )
    V2BB = V2*(B*B.');
  
    b = -LyapProduct(N.',v2,2);
    [v3] = KroneckerSumSolver(Acell(1:3),b,2,3*V2BB);

    S  = Kron2CT(n,3);
    C  = CT2Kron(n,3);
    v3 = C*S*v3;
    
    v{3} = v3;
  end
  
  if ( d>3 )
    V3 = reshape(v3,n,n^2);
    V3BBV3 = (V3.'*B)*(B.'*V3);  % just for now...
    b = -LyapProduct(N.',v3,3) - 9*V3BBV3(:)/4;
    [v4] = KroneckerSumSolver(Acell(1:4),b,3,4*V2BB);

    S  = Kron2CT(n,4);
    C  = CT2Kron(n,4);
    v4 = C*S*v4;
    
    v{4} = v4;
  end
  
  if ( d>4 )
    V4 = reshape(v4,n,n^3);
    V3BBV4 = (V3.'*B)*(B.'*V4);
    V4BBV3 = V3BBV4.';
    b = -LyapProduct(N.',v4,4) - 12*V3BBV4(:)/4 - ...
                                 12*V4BBV3(:)/4;
    [v5] = KroneckerSumSolver(Acell(1:5),b,4,5*V2BB);

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
    b = -LyapProduct(N.',v5,5) - 15*V3BBV5(:)/4 - ...
                                 16*V4BBV4(:)/4 - ...
                                 15*V5BBV3(:)/4;
    [v6] = KroneckerSumSolver(Acell(1:6),b,5,6*V2BB);

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
    b = -LyapProduct(N.',v6,6) - 18*V3BBV6(:)/4 - ...
                                 20*V4BBV5(:)/4 - ...
                                 20*V5BBV4(:)/4 - ...
                                 18*V6BBV3(:)/4;
    [v7] = KroneckerSumSolver(Acell(1:7),b,6,7*V2BB);

    S  = Kron2CT(n,7);
    C  = CT2Kron(n,7);
    v7 = C*S*v7;
    
    v{7} = v7;
  end
  
  if ( d>7 )
    V7 = reshape(v7,n,n^5);
    V3BBV7 = (V3.'*B)*(B.'*V7);
    V4BBV6 = (V4.'*B)*(B.'*V6);
    V5BBV5 = (V5.'*B)*(B.'*V5);
    V6BBV4 = V4BBV6.';
    V7BBV3 = V3BBV7.';
    b = -LyapProduct(N.',v7,7) - 21*V3BBV7(:)/4 - ...
                                 24*V4BBV6(:)/4 - ...
                                 25*V5BBV5(:)/4 - ...
                                 24*V6BBV4(:)/4 - ...
                                 21*V7BBV3(:)/4;
    [v8] = KroneckerSumSolver(Acell(1:8),b,7,8*V2BB);
    v{8} = v8;
  end
  
  if ( d>8 )
    warning('approxPastEnergy: monomial terms higher than degree 8 not computed')
  end
end
