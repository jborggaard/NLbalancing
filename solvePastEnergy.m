function [w] = solvePastEnergy(A,N,B,C,eta,d)
%  Calculates a polynomial approximation to the past energy function
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
    [W2] = icare(-A,B,eta*(C.'*C),R);  % (-A,B) is controllable if (A,B) is
    if ( isempty(W2) )
      warning('solvePastEnergy: icare couldn''t find stabilizing solution')
    end
    
    if ( isempty(W2) )
      fprintf("trying the anti-solution")
      [W2] = icare(-A,B,eta*(C.'*C),R,'anti');
    end
    
    if ( isempty(W2) )
      warning('solvePastEnergy: icare couldn''t find stabilizing solution')
      fprintf('solvePastEnergy: using the hamiltonian')
      [~,W2,E] = hamiltonian(-A,B,eta*(C.'*C),R,true);
    end
     
    if ( isempty(W2) )
      error('Could not find a solution to the ARE')
    end
    %[~,W2,~] = hamiltonian(-A,B,eta*(C.'*C),R);
    %W2 = real(W2);
  
  else % eta==0
    warning('read the paper')
    [W2] = lyap(A,(B*B.'));
    W2 = inv(W2); % yikes!!!!!!!!
    
  end
  
  w2 = W2(:);
  w{2} = w2;

  if ( d>2 )
    W2BB = W2*(B*B.');
  
    b = -LyapProduct(N.',w2,2);
    [w3] = KroneckerSumSolver(Acell(1:3),b,2,3*W2BB);

    S  = Kron2CT(n,3);
    C  = CT2Kron(n,3);
    w3 = C*S*w3;
    
    w{3} = w3;
  end
  
  if ( d>3 )
    W3 = reshape(w3,n,n^2);
    W3BBW3 = (W3.'*B)*(B.'*W3);  % just for now...
    b = -LyapProduct(N.',w3,3) - 9*W3BBW3(:)/4;
    [w4] = KroneckerSumSolver(Acell(1:4),b,3,4*W2BB);

    S  = Kron2CT(n,4);
    C  = CT2Kron(n,4);
    w4 = C*S*w4;
    
    w{4} = w4;
  end
  
  if ( d>4 )
    W4 = reshape(w4,n,n^3);
    W3BBW4 = (W3.'*B)*(B.'*W4);
    W4BBW3 = W3BBW4.';
    b = -LyapProduct(N.',w4,4) - 12*W3BBW4(:)/4 - ...
                                 12*W4BBW3(:)/4;
    [w5] = KroneckerSumSolver(Acell(1:5),b,4,5*W2BB);

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
    b = -LyapProduct(N.',w5,5) - 15*W3BBW5(:)/4 - ...
                                 16*W4BBW4(:)/4 - ...
                                 15*W5BBW3(:)/4;
    [w6] = KroneckerSumSolver(Acell(1:6),b,5,6*W2BB);

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
    b = -LyapProduct(N.',w6,6) - 18*W3BBW6(:)/4 - ...
                                 20*W4BBW5(:)/4 - ...
                                 20*W5BBW4(:)/4 - ...
                                 18*W6BBW3(:)/4;
    [w7] = KroneckerSumSolver(Acell(1:7),b,6,7*W2BB);

    S  = Kron2CT(n,7);
    C  = CT2Kron(n,7);
    w7 = C*S*w7;
    
    w{7} = w7;
  end
  
  if ( d>7 )
    W7 = reshape(w7,n,n^5);
    W3BBW7 = (W3.'*B)*(B.'*W7);
    W4BBW6 = (W4.'*B)*(B.'*W6);
    W5BBW5 = (W5.'*B)*(B.'*W5);
    W6BBW4 = W4BBW6.';
    W7BBW3 = W3BBW7.';
    b = -LyapProduct(N.',w7,7) - 21*W3BBW7(:)/4 - ...
                                 24*W4BBW6(:)/4 - ...
                                 25*W5BBW5(:)/4 - ...
                                 24*W6BBW4(:)/4 - ...
                                 21*W7BBW3(:)/4;
    [w8] = KroneckerSumSolver(Acell(1:8),b,7,8*W2BB);
    w{8} = w8;
  end
  
  if ( d>8 )
    warning('solvePastEnergy: monomial terms higher than degree 8 not computed')
  end
end
