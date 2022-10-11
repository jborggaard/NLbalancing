%%  Testing solution to (52)  (V3 equation)

  addpath('..')
  setKroneckerToolsPath
  
  A = rand(2,2);
  B = rand(2,1);
  C = rand(1,2);
  N = rand(2,4);   N(:,2:3) = [ sum(N(:,2:3),2)/2 sum(N(:,2:3),2)/2 ];

  eta = 1;


  n = size(A,1);
  I = eye(n);

%%  find the V2 solution using Lyap or ARE to get started
  if ( eta == 0 )
      V2 = lyap(A',(C'*C));
  else
      V2 = are(A,eta*(B*B'),(C'*C));
  end

%  Find the V3 solution
  v2 = V2(:);
  rhs = -LyapProduct(N.',v2,2);
  if ( eta == 0 )  
      mat = kron(A,kron(I,I))+kron(I,kron(A,I))+kron(I,kron(I,A));

      V3_1 = mat\rhs;  % direct solution
      
      V3 = KroneckerSumSolver({A,A,A},rhs,2);
      
  else
      mat = kron(A,kron(I,I))+kron(I,kron(A,I))+kron(I,kron(I,A))-3*eta*kron(kron(I,I),V2'*(B*B'));
      V3_1 = mat\rhs;  % direct solution
      
      V3 = KroneckerSumSolver({A,A,A},rhs,2,-3*eta*V2'*(B*B'));
      
  end
  
  fprintf('testKroneckerSumSolver: The error between the kronecker sum solver\n')
  fprintf('and a direct solver is %g\n\n',norm(V3_1-V3,'inf') )
