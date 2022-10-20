%  A script to test the input normal transformation.  Here we use the
%  order 2 example in Example2 and verify that the input-normal conditions
%  are satisfied for a range of z values.
%
%  Part of the NLbalancing repository.
%
%%

setKroneckerToolsPath
addpath('..')
addpath('../examples')
addpath('../utils')

testcase=1;
switch testcase
  case 1    % a random quadratic system
    n = 2;
    m = 2;
    p = 2;
    A = rand(n,n);
    N = rand(n,n^2);
    B = rand(n,m);
    C = rand(p,n);
    degree = 6;
 
    eta = 0.1;    % values should be between -\infty and 1.
    [w] = approxFutureEnergy(A,N,B,C,eta,degree);
    [v] = approxPastEnergy(A,N,B,C,eta,degree);

  case 2
    [v,w] = runExample2(6,false,false,5);

  otherwise
    error('testInputNormal: unknown testcase\n')
end

[sigma,T] = inputNormalTransformation(v,w,5);

%  Testing that 0.5*Eminus(T{1}*z+T{2}*kron(z,z)+...) = 0.5*z.'*z
N = 101;
zmin=-.2; zmax=.2; dz=(zmax-zmin)/N;
maxErr = 0;
degree = 1;
vT = cell(size(v));
for i=2:length(v)
  vT{i} = v{i}.';
end

Eplot = zeros(N,N);
for i=1:N
  for j=1:N
    z = [ zmin+(i-1)*dz; zmin+(j-1)*dz ];
    x = kronPolyEval(T,z,degree);
    LHS = kronPolyEval(vT,x,degree+1);
    RHS = z.'*z;
    error = abs(LHS-RHS);
    Eplot(i,j) = error;
    maxErr = max(maxErr,error);
  end
end

fprintf('testInputNormal: The maximum error is %g\n',maxErr);
% figure(3)
% surf(Eplot)

%  Plot the past energy function in a neighborhood of the origin
figure(4)
Eplot = zeros(N,N);
for i=1:N
  for j=1:N
    z = [ zmin+(i-1)*dz; zmin+(j-1)*dz ];
    Eplot(i,j) = kronPolyEval(vT,z,degree+1);
  end
end
surf(Eplot)
title('Past energy function in original coordinates')

figure(5)
Eplot = zeros(N,N);
for i=1:N
  for j=1:N
    z = [ zmin+(i-1)*dz; zmin+(j-1)*dz ];
    x = kronPolyEval(T,z,degree);
    Eplot(i,j) = kronPolyEval(vT,x,degree+1);
  end
end
surf(Eplot)
title('Past energy function in transformed coordinates')
