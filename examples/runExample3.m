%  This script runs two tests on the Burgers equation system.  
% 
%  The first is a computational performance benchmark with increasing system 
%  order and fixed approximation degree.
%
%  The second uses a fixed system order and tests the accuracy of the
%  energy function approximations as the degree increases.
%
%  This script generates the tables found in Section IV.C. of the paper
%    Nonlinear balanced truncation:  Part 1--Computing energy functions, 
%    by Kramer, Gugercin, and Borggaard, arXiv:2209.07645.
%
%  Part of the NLbalancing repository.
%%

eta =  0.9;

z_factor = 0.001;  % scale factor on the initial conditions.
                   % zInit = 0.5*sin(2*pi*x)^2 on (0,0.5) and 0 otherwise
                   % so z_factor = 0.001 matches the initial condition in IV.C.

m = 4;
p = 4;

epsilon = 0.001;
alpha   = 0.0;     % alpha = 1.0, not controllable at n=8

%%
%  Computational performance of the energy function approximations.
%  Since the initial times are so short, we average nTest times
%
%  This builds TABLE II
%
nTest = 10;
degree = 3;

for n=[8,16,32,64]
  fprintf('%d & ',n)
  [A,B,C,N,zInit] = getSystem3(n,m,p,epsilon,alpha);
  zInit = z_factor*zInit;

  tic; for i=1:nTest, [w] = approxFutureEnergy(full(A),N,B,C,eta,degree); end, tt=toc/nTest;
  fprintf('%10.4e & ',length(w{3}))
  fprintf('%8.2e & ',tt)

  for d=2:degree, w{d} = w{d}.'; end
  wzInit = 0.5*kronPolyEval(w,zInit,degree);
  fprintf('%12.6e \\\\ \n',wzInit)
end

%  Uncomment this block for the full table
nTest = 1;
degree = 3;

for n=[128,256,512,1024]
  fprintf('%d & ',n)

  [A,B,C,N,zInit] = getSystem3(n,m,p,epsilon,alpha);
  zInit = z_factor*zInit;

  tic; for i=1:nTest, [w] = approxFutureEnergy(full(A),N,B,C,eta,degree); end, tt=toc/nTest;
  fprintf('%10.4e & ',length(w{3}))
  fprintf('%8.2e & ',tt)
  
  for d=2:degree, w{d} = w{d}.'; end
  wzInit = 0.5*kronPolyEval(w,zInit,degree);
  fprintf('%12.6e \\\\ \n',wzInit)
end
 
fprintf('\n\n')

%%
%  Test convergence of the energy function at a "point" z_factor*zInit
%  with increasing degree...
%
%  This builds TABLE III
%

n=8;
[A,B,C,N,zInit] = getSystem3(n,m,p,epsilon,alpha);
zInit = z_factor*zInit;
for degree=[2,3,4,5,6]
  fprintf('%d & ',degree)

  [v] = approxPastEnergy(full(A),N,B,C,eta,degree,false);
  for d=2:degree, v{d} = v{d}.'; end
  vzInit = 0.5*kronPolyEval(v,zInit,degree);
  fprintf('%12.6e & ',vzInit)

  [w] = approxFutureEnergy(full(A),N,B,C,eta,degree,false);
  for d=2:degree, w{d} = w{d}.'; end
  wzInit = 0.5*kronPolyEval(w,zInit,degree);
  fprintf('%12.6e \\\\ \n',wzInit)
end

%   %% We can check the eigenvalues of A as one test that matrices are 
%   %  computed correctly.
%   g = sort(eig(A)','descend'); disp(g(1:8))
%   
%   % should converge to the following as n increases:
%   disp(-epsilon*pi^2*(1:8).^2+alpha)
