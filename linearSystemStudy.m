% Scaling study to find the real-world performance of linear system solves
% This helps explain why we perform much better than the O(n^3) estimate
% in our complexity analysis.

addpath('~/KroneckerTools/src')

m=4;
p=4;
epsilon=0.001;
alpha = 0.0;

store = zeros(8,2);
kount = 0;
for n=[8,16,32,64,128,256,512,1024]
  [A,B,C,N,zInit] = getSystem3(n,m,p,epsilon,alpha);

  b = rand(n,1);  V = icare(A,B,C.'*C,eye(m));
  sysMat = A + A(1,1)*eye(n) + V*(B*B.');

  tic
  for k=1:(2^20)/n
    x = sysMat\b;
  end
  time = toc/( (2^20)/n );

  fprintf(' %04d  %12.8e\n',n,time)
  kount = kount + 1;
  store(kount,:) = [n, time];
end

% Time increases for doubling the system sizes
fprintf('The following table shows the time increase between successive \n')
fprintf('powers of 2 in performing linear system solves.  For\n')
fprintf('linear scaling, these should be 2, quadratic, 4, and the usual\n')
fprintf('complexity bound of cubic it would be 8\n')
disp(store(2:end,2)./store(1:end-1,2))
