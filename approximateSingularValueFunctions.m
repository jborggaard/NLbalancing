function [c] = approximateSingularValueFunctions(T,w,sigma,degree)
%Computes the polynomial coefficients of the singular value functions.
%
%  c = approximateSingularValueFunctions(T,w,sigma,degree)
%
%  Computes the coefficients c(i) of the polynomial approximation to the
%  singular value functions
%
%  xi(z) = sigma + c{1}*z + c{2}*z.^2 + ... + c{degree}*z.^degree, 
%
%  for i=1,...,n, where T is the transformation matrices cell array
%  and the Hankel singular values at zero (the constant terms in the
%  singular value functions) are computed as
%
%  [T,sigma] = computeBalancingTransformation(v,w,degree+1)
%
%  and w are the coefficients in the expansion of the past energy
%  function.
%
%  The coefficients are approximated out to _degree_ terms.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Reference:  Nonlinear balanced truncation: Part 2--Model reduction on
%              manifolds, by Kramer, Gugercin, and Borggaard, arXiv
%
%              See Algorithm 2.
%
%  Part of the NLbalancing repository.
%%

  c  = cell(1,degree);

  % Create a vec function for readability
  vec = @(X) X(:);

  Sigma = diag(sigma);
  
  T1 = T{1};
  T2 = T{2};
  
  n  = size(T1,1);
  W2 = reshape(w{2},n,n);
  
  tmp = T2.'*W2*T1;
  rhs = vec(tmp) + 0.5*kroneckerRight(w{3}.',T1).';
  
  Idx = indexSet(n,1);

  c{1} = Sigma\rhs(Idx);
  

  for k=2:degree
    C = vec(T{k+1}.'*W2*T{1}) + vec(T{1}.'*W2*T{k+1});
      
    for i=2:k
      j = k+2-i;
      C = C + vec(T{i}.'*W2*T{j});
    end

    for i=3:k+1
      C = C + calTTv(T,i,k+2,w{i});
    end
 
    Idx = indexSet(n,2);
    rhs = C(Idx);
    for i=1:k-1
      rhs = rhs- c{i}.*c{k-i};
    end

    c{k} = -0.5*Sigma\rhs;
  end

end

function idx = indexSet(n,d)
  % Produces the list of indices 
  idx = zeros(n,1);

  sum = 0;
  for l=1:d+1
    sum = sum + n^l;
  end

  for i=1:n
    idx(i) = i + (i-1)*sum;
  end
end
