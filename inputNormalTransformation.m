function [sigma,T] = inputNormalTransformation(v,w,degree)
%  Computes a polynomial approximation to the input-normal balancing 
%  transformation for a system with polynomial nonlinearities.
%
%     [sigma,T] = approximateBalancingTransformation(v,w,degree)
%
%  using polynomial approximations to the past and future energy functions.
%  v and w contain the coefficients to the past and future energy
%  functions, respectively.  The balancing transformation,
%
%      x = T{1}z + T{2}kron(z,z) + ... + T{degree}kron(kron...,z),z)
%
%  Terms out to v{degree+1} and w{degree+1} must be defined in the input
%
%  The input-normal/output-diagonal balancing transformation is described in
%  Section 4.1 of the reference.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Reference:  Nonlinear balanced truncation model reduction for large-scale
%              polynomial systems, arXiv
%
%              See Algorithm 2.
%
%  Part of the NLbalancing repository.
%%
  fprintf('Computing the %d degree input-normal transformation\n',degree)
  
  validateattributes(v,{'cell'},{})
  validateattributes(w,{'cell'},{})
  dv = length(v);
  dw = length(w);

  if (degree<1)
    error('computeBalancing: degree must be at least 1')
  end
  
  if (dv<degree+1 || dw<degree+1)
    error('computeBalancing: we need degree %d terms in the energy function',...
          degree+1)
  end
  
  n  = sqrt(length(v{2}));
  V2 = reshape(v{2},n,n);
  W2 = reshape(w{2},n,n);
  R  = chol(V2,'upper');
  L  = chol(W2,'upper');
  
  [~,Xi,V] = svd(L.'/R.'); % from Theorem 9
  sigma = diag(Xi);
  
%  T{1} = R\(V/sqrt(Xi));  % T1inv = sqrt(Xi)\(U.'*L.');
  T{1} = R.'\V;  % T1inv = sqrt(Xi)\(U.'*L.');
  
  if (degree>1)
    vT3  = kroneckerRight(v{3}.',T{1});
    T{2} = -0.5*T{1}*reshape(vT3,n^2,n).';
  end
  
  if (degree>2)
    vT4  = kroneckerRight(v{4}.',T{1});
    vL3  = LyapProduct(T{2}.',v{3},3).';
    tmp  = T{2}.'*V2*T{2};  tmp = tmp(:).';
    term = tmp + vL3 + vT4;
    T{3} = -0.5*T{1}*reshape(term,n^3,n).';
  end
  
  if (degree>3)
    vT5  = kroneckerRight(v{5}.',T{1});
    vL4  = LyapProduct(T{2}.',v{4},4).';
    vL3  = LyapProduct(T{3}.',v{3},3).';
    tmp0 = T{3}.'*V2*T{2};   tmp1 = tmp0.';
    tmp0 = tmp0(:).';        tmp1 = tmp1(:).';
    term = tmp0 + tmp1 + vL3 + vL4 + vT5;
    T{4} = -0.5*T{1}*reshape(term,n^4,n).';
  end
  
  if (degree>4)
    vT6  = kroneckerRight(v{6}.',T{1});
    vL5  = LyapProduct(T{2}.',v{5},5).';
    vL4  = LyapProduct(T{3}.',v{4},4).';
    vL3  = LyapProduct(T{4}.',v{3},3).';
    tmp0 = T{2}.'*V2*T{4};  tmp0 = tmp0(:).';
    tmp1 = T{3}.'*V2*T{3};  tmp1 = tmp1(:).';
    tmp2 = tmp0.';          tmp2 = tmp2(:).';
    term = tmp0 + tmp1 + tmp2 + vL3 + vL4 + vL5 + vT6;
    T{5} = -0.5*T{1}*reshape(term,n^5,n).';
  end
end
