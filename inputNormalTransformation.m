function [sigma,T] = inputNormalTransformation(v,w,degree,verbose)
%  Computes a polynomial approximation to the input-normal balancing 
%  transformation for a system with polynomial nonlinearities as
%
%     [sigma,T] = inputNormalTransformation(v,w,degree)
%
%  using polynomial approximations to the past and future energy functions.
%  The variables v and w contain the coefficients to the past and future energy
%  functions, respectively.  Terms out to v{degree+1} and w{degree+1} must be 
%  defined in the input.  Thus,
% 
%      E_past(x) = 1/2 ( v{2}kron(x,x) + ... + v{degree+1}kron(kron...,x),x) )
%                = 0.5*kronPolyEval(v,x,degree+1)
% 
%  and
% 
%      E_future(x) = 1/2 ( w{2}kron(x,x) + ... + w{degree+1}kron(kron...,x),x) )
%                  = 0.5*kronPolyEval(w,x,degree+1)
% 
%  The balancing transformation then has the form
%
%      x = T{1}z + T{2}kron(z,z) + ... + T{degree}kron(kron...,z),z)
%
%  where E_past(x) = 1/2 ( z.'z ) and E_future(x) = 1/2 ( z.'diag(sigma)z )
%  in the z coordinates.  The singular value functions (sigma) are
%  approximated separately using approximateSingularValueFunction.
%
%  The input-normal/output-diagonal balancing transformation is described in
%  Section III.A of the reference.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  License: MIT
%
%  Reference:  Nonlinear balanced truncation: Part 2--Model reduction on
%              manifolds, by Kramer, Gugercin, and Borggaard, arXiv.
%
%              See Algorithm 1.
%
%  Part of the NLbalancing repository.
%%
  if (nargin<4)
    verbose = false;
  end

  if (verbose)
    fprintf('Computing the degree %d input-normal balancing transformation\n',degree)
  end

  % Create a vec function for readability
  vec = @(X) X(:);

  validateattributes(v,{'cell'},{})
  validateattributes(w,{'cell'},{})
  dv = length(v);
  dw = length(w);

  if (degree<1)
    error('inputNormalTransformation: degree must be at least 1')
  end
  
  if (dv<degree+1 || dw<degree+1)
    error('inputNormalTransformation: we need degree %d terms in the energy function',...
          degree+1)
  end

  % preallocate storage for the output T.
  T  = cell(1,degree);

  n  = sqrt(length(v{2}));
  V2 = reshape(v{2},n,n);
  W2 = reshape(w{2},n,n);
  R  = chol(V2,'lower');    % V2 = R*R.'
  L  = chol(W2,'lower');    % W2 = L*L.'
  
  [~,Xi,V] = svd(L.'/R.');  % from Theorem 2
  sigma = diag(Xi);
  
  T{1} = R.'\V;

  if (degree>1)
%    k = 2;
    vT3  = kroneckerRight(v{3}.',T{1});
%     A = inv(T{1});
%     AA = zeros(n^(k+1),n^(k+1));
%     for i=1:n^k
%       AA(n*i-1:n*i,n*i-1:n*i) = A;
%     end
%     for j=1:n
%       for i=1:n^k
%         row  = i+(j-1)*n^k;
%         cols = 1+(i-1)*n:n+(i-1)*n;
%         AA(row,cols) = AA(row,cols) + A(j,:);
%       end
%     end
%     Tk = -AA\vT3.';
%     T{k} = reshape(Tk,n,n^k);
% 
%     temp1 = T{k};
    T{2} = -0.5*T{1}*reshape(vT3,n,n^2);
%    norm(T{2}-temp1)
  end
  
  if (degree>2)
    for k=3:degree
      % compute the j=2 term
      TVT = vec(T{k-1}.'*V2*T{2});

      % compute the remaining terms in the sum of vecs
      for j=3:k-1
        TVT = TVT + vec(T{k-j+1}.'*V2*T{j});  % half of this work is redundant
      end

      % compute the j=3 term
      Tv = calTTv(T,3,k+1,v{3});

      % compute the remaining terms in the sum of calTTv
      for j=4:k+1
        temp = calTTv(T,j,k+1,v{j});
        Tv = Tv + temp; 
      end

      term = TVT + Tv;
      
      T{k} = -0.5*T{1}*reshape(term,n^k,n).';

      if (verbose)
        fprintf('The residual error for T{%d} is %g\n',k,...
                norm( (kron(T{1}.',T{k}.')+kron(T{k}.',T{1}.'))*v{2} + term ) )
      end
    end
  end


end  
