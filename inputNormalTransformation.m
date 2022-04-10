function [sigma,T] = inputNormalTransformation(v,w,degree)
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
%  in the z coordinates.
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
  fprintf('Computing the degree %d input-normal balancing transformation\n',degree)
  
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
  R  = chol(V2,'lower');
  L  = chol(W2,'lower');
  
  [~,Xi,V] = svd(L.'/R.'); % from Theorem 9
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
      % compute j=2 term
      temp = T{k-1}.'*V2*T{2};
      TVT = temp(:);

      for j=3:k-1
        temp = T{k-j+1}.'*V2*T{j};
        TVT = TVT + temp(:);
      end

      Tv = calTTv(T,3,k+1,v{3});
      for j=4:k+1
        temp = calTTv(T,j,k+1,v{j});
        Tv = Tv + temp; 
      end

      term = TVT + Tv;
      
%       % option 1
%       S1 = perfectShuffle(n,n^k);
%       term = (1.001*S1+eye(n^(k+1)))\term;
%       T{k} = -T{1}*reshape(term,n,n^k);

%       % option 2
%       A = inv(T{1});
%       AA = zeros(n^(k+1),n^(k+1));
%       for i=1:n^k
%         AA(n*i-1:n*i,n*i-1:n*i) = A;
%       end
%       for j=1:n
%         for i=1:n^k
%           row  = i+(j-1)*n^k;
%           cols = 1+(i-1)*n:n+(i-1)*n;
%           AA(row,cols) = AA(row,cols) + A(j,:);
%         end
%       end
%       [UU,SS,VV] = svd(AA);
%       r = rank(AA);
%       fprintf(' is rhs in span of AA?: %g\n',norm( term - UU(:,1:r)*UU(:,1:r).'*term ))
%       Tk = -AA\term;
%       T{k} = reshape(Tk,n,n^k);
      
%       % checking rhs calculation
%       if (k==3)
%         a1 = v{2}.'*kron(T{2},T{2}) + ...
%          v{3}.'*(kron(T{1},kron(T{1},T{2})) + kron(T{1},kron(T{2},T{1})) + ...
%                  kron(T{2},kron(T{1},T{1}))) + ...
%          v{4}.'*(kron(kron(kron(T{1},T{1}),T{1}),T{1}) )   ;
%         fprintf('the error in the rhs for T{3} calculation is %g\n',norm(term-a1.'))
%       end

      % option 3
      T{k} = -0.5*T{1}*reshape(term,n,n^k);

      fprintf('The residual error for T{%d} is %g\n',k,...
              norm( (kron(T{1}.',T{k}.')+kron(T{k}.',T{1}.'))*v{2} + term ) )
    end
  end


end  
