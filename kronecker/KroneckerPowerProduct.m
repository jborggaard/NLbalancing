function [Mv] = KroneckerPowerProduct(M,v,d,side)
%LyapProduct Computes the Kronecker powers of a matrix multiplied by a vector v.
%     Given a matrix M, compute the product with the column vector v.
%     Optionally, this does the left multiplication by a row vector v.
%
%   Usage:  [P] = LyapProduct(M,v,d,side);
%
%     P = (M x M x ... x M)*v    or    P = v.'*(M x M x ... x M)
%          |-- d terms --|                      |-- d terms --|            
%                                   with the optional argument side='left' 
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Reference:  Nonlinear balanced trunction model reduction for
%  large-scale polynomial systems, arXiv
%
%
%  Part of the NLbalancing library.
%%

  if (nargin==4 && side=='left')
    Mv = KroneckerPowerProduct(M.',v,d).';

    warning('need to rearrange')
  end

  [m,n] = size(M);
  t     = size(v,1);  % right now, we are assuming v is a single column
  if ( t ~= n^d )
    error('The dimensions of v do not match the degree of the multiLyapunov matrix')
  end
  
  if (d<2)
    Mv = M*v;
  
  else
    V = reshape(v,n^(d-1),n);
    Mv = reshape(V*M.',m*n^(d-1),1);
  
    V = reshape(v,n,n^(d-1));
    Mv = Mv + reshape(M*V,m*n^(d-1),1);
  
    for l=1:d-2
      V1 = reshape(v,n^(d-l),n^l);
    
      mat = zeros(m*n^(d-l-1),n^l);
      for i=1:n^l
        vi = V1(:,i);
        mat(:,i) = reshape( reshape(vi,n^(d-l-1),n)*M.', m*n^(d-l-1),1);
      end
      Mv = Mv + reshape(mat,m*n^(d-1),1);
    end
  end
  
  Mv = Mv(:);
end

