function indexSet = mkIndices(m,k)
%  Returns all combinations of m natural numbers that sum to k.
%
%  indexSet = mkIndices(m,k); 
%
%  Used to compute calT_{m,k} in the function calTtimesv
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Reference:  Nonlinear balanced truncation model reduction for large-scale
%              polynomial systems, arXiv
%
%              See Algorithm 3.
%
%  Part of the NLbalancing repository.
%%
  indexSet = [];

  if m==1
    indexSet = k;

  else % fill the remaining indices through recursion

    for i=1:k-m+1
      J = mkIndices(m-1,k-i);
      rJ = size(J,1);

      indexSet = [indexSet; i*ones(rJ,1) J];
    end

  end

end