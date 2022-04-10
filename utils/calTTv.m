function [Tv] = calTTv(T,m,k,v)
%calTTv Calculates the term \cal{T}_{m,k}.'v
%
%  Note the additional T in calTTv ensures the transpose of T is taken.
%
%         Tv = \|i\|=k kron(T_i1.',kron(T_i2.',kron(T_i3.',...,T_im.'))*v
%
%  Usage:  Tv = calTTv(T,m,k,v)
%
%  Variables:   T  a cell array of matrices.  The pth entry has size n times n^p
%               v  a vector of dimension (n^m)
%
%  This is multiplication performed recursively using Kronecker product rules.
%
%  This function assumes the functions mkIndices and kroneckerRight have 
%  been imported from the KroneckerTools repository (and in the Matlab path).
%
%  Author: Jeff Borggaard
%
%  License: MIT
%
%  Part of the NLbalancing repository.
%%
  
  % Get a list of indices
  indexSet = mkIndices(m,k);
  nTerms = size(indexSet,1);

  n  = size(T{1},2);

  Tv = zeros(1,n^k);
  for i=1:nTerms
    Tv = Tv + kroneckerRight(v.',T(indexSet(i,:)));
  end
  
  Tv = Tv.';
end

