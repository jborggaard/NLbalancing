function [etaMin,etaMax] = estimateEta(A,B,C)
%estimateEta Estimates the best value of eta using Mustafa and Glover.
%
%   Usage:
%        [etaMin,etaMax] = estimateEta(A,B,C)
%
%   calculates gammaHat = sqrt(1+lambdaMax(Xinf*Yinf) then returns the
%   range of eta values  etaMin <= eta <= etaMax computed
%   using the estimate
%          gammaMin = gammaHat-1 <= gamma <= gammaHat+1
%   then converts this to a range over eta using
%            eta = 1 - gamma^(-2) via equation (40).
%
%  Details are in Remark 1 of the reference.
%
%  Author: Jeff Borggaard, Virginia Tech
%
%  Licence: MIT
%
%  Reference: Nonlinear balanced truncation: Part 1--Computing energy functions,
%             Kramer, Gugercin, and Borggaard, arXiv.
%
%             See Remark 1.
%
%  Part of the NLbalancing repository.
%%

  m = size(B,2);
  p = size(C,1);

  %  Compute the solution to the Hinf Riccati equation for the future energy 
  %  function of the linear problem with "gamma=inf".
  if ( issparse(A) || issparse(B) || issparse(C) )
    Xinf = icare(full(A),full(B),full(C.'*C),eye(m));
  else
    Xinf = icare(A,B,C.'*C,eye(m));
  end

  %  Compute the solution to the Hinf Riccati equation for the past energy
  %  function of the linear problem with "gamma=inf".
  if ( issparse(A) || issparse(B) || issparse(C) )
    Yinf = icare(full(A.'),full(C.'),full(B*B.'),eye(p));
  else
    Yinf = icare(A.',C.',B*B.',eye(p));
  end

  %  Now compute the estimate
  lambdaMax = eigs(Xinf*Yinf,1);
  gammaHat  = sqrt( 1 + lambdaMax );
  gammaMin  = gammaHat-1;
  gammaMax  = gammaHat+1;

  etaMin = 1-gammaMin^(-2);
  etaMax = 1-gammaMax^(-2);
end
