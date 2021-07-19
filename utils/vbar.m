function VBar = vbar(v,x)
%VBar Compute the state-dependent quadratic energy matrix from coefficients
%     of the cell array, v and the state, x.
%
%   Given an energy function of the form
%     v(x) = v2.'*kron(x,x) + ...
%            v3.'*kron(kron(x,x),x) + ...
%            v4.'*kron(kron(kron(x,x),x),x) + ...     
% 
%   the coefficients of v(x) are given in a cell array where
%     v{2} = v2, v{3} = v3, etc.
%
%   This is to be written as
%     v(x) = vBar(x)*kron(x,x) = x.'*VBar*x.
%
%   This function computes the matrix VBar.
%
%   Usage:  [VBar] = vbar(v,x)
%
%   where 
%         v = {v2,v3,v4}
%   and
%         x   is the value of the state to construct this approximation at.
%
%   Author: Jeff Borggaard, Virginia Tech
%
%   Part of the NLbalancing library.
%
%%    
  n   = length(x);
  deg = length(v);

  vBar = v{2}.';

  if (deg>2)
    xd   = x;
    vBar = vBar + xd.'*reshape(v{3},n,n^2);

    for d=4:deg
      xd = kron(xd,x);
      vBar = vBar + xd.'*reshape(v{d},n^(d-2),n^2);
    end
  end

  VBar = reshape(vBar,n,n);
end
