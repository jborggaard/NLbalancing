function [v] = approxPastEnergy(A, N, g, C, eta, d, verbose)
%  Calculates a polynomial approximation to the past energy function
%  for a quadratic drift, polynomial input system. The default usage is
%
%  v = approxPastEnergy(A,N,g,C,eta,d,verbose)
%
%  where 'verbose' is an optional argument. If the system has a constant
%  input vector field Bu, the matrix B may be passes in place of the cell
%  array 'g'. The cell array 'g' should be g{1} = B = G_0, g{2} = G1, g{3} = G2 ...
%
%  Computes a degree d polynomial approximation to the past energy function
%
%          E^-(x) = 1/2 ( v{2}'*kron(x,x) + ... + v{d}'*kron(.. x) )
%
%  for the polynomial system
%
%    \dot{x} = Ax + Bu + N*kron(x,x) + G1*kron(x,u) + G2*kron(x,x,u) + ...
%          y = Cx
%
%  where eta = 1-gamma^(-2), gamma is the parameter in the algebraic Riccati
%  equation
%
%    A'*V2 + V2*A + V2*B*B'*V2 - eta*C'*C = 0.
%
%  and in the subsequent linear systems arising from the Past H-infinity
%  Hamilton-Jacobi-Bellman Partial Differential Equation.
%
%  Note that v{2} = vec(V2) = V2(:).  Details are in Section III.C of reference [1].
%
%  Requires the following functions from the KroneckerTools repository
%      KroneckerSumSolver
%      kronMonomialSymmetrize
%      LyapProduct
%
%  Authors: Jeff Borggaard, Virginia Tech
%           Nick Corbin, UCSD
%
%  License: MIT
%
%  Reference: [1] Nonlinear balanced truncation: Part 1--Computing energy
%             functions, by Kramer, Gugercin, and Borggaard, arXiv:2209.07645.
%             [2] Nonlinear balanced truncation for quadratic bilinear
%             systems, by Nick Corbin ... (in progress).
%
%             See Algorithm 1 in [1].
%
%  Part of the NLbalancing repository.
%%

if (nargin < 7)
  verbose = false;
end

% Create pointer/shortcuts for dynamical system polynomial coefficients
if iscell(g)
  % QB or polynomial input balancing
  B = g{1};
  l = length(g) - 1;
else
  % Will reduce to Jeff's original code
  B = g;
  l = 0;
  g = {B};
end

n = size(A, 1); % A should be n-by-n
m = size(B, 2); % B should be n-by-m
p = size(C, 1); % C should be p-by-n

% Create a vec function for readability
vec = @(X) X(:);

%% k=2 case
R = eye(m);

if (eta ~= 0)
  %  We multiply the ARE by -1 to put it into the standard form in icare,
  %  and know (-A,-B) is a controllable pair if (A,B) is.
  [V2] = icare(-A, -B, eta * (C.' * C), R);

  if (isempty(V2) && verbose)
    warning('approxPastEnergy: icare couldn''t find stabilizing solution')
  end

  if (isempty(V2))

    if (verbose)
      warning('approxPastEnergy: using matrix inverse')
    end

    [Yinf] = icare(A.', C.', B * B.', eye(p) / eta);

    if (isempty(Yinf))
      V2 = [];
    else
      V2 = inv(Yinf); % yikes!!!
    end

  end

  if (isempty(V2))
    if (verbose), fprintf("trying the anti-solution\n"), end
    [V2] = icare(-A, B, eta * (C.' * C), R, 'anti');
  end

  if (isempty(V2))

    if (verbose)
      warning('approxPastEnergy: icare couldn''t find stabilizing solution')
      fprintf('approxPastEnergy: using the hamiltonian\n')
    end

    [~, V2, ~] = hamiltonian(-A, B, eta * (C.' * C), R, true);
    V2 = real(V2);
  end

  if (isempty(V2))
    error('Could not find a solution to the ARE, adjust the eta parameter')
  end

else % eta==0
  %  This case is described in Section II.B of the paper and requires a
  %  matrix inverse to calculate E_c.
  [V2] = lyap(A, (B * B.'));
  V2 = inv(V2); % yikes!!!!!!!!

  %  To do: look at approximating this by [V2] = icare(-A,-B,eta*(C.'*C),R)
  %  with a small value of eta (and perhaps other choices for C.'*C)

end

%  Reshape the resulting quadratic coefficients
v{2} = vec(V2);

%% k=3 case
if (d > 2)
  GaVb = cell(2 * l + 1, d - 1); % Pre-compute G_a.'*V_b, etc for all the a,b we need
  GaVb{1, 2} = B.' * V2; % Recall g{1} = B
  % set up the generalized Lyapunov solver
  [Acell{1:d}] = deal(A.' + V2 * (B * B.'));

  b = -LyapProduct(N.', v{2}, 2);

  if l > 0 % New for QB/polynomial input
    Im = speye(m);
    GaVb{2, 2} = g{2}.' * V2;
    b = b - 2 * kron(speye(n ^ 3), vec(Im).') * vec(kron(GaVb{2, 2}, GaVb{1, 2}));
  end

  [v{3}] = KroneckerSumSolver(Acell(1:3), b, 3);
  [v{3}] = kronMonomialSymmetrize(v{3}, n, 3);

  %% k>3 case (up to d)
  for k = 4:d
    GaVb{1, k - 1} = B.' * reshape(v{k - 1}, n, n ^ (k - 2));

    b = -LyapProduct(N.', v{k - 1}, k - 1);

    for i = 3:(k + 1) / 2
      j = k + 2 - i;
      tmp = GaVb{1, i}.' * GaVb{1, j};
      b = b - 0.25 * i * j * (vec(tmp) + vec(tmp.'));
    end

    if ~mod(k, 2) % k is even
      i = (k + 2) / 2;
      j = i;
      tmp = GaVb{1, i}.' * GaVb{1, j};
      b = b - 0.25 * i * j * vec(tmp);
    end

    % Now add the higher order polynomial terms by iterating through the sums
    [g{l + 2:2 * l + 1}] = deal(0); % Need an extra space in g because of GaVb indexing

    for o = 1:2 * l
      GaVb{o + 1, k - 1} = g{o + 1}.' * reshape(v{k - 1}, n, n ^ (k - 2));

      for p = max(0, o - l):min(o, l)

        for i = 2:k - o
          q = o - p;
          j = k - o - i + 2;
          tmp = kron(speye(n ^ p), vec(Im).') ...
            * kron(vec(GaVb{q + 1, j}).', kron(GaVb{p + 1, i}, Im)) ...
            * kron(speye(n ^ (j - 1)), kron(perfectShuffle(n ^ (i - 1), n ^ q * m), Im)) ...
            * kron(speye(n ^ (k - p)), vec(Im));
          b = b - 0.25 * i * j * vec(tmp);
        end

      end

    end

    [v{k}] = KroneckerSumSolver(Acell(1:k), b, k);
    [v{k}] = kronMonomialSymmetrize(v{k}, n, k);
  end

end

end
