function [w] = approxFutureEnergy(A, N, g, C, eta, d, verbose)
%  Calculates a polynomial approximation to the future energy function
%  for a quadratic drift, polynomial input system. The default usage is
%
%  w = approxFutureEnergy(A,N,g,C,eta,d,verbose)
%
%  where 'verbose' is an optional argument. If the system has a constant
%  input vector field Bu, the matrix B may be passes in place of a cell
%  array 'g'.The cell array 'g' should be g{1} = B = G_0, g{2} = G1, g{3} = G2 ...
%
%  Computes a degree d polynomial approximation to the future energy function
%
%          E^+(x) = 1/2 ( w{2}'*kron(x,x) + ... + w{d}'*kron(.. x) )
%
%  for the polynomial system
%
%    \dot{x} = Ax + Bu + N*kron(x,x) + G1*kron(x,u) + G2*kron(x,x,u) + ...
%          y = Cx
%
%  where eta = 1-gamma^(-2), gamma is the parameter in the algebraic Riccati
%  equation
%
%    A'*W2 + W2*A - eta*W2*B*B'*W2 + C'*C = 0,
%
%  and in the subsequent linear systems arising from the Future H-infinity
%  Hamilton-Jacobi-Bellman Partial Differential Equation.
%
%  Note that w{2} = vec(W2) = W2(:).  Details are in Section III.B of reference [1].
%
%  Requires the following functions from the KroneckerTools repository:
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
if ismatrix(g)
  % Will reduce to Jeff's original code
  B = g;
  l = 0;
  g = {B};
else
  % QB or polynomial input balancing
  B = g{1};
  l = length(g) - 1;
end

n = size(A, 1); % A should be n-by-n
m = size(B, 2); % B should be n-by-m
p = size(C, 1); % C should be p-by-n

% Create a vec function for readability
vec = @(X) X(:);

%% k=2 case
R = eye(m) / eta;

if (eta > 0)
  [W2] = icare(A, B, (C.' * C), R);

  if (isempty(W2) && verbose)
    warning('approxFutureEnergy: icare couldn''t find stabilizing solution')
  end

elseif (eta < 0)
  [W2] = icare(A, B, (C.' * C), R, 'anti');

  if (isempty(W2) && verbose)
    warning('approxFutureEnergy: icare couldn''t find stabilizing solution')
    warning('approxFutureEnergy: using the hamiltonian')
    [~, W2, ~] = hamiltonian(A, B, C.' * C, R, true);
  end

else % eta==0
  [W2] = lyap(A.', (C.' * C));

end

if (isempty(W2))
  error('approxFutureEnergy: Can''t find a stabilizing solution')
end

%  Check the residual of the Riccati/Lyapunov equation
if (verbose)
  RES = A' * W2 + W2 * A - eta * (W2 * B) * (B' * W2) + C' * C;
  fprintf('The residual of the Riccati equation is %g\n', norm(RES, 'inf'));
end

%  Reshape the resulting quadratic coefficients
w{2} = vec(W2);

%% k=3 case
if (d > 2)
  GaWb = cell(2 * l + 1, d - 1); % Pre-compute G_a.'*W_b, etc for all the a,b we need
  GaWb{1, 2} = B.' * W2;
  % set up the generalized Lyapunov solver
  [Acell{1:d}] = deal(A.' - eta * W2 * (B * B.'));

  b = -LyapProduct(N.', w{2}, 2);

  if l > 0 % New for QB/polynomial input
    Im = speye(m);
    GaWb{2, 2} = g{2}.' * W2;
    b = b + 2 * eta * kron(speye(n ^ 3), vec(Im).') * vec(kron(GaWb{2, 2}, GaWb{1, 2}));
  end

  [w{3}] = KroneckerSumSolver(Acell(1:3), b, 3); % Solve Ax=b for k=3
  [w{3}] = kronMonomialSymmetrize(w{3}, n, 3); % Symmetrize w3

  %% k>3 cases (up to d)
  for k = 4:d
    GaWb{1, k - 1} = B.' * reshape(w{k - 1}, n, n ^ (k - 2));

    b = -LyapProduct(N.', w{k - 1}, k - 1); % Pre-compue all the L(N') terms

    % Now add all the terms from the 'B' sum by looping through the i and j
    for i = 3:(k + 1) / 2 % i+j=k+2
      j = k + 2 - i;
      tmp = GaWb{1, i}.' * GaWb{1, j};
      b = b + 0.25 * eta * i * j * (vec(tmp) + vec(tmp.'));
    end

    if ~mod(k, 2) % k is even
      i = (k + 2) / 2;
      j = i;
      tmp = GaWb{1, i}.' * GaWb{1, j};
      b = b + 0.25 * eta * i * j * vec(tmp);
    end

    % Now add the higher order polynomial terms by iterating through the sums
    [g{l + 2:2 * l + 1}] = deal(0); % Need an extra space in g because of NaVb indexing

    for o = 1:2 * l
      GaWb{o + 1, k - 1} = g{o + 1}.' * reshape(w{k - 1}, n, n ^ (k - 2));

      for p = max(0, o - l):min(o, l)

        for i = 2:k - o
          q = o - p;
          j = k - o - i + 2;
          tmp = kron(speye(n ^ p), vec(Im).') ...
            * kron(vec(GaWb{q + 1, j}).', kron(GaWb{p + 1, i}, Im)) ...
            * kron(speye(n ^ (j - 1)), kron(perfectShuffle(n ^ (i - 1), n ^ q * m), Im)) ...
            * kron(speye(n ^ (k - p)), vec(Im));
          b = b + 0.25 * eta * i * j * vec(tmp);
        end

      end

    end

    [w{k}] = KroneckerSumSolver(Acell(1:k), b, k);
    [w{k}] = kronMonomialSymmetrize(w{k}, n, k);
  end

end

end
