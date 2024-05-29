% Validation tests for energy function approximations.
% 
% Using our examples, we choose a point x0 in state-space, then
%
% 1.) Compute polynomial approximations to the energy functions using
%     approxFutureEnergy and approxPastEnergy for eta=0
%
% 2.) Numerically integrate the controlled system backward and uncontrolled
%     system forward and record the appropriate integrals.
%
% Note: the numerical integration goes to a relatively large number, but
% will ultimately underpredict the actual value of the energy functions
% (when neglecting numerical integration error).

setKroneckerToolsPath
addpath('..')
addpath('../examples')
addpath('../utils')

%% set up the test
testCase = 2;
testFutureEnergyFunction = true;
testPastEnergyFunction   = false;
visualize = false;

%% set up the problem
eta = 0.0;  % should be 0 if we want the integrated quantities to approx e.f.s

switch (testCase)
  case 2
    %  test for example 2
    [A,B,C,N] = getSystem2();
    x0 = [0.25;-0.25];

    %  Set the options for ode integration
    T   = 400;  % final time for integration of the future energy
    options = odeset('AbsTol',1e-15);%,'NonNegative',n+1);

    d   = 7;    % maximum degree of energy function approximations

    feedbackDegree = 7;  % degree of energy function approximation for 
                         % defining the control input.
    
    % requires eta=0
    analyticSolutionFEF = (5*x0(2)^2)/8 - (11*x0(2)^3)/36 + x0(2)^4/24 + ...
                          (3*x0(1)*x0(2))/4 - (x0(1)*x0(2)^2)/6 + x0(1)^2/4;
    fprintf('For example 2, the future energy function is %12.8e\n',...
            analyticSolutionFEF)
  case 3
    %  test for example 3: nominal T=4, n=16, m=4, p=2, d=feedback=5, tol=10^-15
    %                      epsilon=0.1, alpha=0
    epsilon = 0.1;
    alpha = 0;
    n = 16;
    m = 4;
    p = 2;
    [A,B,C,N,zInit,M] = getSystem3(n,m,p,epsilon,alpha);
    %x0 = 0.02*zInit;    % gives a very interesting convergence pattern
    %x0 = 0.001*zInit;   % matches the initial condition studied in the paper
    x0 = 0.0001*zInit;   % gives reasonable agreement with past energy
                         

    %  Set the options for ode integration
    T   = 8.00;  % final time for integration of the future energy
    options = odeset('AbsTol',1e-15);%,'NonNegative',n+1);

    d   = 5;    % maximum degree of energy function approximations

    feedbackDegree = 5;  % degree of energy function approximation for 
                         % defining the control input.

  case 4
    %  test for example 4: nominal n=8, m=4, p=4, L=6, T=4, tol=1e-15,
    %  z_factor=0.01, d=5, feedback=0, z0=sin(4\pi x)
    n = 8;
    m = 4;
    p = 4;
    L = 6; %13.0291;

    z_factor = 0.01;  % scale factor on the initial conditions.

    [A,B,C,N,zInit,M] = getSystem4(n,m,p,1/L^2);
    x0 = z_factor*zInit;    % gives a very interesting convergence pattern

    %  Set the options for ode integration
    T   = 40.0;  % final time for integration of the future energy
    options = odeset('AbsTol',1e-9);%,'NonNegative',n+1);

    d   = 5;             % maximum degree of energy function approximations

    feedbackDegree = 2;  % degree of energy function approximation for 
                         % defining the control input.
    testPastEnergyFunction   = false; % the past energy function doesn't 
                                      % generate a stabilizing control for
                                      % this n,m,p,L

  otherwise
    error('testCase is not implemented')

end

n = size(A,1);


%% perform future energy function test
if (testFutureEnergyFunction)
  [w] = approxFutureEnergy(A,N,B,C,eta,d);  %#ok
  
  rhs_o = @(t,x) [ A*x(1:n)+N*kron(x(1:n),x(1:n)); 0.5*x(1:n).'*(C.'*C)*x(1:n) ];
%  [time,X] = ode45(rhs_o,[0 T],[x0;0],options);
  [time,X] = ode23s(rhs_o,[0 T],[x0;0],options);
  
  integratedOutputEnergy = X(end,end);
  for i=1:d
    w{i} = w{i}.';
  end
  fprintf('The integral approximation to (9) is %12.8e\n',integratedOutputEnergy)
  for degree=2:d
    approximatedOutputEnergy = 0.5*kronPolyEval(w,x0,degree);
    fprintf('The degree %d approximation to (9) is %12.8e\n',degree,approximatedOutputEnergy)
  end
  fprintf('\n\n')
end

%% perform past energy function test
%
%  This test requires the exact past energy function to determine the control
%  input.  Since that is not available, here we are replacing it with a 
%  low-degree polynomial approximation.  Hence, we test initial points that
%  are even closer to the origin.
if (testPastEnergyFunction)
  [v] = approxPastEnergy(A,N,B,C,eta,d,true);  %#ok
  
  % these controls are set up approximately from polynomial approx to Ec
  switch feedbackDegree

    case 2
      u2 = @(x,v,B) 0.5*( 2*v{2}.'*kron(B,x) ).';
      rhs2_c = @(t,x) [-A*x(1:n)-N*kron(x(1:n),x(1:n))-B*u2(x(1:n),v,B); 0.5*norm(u2(x(1:n),v,B))^2];
      [time,X] = ode45(rhs2_c,[0 T],[x0;0],options);

    case 3
      u3 = @(x,v,B) 0.5*( 2*v{2}.'*kron(B,x) + ...
                          3*v{3}.'*kron(B,kron(x,x)) ).';
      rhs3_c = @(t,x) [-A*x(1:n)-N*kron(x(1:n),x(1:n))-B*u3(x(1:n),v,B); 0.5*norm(u3(x(1:n),v,B))^2];
      [time,X] = ode45(rhs3_c,[0 T],[x0;0],options);

    case 4
      u4 = @(x,v,B) 0.5*( 2*v{2}.'*kron(B,x) + ...
                          3*v{3}.'*kron(B,kron(x,x)) + ...
                          4*v{4}.'*kron(B,kron(x,kron(x,x))) ).';
      rhs4_c = @(t,x) [-A*x(1:n)-N*kron(x(1:n),x(1:n))-B*u4(x(1:n),v,B); 0.5*norm(u4(x(1:n),v,B))^2];
      [time,X] = ode45(rhs4_c,[0 T],[x0;0],options);

    case 5
      u5 = @(x,v,B) 0.5*( 2*v{2}.'*kron(B,x) + ...
                          3*v{3}.'*kron(B,kron(x,x)) + ...
                          4*v{4}.'*kron(B,kron(x,kron(x,x))) + ...
                          5*v{5}.'*kron(B,kron(x,kron(x,kron(x,x)))) ).';
      rhs5_c = @(t,x) [-A*x(1:n)-N*kron(x(1:n),x(1:n))-B*u5(x(1:n),v,B); 0.5*norm(u5(x(1:n),v,B))^2];
      [time,X] = ode45(rhs5_c,[0 T],[x0;0],options);

    case 6
      u6 = @(x,v,B) 0.5*( 2*v{2}.'*kron(B,x) + ...
                          3*v{3}.'*kron(B,kron(x,x)) + ...
                          4*v{4}.'*kron(B,kron(x,kron(x,x))) + ...
                          5*v{5}.'*kron(B,kron(x,kron(x,kron(x,x)))) + ...
                          6*v{6}.'*kron(B,kron(x,kron(x,kron(x,kron(x,x)))))).';
      rhs6_c = @(t,x) [-A*x(1:n)-N*kron(x(1:n),x(1:n))-B*u6(x(1:n),v,B); 0.5*norm(u6(x(1:n),v,B))^2];
      [time,X] = ode45(rhs6_c,[0 T],[x0;0],options);
  
    case 7
      u7 = @(x,v,B) 0.5*( 2*v{2}.'*kron(B,x) + ...
                          3*v{3}.'*kron(B,kron(x,x)) + ...
                          4*v{4}.'*kron(B,kron(x,kron(x,x))) + ...
                          5*v{5}.'*kron(B,kron(x,kron(x,kron(x,x)))) + ...
                          6*v{6}.'*kron(B,kron(x,kron(x,kron(x,kron(x,x))))) + ...
                          7*v{7}.'*kron(B,kron(x,kron(x,kron(x,kron(x,kron(x,x))))))).';
      rhs7_c = @(t,x) [-A*x(1:n)-N*kron(x(1:n),x(1:n))-B*u7(x(1:n),v,B); 0.5*norm(u7(x(1:n),v,B))^2];
      [time,X] = ode45(rhs7_c,[0 T],[x0;0],options);
      
    otherwise
      warning('not implemented, without v{8} symmetry, this is a long expr.')

  end


  integratedControlEnergy = X(end,end);
  for i=1:d
    v{i} = v{i}.';
  end
  fprintf('The integral approximation to (8) is %12.8e\n',integratedControlEnergy)
  for degree=2:d
    approximatedControlEnergy = 0.5*kronPolyEval(v,x0,degree);
    fprintf('The degree %d approximation to (8) is %12.8e\n',degree,approximatedControlEnergy)
  end
  fprintf('\n\n')
end

%%  Visualize the controlled solution if the past energy function was tested
if (testPastEnergyFunction && visualize)
  XX = sqrtm(full(M))\X(:,1:n).';
  mesh(XX(1:2:n,:))
end
