%  A script to run the examples in Nonlinear Balanced Truncation: 
%  Part 2-Model Reduction on Manifolds, Kramer, Gugercin, and Borggaard.
%
%  testcase=1 produces Table 1.

setKroneckerToolsPath
addpath('examples')
addpath('utils')

testcase = 1;

%% Example 1
switch testcase
  case 0 % an interesting random 2-dof example
    % set up the test parameters
    validateInputNormalPastEnergy = true;
    makeEnergyPlots = false;

    maxDegree = 7;

    eta = 0.2;    % values should be between -\infty and 1.

    zmin =-0.2;%-0.1; 
    zmax = 0.2;%0.1; 
    Npts = 51;

    rng(20,'twister')
    n = 2;
    m = 2;
    p = 2;
    A = rand(n,n);
    N = rand(n,n^2);
    B = rand(n,m);
    C = rand(p,n);

    %
    % approximate the energy functions 
    [w] = approxFutureEnergy(A,N,B,C,eta,maxDegree+1);
    [v] = approxPastEnergy(A,N,B,C,eta,maxDegree+1);

    nSingValFcns = 2;

  case 1 % Corresponds to example 1 from the second paper
    % set up the test
    validateInputNormalPastEnergy = true;
    makeEnergyPlots = true;

    maxDegree = 8;
    
    eta = 0.1;    % eta = 0.1 matches Example 2 in the first paper
                  % values should be between -\infty and 1.

    zmin =-0.2; 
    zmax = 0.2; 
    Npts = 51;

    % get the system
    [A,B,C,N] = getSystem2();
    n = size(A,1);

    %
    % approximate the energy functions 
    [w] = approxFutureEnergy(A,N,B,C,eta,maxDegree+1);
    [v] = approxPastEnergy(A,N,B,C,eta,maxDegree+1);

    nSingValFcns = 2;

  case 2
    % set up the test
    validateInputNormalPastEnergy = false;
    makeEnergyPlots = false;

    maxDegree = 5;

    eta =  0.99;     % these parameter values match part 1

    zmin =-0.0002; 
    zmax = 0.0002; 
    Npts = 51;

    m = 4;
    p = 4;

    n = 16;

    epsilon = 0.001;
    alpha   = 0.0;

    % get the system
    [A,B,C,N,zInit] = getSystem3(n,m,p,epsilon,alpha);

    %
    % approximate the energy functions 
    [w] = approxFutureEnergy(A,N,B,C,eta,maxDegree+1);
    [v] = approxPastEnergy(A,N,B,C,eta,maxDegree+1);
    
    nSingValFcns = 10;

  otherwise

end

  %
  % compute the input-normal transformation approximation
  verbose = true;
  [sigma,T] = inputNormalTransformation(v,w,maxDegree,verbose);


  %
  %% Some convenient definitions
  dz   = (zmax-zmin)/(Npts-1);
  zRange = linspace(zmin,zmax,Npts);

  vT = cell(size(v));  % calculate transposes for convenience
  for i=2:length(v)
    vT{i} = v{i}.';
  end

  wT = cell(size(w));  % calculate transposes for convenience
  for i=2:length(w)
    wT{i} = w{i}.';
  end
  


  if (validateInputNormalPastEnergy || makeEnergyPlots)
    %
    %  Plot the past energy function in a neighborhood of the origin
    %
  
    %  first in the original coordinates
    Eplot = zeros(Npts,Npts);
    for i=1:Npts
      for j=1:Npts
        z = [ zRange(i); zRange(j) ];
        Eplot(i,j) = 0.5*kronPolyEval(vT,z,maxDegree+1);
      end
    end
    if (makeEnergyPlots)
      figure(10)
      surf(zRange,zRange,Eplot)
      title('Past energy function in original coordinates')
    end
  
    %  the ideal past energy function for input-normal form should be
    %  0.5\|z\|.  Thus,
    Epast = zeros(Npts,Npts);
    for i=1:Npts
      for j=1:Npts
        z = [ zRange(i); zRange(j) ];
        Epast(i,j) = 0.5*(z.'*z);
      end
    end
  
    %  then in the balanced coordinates
    for degree=1:maxDegree
      Eplot = zeros(Npts,Npts);
      for i=1:Npts
        for j=1:Npts
          z = [ zRange(i); zRange(j) ];
          x = kronPolyEval(T,z,degree);
          Eplot(i,j) = 0.5*kronPolyEval(vT,x,maxDegree+1);
        end
      end
      if ( makeEnergyPlots )
        figure(degree)
        surf(zRange,zRange,Eplot)
        string = sprintf('Past energy function using degree %d transformation',degree);
        title(string)
      end
 
      if ( validateInputNormalPastEnergy )
        EminusError = max(max(abs(Eplot-Epast)));
        fprintf('Transformed past energy function error at degree %d is %g\n',...
                degree,EminusError);
        %[g,i] = max(max(abs(Eplot-Epast)))
      end
    end
  
    %
    %  Plot the future energy function in a neighborhood of the origin
    %
  
    %  first in the original coordinates
    Eplot = zeros(Npts,Npts);
    for i=1:Npts
      for j=1:Npts
        z = [ zRange(i); zRange(j) ];
        Eplot(i,j) = 0.5*kronPolyEval(wT,z,maxDegree+1);
      end
    end
    if (makeEnergyPlots)
      figure(20)
      surf(zRange,zRange,Eplot)
      title('Future energy function in original coordinates')
    end
  
    %  then in the balanced coordinates
    for degree=1:maxDegree
      Eplot = zeros(Npts,Npts);
      for i=1:Npts
        for j=1:Npts
          z = [ zRange(i); zRange(j) ];
          x = kronPolyEval(T,z,degree);
          Eplot(i,j) = 0.5*kronPolyEval(wT,x,maxDegree+1);
        end
      end
      if (makeEnergyPlots)
        figure(10+degree)
        surf(zRange,zRange,Eplot)
        string = sprintf('Future energy function using degree %d transformation',degree);
        title(string)
      end
    end
  end

  
  %
  %% Approximate the singular value functions using Algorithm 2.
  [c] = approximateSingularValueFunctions(T,w,sigma,maxDegree-1);

  %
  %% Generate data for plots of singular value functions
  zRange = linspace(-1e-5,1e-5,51);
  plotSingularValueFunctions(sigma,c,zRange,nSingValFcns)


