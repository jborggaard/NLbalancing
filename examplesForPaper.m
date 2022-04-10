function [] = examplesForPaper(ex_flag,plot_flag)
%  examples for paper
  setKroneckerToolsPath
  addpath('examples')
  addpath('utils')
  
  if ( nargin<1 )
    ex_flag = 1;
  end
  if ( nargin<2 )
    plot_flag = false;
  end


  if ( ex_flag==1 )
    %---------------------------------------------------------------------------
    %  Example 1
    %---------------------------------------------------------------------------
    A = -2;
    B = 2;
    C = 2;
    N = 1;
    eta = 1;    % values should be between -\infty and 1.
    %eta = 1 - 1/gamma^2;

    [v] = solveFutureEnergy(A,N,B,C,eta,8);
    v2 = v{2}; v3 = v{3}; v4 = v{4}; v5 = v{5}; v6 = v{6}; v7 = v{7}; v8 = v{8};
    
    % division by 2 to match equation (40) in our definition of E^+
    Eplus = Polynom([v8 v7 v6 v5 v4 v3 v2 0 0]/2);
    disp('Future energy function for example 1')
    disp(Eplus)

    if ( plot_flag )
      figure(11)
      hold on
      plot(Eplus)
     % set(gca,'yscale','log')
      legend('Analytical','Polynomial')
    end
    [w] = solvePastEnergy(A,N,B,C,eta,8);
    w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};
    
    Eminus = Polynom([w8 w7 w6 w5 w4 w3 w2 0 0]/2);
    disp('Past energy function for example 1')
    disp(Eminus)

    if ( plot_flag )
      figure(12)
      hold on
      plot(Eminus)
      legend('Analytical','Polynomial')
    end
    
  elseif ( ex_flag==2 )
    %---------------------------------------------------------------------------
    %  Example 2 (made up example in examples_energyFcts_1d.m)
    %---------------------------------------------------------------------------
    A = -1;
    B = 1;
    C = 1;
    N = 2;
    eta = 1;    % values should be between -\infty and 1.
    %eta = 1 - 1/gamma^2;

    [v] = solveFutureEnergy(A,N,B,C,eta,8);
    v2 = v{2}; v3 = v{3}; v4 = v{4}; v5 = v{5}; v6 = v{6}; v7 = v{7}; v8 = v{8};
    
    % division by 2 to match equation (40) in our definition of E^+
    Eplus = Polynom([v8 v7 v6 v5 v4 v3 v2 0 0]/2);
    disp('Future energy function for example 2')
    disp(Eplus)

    [w] = solvePastEnergy(A,N,B,C,eta,8);
    w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};
    
    Eminus = Polynom([w8 w7 w6 w5 w4 w3 w2 0 0]/2);
    disp('Past energy function for example 2')
    disp(Eminus)

  elseif ( ex_flag==3 )
    %---------------------------------------------------------------------------
    %  Example 3
    %---------------------------------------------------------------------------
    A = [-1 1;0 -1];
    N = [0 0 0 -1;0 0 0 0];
    B = [1;1];
    C = [1 1];
    eta = 0.1;
    fprintf('Simulating for eta=%g (gamma=%g)\n',eta,1/sqrt(1-eta))
    n=2;
    
    [v] = solveFutureEnergy(A,N,B,C,eta,4);
    v2=v{2}/2;
    v3=v{3}/2;
    v4=v{4}/2;

    %  Plot the future energy function for this example
    disp('Future energy function for example 3')    

    nX = 101; nY = 101;
    xPlot = linspace(-1,1,nX);
    yPlot = xPlot;
    eFuture = zeros(nX,nY);
    [X,Y] = meshgrid(xPlot,yPlot);
    for i=1:nX
      for j=1:nY
        x = [X(i,j);Y(i,j)];
        VBar = vbar({[],v2,v3,v4},x);
        eFuture(i,j) = x.'*VBar*x;
      end
    end
    figure(1)
    contourf(X,Y,eFuture)
    xlabel('x_1'); ylabel('x_2'); colorbar
    
    [w] = solvePastEnergy(A,N,B,C,eta,4);
    w2=w{2}/2;
    w3=w{3}/2;
    w4=w{4}/2;
    
    %  Plot the past energy function for this example
    disp('Past energy function for example 3')    

    nX = 101; nY = 101;
    xPlot = linspace(-1,1,nX);
    yPlot = xPlot;
    ePast = zeros(nX,nY);
    [X,Y] = meshgrid(xPlot,yPlot);
    
    
    for i=1:nX
      for j=1:nY
        x = [X(i,j);Y(i,j)];
        WBar = vbar({[],w2,w3,w4},x);
        ePast(i,j) = x.'*WBar*x;
      end
    end
    figure(2)
    contourf(X,Y,ePast)
    xlabel('x_1'); ylabel('x_2'); colorbar
    
  elseif ( ex_flag==4 )
    QQRfolder = '/Volumes/borggaard/Software/MyPublicSoftware/QQR/';
    addpath([QQRfolder,'examples/Burgers1DControl'])
    
    %  Set problem dimensions
    n = 60;
    m = 2;
    
    [M,A,B,N,zInit] = BurgersFEMControl(n,m);
    
    % add linear reaction term
    alpha = 0.0;
    epsilon = 0.001;
    
    A = epsilon*A + alpha*M;
    Q = M;
    
    C = (1/n)*ones(1,n);

    %  Perform a change of variables to eliminate the positive definite
    %  mass matrix.  M^(1/2)z -> z
    %
    %    \dot{z} = Az+Bu+N*kron(z,z),
    %
    %  This is required until we extend qqr to handle "mass matrices"
    scaling = true;
  
    if ( scaling )
      sqM = sqrtm(full(M));
      sqMinv = inv(sqM);

      Ac = sqMinv*A*sqMinv;                  %#ok
      Bc = sqMinv*B;                         %#ok
      Nc = kroneckerRight(sqMinv*N,sqMinv);  %#ok
      Cc = C*sqMinv;                         %#ok
      Qc = eye(n);  R = eye(m);
    else
      %  use a clunky way to put it into standard form
      Ac = M\A;
      Bc = M\B;
      Nc = M\N;
      Cc = C;
      Qc = M;       R = eye(m);
    end
    
    eta = 0.1;
    fprintf('Simulating for eta=%g (gamma=%g)\n',eta,1/sqrt(1-eta))
    n=2;
    
    tic;
    [v] = solveFutureEnergy(Ac,Nc,Bc,Cc,eta,4);
    t = toc;
    fprintf('CPU time is %g\n',t)
    v2=v{2}/2;
    v3=v{3}/2;
    v4=v{4}/2;


    
  else
    error('example is not set up')
    
  end

end
