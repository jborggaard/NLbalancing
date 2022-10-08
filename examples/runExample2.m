function [v,w] = runExample2(degree,plotEnergy,plotBalancing,balancingDegree)
%EXAMPLE2 Runs the second example from the paper
%
%   Usage:  [v,w] = runExample2(degree,plotEnergy,plotBalancing,balancingDegree)
%
%   where
%         degree          is the degree of energy function approximations
%         plotEnergy      is a logical variable to determine if a plot is made.
%
%         plotBalancing   is a logical variable to determine if a plot is made.
%         balancingDegree is a small integer setting the degree of approximation
%                         in the balancing transformation. Must be < degree.
%                         (used if plotBalancing = true)
%
%         v,w             are coefficients of the past and future energy
%                         function approximations, respectively.
%
%   The value of eta is set below.
%
%   Reference: Nonlinear Balanced Truncation Model Reduction: 
%        Part 1-Computing Energy Functions, by Kramer, Gugercin, and Borggaard.
%        arXiv:2209.07645.
%
%   This example is motivated by Kawano and Scherpen, IEEE Transactions
%   on Automatic Control, 2016.  Here we ignore the bilinear term 2*x_2*u.
%   
%   Part of the NLbalancing repository.
%%

  fprintf('Running Example2\n')
  
  eta = 0.1;    % values should be between -\infty and 1.
                % eta=0.1 corresponds to gamma= 1.0541...
                % since eta = 1 - 1/gamma^2;

  fprintf('Simulating for eta=%g (gamma=%g)\n',eta,1/sqrt(1-eta))


  if (nargin<1)
    degree = 4;
    plotEnergy = true;

    plotBalancing = false;
    balancingDegree = 3;
  end


  if ( plotBalancing )
    dataRange = 0.35; %6.0
  else
    dataRange = 0.35; %0.75;
  end
  
  [A,B,C,N] = getSystem2();
 
  %  Compute the polynomial approximations to the future energy function
  [w] = approxFutureEnergy(A,N,B,C,eta,degree,true);
  futureEnergy{degree} = [];
  for k=2:degree
    futureEnergy{k} = w{k}.'/2;
  end

  %  Plot the future energy function for this example
  %disp('Future energy function for example 2')    

  if (plotEnergy || plotBalancing)
    nX = 101; nY = 101;
    xPlot = linspace(-dataRange,dataRange,nX);
    yPlot = xPlot;
    eFuture = zeros(nX,nY);
    [X,Y] = meshgrid(xPlot,yPlot);
    for i=1:nX
      for j=1:nY
        x = [X(i,j);Y(i,j)];
%         WBar = vbar(futureEnergy,x);
%         eFuture(i,j) = x.'*WBar*x;
        eFuture(i,j) = kronPolyEval(futureEnergy,x,degree);
      end
    end
    figure(1)
    contourf(X,Y,eFuture)
    xlabel('$x_1$','interpreter','latex'); 
    ylabel('$x_2$','interpreter','latex'); 
    colorbar('FontSize',16)
    ax = gca;
    set(gca,'FontSize',20)
    exportgraphics(ax,'FEF_p0_1.pdf','ContentType','vector');
    title('Future Energy Function')
  end

  [v] = approxPastEnergy(A,N,B,C,eta,degree,true);
  pastEnergy{degree} = [];
  for k=2:degree
    pastEnergy{k} = v{k}.'/2;
  end
    
  %  Plot the past energy function for this example
  %disp('Past energy function for example 2')    

  if (plotEnergy || plotBalancing)
    nX = 101; nY = 101;
    xPlot = linspace(-dataRange,dataRange,nX);
    yPlot = xPlot;
    ePast = zeros(nX,nY);
    [X,Y] = meshgrid(xPlot,yPlot);
  
  
    for i=1:nX
      for j=1:nY
        x = [X(i,j);Y(i,j)];
%         VBar = vbar(pastEnergy,x);
%         ePast(i,j) = x.'*VBar*x;
        ePast(i,j) = kronPolyEval(pastEnergy,x,degree);
      end
    end
    figure(2)
    contourf(X,Y,ePast)
%    mesh(X,Y,ePast)
    xlabel('$x_1$','interpreter','latex'); 
    ylabel('$x_2$','interpreter','latex');
    colorbar('FontSize',16)
    ax = gca;
    set(gca,'FontSize',20)
    exportgraphics(ax,'PEF_p0_1.pdf','ContentType','vector');
    title('Past Energy Function')
  end

  save('Ex2_RawData.mat','v','w')


  if ( plotBalancing )
    [sigma,T] = inputNormalTransformation(v,w,balancingDegree);
    nPts = 201;
    s = linspace(-2,2,nPts);
    lin = T{1}(:,1)*s;
    
    coord = lin;
    for k=2:balancingDegree
      coord = coord + T{k}(:,1)*s.^k;
    end
    
    idxLin = zeros(1,nPts);
    linCount = 0;
    for i=1:nPts
      if (norm(lin(:,i),inf)<dataRange)
        linCount = linCount + 1;
        idxLin(linCount) = i;
      end
    end
    idxLin = idxLin(1:linCount);
    
    idxCoord = zeros(1,nPts);
    coordCount = 0;
    for i=1:nPts
      if (norm(coord(:,i),inf)<dataRange)
        coordCount = coordCount + 1;
        idxCoord(coordCount) = i;
      end
    end
    idxCoord = idxCoord(1:coordCount);
    
    figure(1); hold on
    plot(lin(1,idxLin),lin(2,idxLin),'w+')
    plot(coord(1,idxCoord),coord(2,idxCoord),'r+')

    figure(2); hold on
    plot(lin(1,idxLin),lin(2,idxLin),'w+')
    plot(coord(1,idxCoord),coord(2,idxCoord),'r+')
  end
end


