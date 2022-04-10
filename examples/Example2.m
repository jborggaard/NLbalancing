function [v,w] = Example2(degree,plotEnergy,plotBalancing,balancingDegree)
%EXAMPLE2 Runs the second example from the paper
%
%   Usage:   [v,w] = Example2(degree,plotEnergy,plotBalancing,balancingDegree)
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
%   Nonlinear Balanced Truncation Model Reduction for Large-Scale 
%   Polynomial Systems, by Kramer, Borggaard, and Gugercin, arXiv.
%
%   This example is motivated by Kawano and Scherpen, IEEE Transactions
%   on Automatic Control, 2016.  Here we ignore the bilinear term 2*x_2*u.
%   
%   Part of the NLbalancing repository.
%%

  if (nargin<1)
    degree = 6;
    plotEnergy = true;

    plotBalancing = false;
    balancingDegree = 5;
  end


  if ( plotBalancing )
    dataRange = 6;
  else
    dataRange = 0.75;
  end
  
  A = [-1 1;0 -1];
  N = [0 0 0 -1;0 0 0 0];
  B = [1;1];
  C = [1 1];
 
  eta = 0.1;    % values should be between -\infty and 1.
                % eta=0.1 corresponds to gamma= 1.0541...
                % since eta = 1 - 1/gamma^2;

  fprintf('Simulating for eta=%g (gamma=%g)\n',eta,1/sqrt(1-eta))

  %  Compute the polynomial approximations to the future energy function
  [w] = approxFutureEnergy(A,N,B,C,eta,degree);
  futureEnergy{degree} = [];
  for k=2:degree
    futureEnergy{k} = w{k}/2;
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
        WBar = vbar(futureEnergy,x);
        eFuture(i,j) = x.'*WBar*x;
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

  [v] = approxPastEnergy(A,N,B,C,eta,degree);
  pastEnergy{degree} = [];
  for k=2:degree
    pastEnergy{k} = v{k}/2;
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
        VBar = vbar(pastEnergy,x);
        ePast(i,j) = x.'*VBar*x;
      end
    end
    figure(2)
    contourf(X,Y,ePast)
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


