function [v,w] = runExample1()
%runExample1 Runs the first example from the paper
%   Usage:  [] = runExaple1()
%
%   
%
%   Reference: Nonlinear Balanced Truncation Model Reduction: 
%        Part 1-Computing Energy Functions, by Kromer, Gugercin, and Borggaard.
%        arXiv.
%
%   Part of the NLbalancing repository.
%%

  [A,B,C,N] = getSystem1();

  eta = 0.5;    % values should be between -\infty and 1.
                % eta=0.5 corresponds to gamma= sqrt(2)
                % since eta = 1 - 1/gamma^2;

  %  Compute the polynomial approximations to the future energy function
  [w] = approxFutureEnergy(A,N,B,C,eta,8);
  w2 = w{2}; w3 = w{3}; w4 = w{4}; w5 = w{5}; w6 = w{6}; w7 = w{7}; w8 = w{8};

  x = linspace(-6,6,121);
    
  Ef2 =       0.5*w2*x.^2;
  Ef3 = Ef2 + 0.5*w3*x.^3;
  Ef4 = Ef3 + 0.5*w4*x.^4;
  Ef5 = Ef4 + 0.5*w5*x.^5;
  Ef6 = Ef5 + 0.5*w6*x.^6;
  Ef7 = Ef6 + 0.5*w7*x.^7;
  Ef8 = Ef7 + 0.5*w8*x.^8;
    
  %  Compute the analytical solution for comparison
  EPlusAnalytic = EgammaPlus(x,A,B,C,N,eta);

  plot(x,EPlusAnalytic,'+',...
       x,Ef2,...
       x,Ef4,...
       x,Ef6,...
       x,Ef8,...
       'LineWidth',2)
  legend('analytic',...
         'degree 2',...
         'degree 4',...
         'degree 6',...
         'degree 8',...
         'Location','northwest')
  xlabel('$x$',...
         'interpreter','latex',...
         'FontSize', 20,...
         'fontweight','bold')
  ylabel('$\mathcal{E}_\gamma^+$',...
         'interpreter','latex',...
         'FontSize', 20,...
         'fontweight','bold')

  %  Save data to generate tikz plots for the paper
  fid = fopen('plots/ex1_future_a.txt','w');
  fprintf(fid,'%g %g\n',[x;EPlusAnalytic]);
  fclose(fid);

  fid = fopen('plots/ex1_future_2.txt','w');
  fprintf(fid,'%g %g\n',[x;Ef2]);
  fclose(fid);

  fid = fopen('plots/ex1_future_4.txt','w');
  fprintf(fid,'%g %g\n',[x;Ef4]);
  fclose(fid);

  fid = fopen('plots/ex1_future_6.txt','w');
  fprintf(fid,'%g %g\n',[x;Ef6]);
  fclose(fid);

  fid = fopen('plots/ex1_future_8.txt','w');
  fprintf(fid,'%g %g\n',[x;Ef8]);
  fclose(fid);


  %  Compute the polynomial approximations to the past energy function
  [v] = approxPastEnergy(A,N,B,C,eta,8);
  v2 = v{2}; v3 = v{3}; v4 = v{4}; v5 = v{5}; v6 = v{6}; v7 = v{7}; v8 = v{8};

  x = linspace(-6,6,121);
    
  Ep2 =       0.5*v2*x.^2;
  Ep3 = Ep2 + 0.5*v3*x.^3;
  Ep4 = Ep3 + 0.5*v4*x.^4;
  Ep5 = Ep4 + 0.5*v5*x.^5;
  Ep6 = Ep5 + 0.5*v6*x.^6;
  Ep7 = Ep6 + 0.5*v7*x.^7;
  Ep8 = Ep7 + 0.5*v8*x.^8;

  %  Compute the analytical solution for comparison
  EMinusAnalytic = EgammaMinus(x,A,B,C,N,eta);

  figure
  plot(x,EMinusAnalytic,'+',...
       x,Ep2,...
       x,Ep4,...
       x,Ep6,...
       x,Ep8,...
       'LineWidth',2)
  legend('analytic',...
         'degree 2',...
         'degree 4',...
         'degree 6',...
         'degree 8',...
         'Location','northeast')
  xlabel('$x$',...
         'interpreter','latex',...
         'FontSize', 20,...
         'fontweight','bold')
  ylabel('$\mathcal{E}_\gamma^-$',...
         'interpreter','latex',...
         'FontSize', 20,...
         'fontweight','bold')

  %
  fid = fopen('plots/ex1_errorpe_2.txt','w');
  fprintf(fid,'%g %g\n',[x;abs(Ep2-EMinusAnalytic)]);
  fclose(fid);

  fid = fopen('plots/ex1_errorpe_4.txt','w');
  fprintf(fid,'%g %g\n',[x;abs(Ep4-EMinusAnalytic)]);
  fclose(fid);

  fid = fopen('plots/ex1_errorpe_6.txt','w');
  fprintf(fid,'%g %g\n',[x;abs(Ep6-EMinusAnalytic)]);
  fclose(fid);

  fid = fopen('plots/ex1_errorpe_8.txt','w');
  fprintf(fid,'%g %g\n',[x;abs(Ep8-EMinusAnalytic)]);
  fclose(fid);

  %
  fid = fopen('plots/ex1_past_a.txt','w');
  fprintf(fid,'%g %g\n',[x;EMinusAnalytic]);
  fclose(fid);

  fid = fopen('plots/ex1_past_2.txt','w');
  fprintf(fid,'%g %g\n',[x;Ep2]);
  fclose(fid);

  fid = fopen('plots/ex1_past_4.txt','w');
  fprintf(fid,'%g %g\n',[x;Ep4]);
  fclose(fid);

  fid = fopen('plots/ex1_past_6.txt','w');
  fprintf(fid,'%g %g\n',[x;Ep6]);
  fclose(fid);

  fid = fopen('plots/ex1_past_8.txt','w');
  fprintf(fid,'%g %g\n',[x;Ep8]);
  fclose(fid);


  save('Ex1_RawData.mat','x',...
       'EPlusAnalytic', 'Ef2','Ef3','Ef4','Ef5','Ef6','Ef7','Ef8',...
       'EMinusAnalytic','Ep2','Ep3','Ep4','Ep5','Ep6','Ep7','Ep8')
end


function [E] = EgammaPlus(x,a,b,c,n,eta)
%  Taking the upper bound of E_gamma^+
%  E =  1/(b^2 * eta)*( (a*b^2*c^2*eta*sqrt(x.^2.*((a + n*x).^2 + b^2*c^2*eta)).* ...
%     log(sqrt((a + n*x).^2 + b^2*c^2*eta) + a + n*x))./ ...
%     (2*n^2*x.*sqrt((a + n*x).^2 + b^2*c^2*eta)) ...
%    - ( sqrt(x.^2.*((a + n*x).^2 + b^2*c^2*eta)) .* ...
%    ((a + n*x).^2/(3*n) - a*(a + n*x)/(2*n) + b^2*c^2*eta/(3*n)))./(n*x)...
%    + (a*x.^2)/2 + (n*x.^3)/3 ); 

 E= 1/(b^2 * eta)*(a*b^2*c^2*eta*log(a + sqrt(a^2 + b^2*c^2*eta)))/(2*n^2) + ...
    1/(b^2 * eta)*(sqrt(a^2 + b^2*c^2*eta)*(a^2 - 2*b^2*c^2*eta))/(6*n^2) + ...
    1/(b^2 * eta)*(sqrt(a^2 + b^2*c^2*eta + 2*a*n*x + n^2*x.^2).*(-a^2 + 2*b^2*c^2*eta + a*n*x + 2*n^2*x.^2))/(6*n^2) - ...
    1/(b^2 * eta)*(a*b^2*c^2*eta*log(a + n*x + sqrt(a^2 + b^2*c^2*eta + 2*a*n*x + n^2*x.^2)))/(2*n^2) + ...
    1/(b^2 * eta)*((a*x.^2)/2 + (n*x.^3)/3 ); 
end


function [E] = EgammaMinus(x,a,b,c,n,eta)

 E = (log(sqrt(a^2 + b^2*c^2*eta + 2*a*n*x + n^2*x.^2) + (a + n*x)).*(8*a^3*n^3 - 8*a*n^3*(a^2 + b^2*c^2*eta)))/(16*b^2*n^5) - ...
     (n*x.^3)/(3*b^2) - (a*x.^2)/(2*b^2) - ...
     (log(a + sqrt(a^2 + b^2*c^2*eta))*(8*a^3*n^3 - 8*a*n^3*(a^2 + b^2*c^2*eta)))/(16*b^2*n^5) + ...
     (((n^2*(a^2 + b^2*c^2*eta + n^2*x.^2))/3 - (a^2*n^2)/2 + (a*n^3*x)/6).*sqrt(a^2 + b^2*c^2*eta + 2*a*n*x + n^2*x.^2))/(b^2*n^4) - ...
     (sqrt(a^2 + b^2*c^2*eta)*((n^2*(a^2 + b^2*c^2*eta))/3 - (a^2*n^2)/2))/(b^2*n^4);
 
end
