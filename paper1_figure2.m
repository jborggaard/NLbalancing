%  A script to run an example from ``Scalable Computation of Energy Functions
%  for Nonlinear Balanced Truncation,'' CMAME 427(117011), 2024.
%  by Kramer, Gugercin, Borggaard, and Balicki.

setKroneckerToolsPath
addpath('examples')
addpath('utils')

%%  runExample2 produces the plots for Fig. 2.  
%  For Table 1, please use energyFunctionValidation in the tests directory.
[v,w] = runExample2(4,true,false,3);  

%  The plots in the paper are actually for degree 4 approximations to the energy
%  functions.  To plot the degree 6 approximations referenced in the paper, use
%  the following line instead:
%
%  [v,w] = runExample2(6,true,false,5);
