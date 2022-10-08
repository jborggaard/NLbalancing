%  A script to run the examples in arXiv:2209.07645
setKroneckerToolsPath
addpath('examples')
addpath('utils')

%%  runExample1 produces data for Fig. 1.
%runExample1

%%  runExample2 produces the plots for Fig. 2.  
%  For Table 1, please use energyFunctionValidation in the tests directory.
%[v,w] = runExample2(7,true,true,6);

%%  runExample3 produces Tables II and III.
%runExample3

%%  runExample4 produces Tables IV and V.
runExample4
