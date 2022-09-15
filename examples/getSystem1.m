function [A,B,C,N] = getSystem1()
%getSystem1  Generates a simple quadratic system for testing energy functions.
%
%   Usage:  [A,B,C,N] = runExaple1()
%
%   The "matrices" correspond to the quadratic input-output system
%
%       \dot(x) = -2x + x^2 + 2u
%             y = 2x
%
%   for which there is an analytic solution to the past and future energy
%   functions.
%
%   Reference: Nonlinear Balanced Truncation Model Reduction: 
%        Part 1-Computing Energy Functions, by Kromer, Gugercin, and Borggaard.
%        arXiv.
%
%   Part of the NLbalancing repository.
%%

    A = -2;
    B =  2;
    C =  2;
    N =  1;

end
