# ILEQS - Iterative optimization solution of linear equations
%
% Input:
%   Ab       - Augmented matrix of linear equations, N * (N + 1) double matrix
%   method   - Method of iteration: 'Jacob', 'G-S', 'SOR' or 'Adaptive'
%   options  - Options of algorithm(Struct data), include: 
%            - options.Ifsort: input 0(default) or 1, if 1, swap rows and
%                           colums to ensure the diagonal elements maximum
%            - options.p: use p-norm
%            - options.X0: Intial value of algorithm
%            - options.TolFun: olerance of equation
%            - options.Maxiter: Maximum number of iterations
%            - options.Display: If display the final iterations number
%            - options.omega: omega is the relaxation factor of 'SOR' method
%            - options.PlotFcns: Draw the error value of each iteration 
%                                during the execution of the algorithm
%
% Output:
%   Xstar    - Iterative solutions of equations
%   Gvrho    - Spectral radius of Iterative matrix
%   X_iter   - Iterative solutions in every iteration
%   errorFun - Error of linear equations in each iteration
%
% Usage:
%   [ ] = ILEQS(Ab) uses the default settings soluting linear equations
%   [ ] = ILEQS(Ab, method) uses the input method to iterate linear equatioins
%   [ ] = ILEQS(Ab, [],options) uses the options' settings soluting linear equatioins
%   Xstar = ILEQS( ... ) returns the iterative solution of linear equations
%   [~, Gvrho] = ILEQS( ... ) returns the spectral radius of iterative matrix
%   [~, ~, X_iter] = ILEQS( ... ) returns the solutions in each iterations
%   [~, ~, ~, error] = ILEQS( ... ) returns the error of linear equations
%       in each iterations

% Author  : ZH.Yuan
% Update  : 2021/10/20 (First Version: 2021/10/19)
% Email   : 937341489@qq.com (If any suggestions or questions)
