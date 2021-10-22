function [Xstar, Gvrho, X_iter, errorFun] = ileqs(Ab, method, options)
% ILEQS - Iterative optimization solution of linear equations
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
%   [] = ILEQS(Ab) uses the default settings soluting linear equations
%   [] = ILEQS(Ab, method) uses the input method to iterate linear equatioins
%   [] = ILEQS(Ab, [],options) uses the options' settings soluting linear equatioins
%   Xstar = ILEQS( ... ) returns the iterative solution of linear equations
%   [~, Gvrho] = ILEQS( ... ) returns the spectral radius of iterative matrix
%   [~, ~, X_iter] = ILEQS( ... ) returns the solutions in each iterations
%   [~, ~, ~, error] = ILEQS( ... ) returns the error of linear equations
%       in each iterations

% Author  : ZH.Yuan
% Update  : 2021/10/20 (First Version: 2021/10/19)
% Email   : 937341489@qq.com (If any suggestions or questions)

% Check augmented matrix's format
[N, N1] = size(Ab);
if N1 - N ~= 1
    error('Format of A is Wrong, please check it!')
end
A = Ab(:, 1 : end - 1);
b = Ab(:, end);
if det(A) == 0
    error('The coefficient matrix is a singular matrix!')
end

% Set default value of value
if ~exist('method', 'var') || isempty(method)
    method = 'Adaptive';
end

% Set default value of options
if ~exist('options', 'var') || isempty(options)
    options.Ifsort = 1;
    options.p = 2;
    options.X0 = zeros(N, 1);
    options.TolFun = 1e-08;
    options.Maxiter = max([200 N]);
    options.Display = 0;
    options.omega = [];
    options.PlotFcns = 'off';
end

% Set default value of options.Ifsort
if ~isfield(options, 'Ifsort')
    options.Ifsort = 1;
elseif isempty(options.Ifsort)
    options.Ifsort = 1;
end

% Set default value of options.p
if ~isfield(options, 'p')
    options.p = 2;
elseif isempty(options.p)
    options.p = 2;
end

% Set default value of options.X0
if ~isfield(options, 'X0')
    options.X0 = zeros(N, 1);
elseif isempty(options.X0)
    options.X0 = zeros(N, 1);
end

% Set default value of options.TolFun
if ~isfield(options, 'TolFun')
    options.TolFun = 1e-08;
elseif isempty(options.TolFun)
    options.TolFun = 1e-08;
end

% Set default value of options.Maxiter
if ~isfield(options, 'Maxiter')
    options.Maxiter = max([200 N]);
elseif isempty(options.Maxiter)
    options.Maxiter = max([200 N]);
end

% Set default value of options.Display
if ~isfield(options, 'Display')
    options.Display = 0;
elseif isempty(options.Display)
    options.Display = 0;
end

% Set default value of options.omega
if ~isfield(options, 'omega')
    options.omega = [];
elseif isempty(options.omega)
    options.omega = [];
end

% Set default value of options.PlotFcns
if ~isfield(options, 'PlotFcns')
    options.PlotFcns = 'off';
elseif isempty(options.PlotFcns)
    options.PlotFcns = 'off';
end

% Sort the augmented matrix
if options.Ifsort
    [A, b, S] = smatrix(A, b);
else
    S = 1 : N;
end

% Iterate linear equations solution by input method
switch method
    case 'Jacob'
        [G, f] = JacobLE(A, b);
    case 'G-S'
        [G, f] = GSLE(A, b);
    case 'SOR'
        [G, f] = SORLE(A, b, options.omega);
end

if ~strcmp(method, 'Adaptive')
    Gvrho = max(abs(eig(G)));
    if Gvrho >= 1
        warning([method '''s L' num2str(options.p) '-norm is larger' ...
            ' than 1, follow soultions use the adaptive method'])
        method = 'Adaptive';
    end
end

if strcmp(method, 'Adaptive')
    [G, f, method] = AdaptiveLE(A, b, options.omega);
    Gvrho = max(abs(eig(G)));
end

iter = 0;
X0 = options.X0;

while iter < options.Maxiter
    iter = iter + 1;
    X_new = G * X0 + f;
    X_iter(:, iter) = X_new;
    errorFun(iter) = norm(A * X_new - b, options.p);
    if errorFun(iter) <= options.TolFun
        break
    end
    X0 = X_new;
end

if options.Display
    fprintf(['Algorithm-' mfilename ' stop at the %d-th iteration by ' ...
        method ' method.\n'], iter);
end

Xstar(S) = X_new;
X_iter(S, :) = X_iter;

if strcmp(options.PlotFcns, 'on')
    plot(errorFun, 'LineWidth', 2)
    title('Error of equations in each iteration')
    xlabel('Iteration number')
    ylabel('Error of equations')
end

end

%% Algorithm for replacing the main element of the coefficient matrix
function [A, b, Scol] = smatrix(A, b)
% replacing diagonal element of the coefficient matrix
%
% Input:
%   A  - Coefficient matrix
%   b  - Constant column vector
%
% Output:
%   AS - Coefficient matrix after sort
%   bS - Constant column vector after sort
%   Scol - Sort of column

% Author  : ZH.Yuan
% Update  : 2021/10/19 (First Version: 2021/10/19)
% Email   : 937341489@qq.com (If any suggestions or questions)

N = size(A, 1);
Scol = 1 : N;
k = 1;
ks = 1;
while k <= N
    Aiter = abs(A(k : N, k : N));
    Adet = Aiter;
    [~, Sind] = maxk(abs(Aiter(:)), ks);
    [Maxrow, Maxcol] = ind2sub([N + 1 - k, N + 1 - k], Sind(end));
    Adet(Maxrow, :) = [];
    Adet(:, Maxcol) = [];
    if det(Adet) == 0 
        ks = ks + 1;
    else
        ks = 1;
        A([k, Maxrow + k - 1], :) = A([Maxrow + k - 1, k], :);
        A(:, [k, Maxcol + k - 1]) = A(:, [Maxcol + k - 1, k]);
        b([k, Maxrow + k - 1], :) = b([Maxrow + k - 1, k], :);
        Scol(:, [k, Maxcol + k - 1]) = Scol(:, [Maxcol + k - 1, k]);
        k = k + 1;
    end
end

end

%% Iterative matrix and residual vector of Jacobian method
function [G, f] = JacobLE(A, b)
% Compute for iterative matrix and residual vector of Jacobian method
%
% Input:
%   A  - Coefficient matrix
%   b  - Constant column vector
%
% Output:
%   G  - Iterative matrix
%   f  - Residual vector

% Author  : ZH.Yuan
% Update  : 2021/10/19 (First Version: 2021/10/19)
% Email   : 937341489@qq.com (If any suggestions or questions)


G = eye(size(A, 1)) - diag(1./ diag(A)) * A;
f = diag(1./ diag(A)) * b;

end

%% Iterative matrix and residual vector of Gauss-Seidel method
function [G, f] = GSLE(A, b)
% Compute for iterative matrix and residual vector of Gauss-Seidel method
%
% Input:
%   A  - Coefficient matrix
%   b  - Constant column vector
%
% Output:
%   G  - Iterative matrix
%   f  - Residual vector

% Author  : ZH.Yuan
% Update  : 2021/10/19 (First Version: 2021/10/19)
% Email   : 937341489@qq.com (If any suggestions or questions)


U = - triu(A, 1);
L = - tril(A, -1);
D = diag(diag(A));

G = (D - L)^-1 * U;
f = (D - L)^-1 * b;

end

%% Iterative matrix and residual vector of Successive over relaxation method
function [G, f] = SORLE(A, b, omega)
% Compute for iterative matrix and residual vector of SOR method
%
% Input:
%   A  - Coefficient matrix
%   b  - Constant column vector
%   omega - Relaxation factor
%
% Output:
%   G  - Iterative matrix
%   f  - Residual vector

% Author  : ZH.Yuan
% Update  : 2021/10/19 (First Version: 2021/10/19)
% Email   : 937341489@qq.com (If any suggestions or questions)


U = - triu(A, 1);
L = - tril(A, -1);
D = diag(diag(A));
Gfun = @(omega) (D - omega * L)^-1 * ((1 - omega) * D + omega * U);
Gvrhofun = @(omega) max(abs(eig(Gfun(omega))));

if ~exist('omega', 'var') || isempty(omega)
    [Alower, Aupper] = bandwidth(A);
    if Alower + Aupper == 2
        pho = max(abs(eig(JacobLE(A, b))));
        omega = 2 / (1 + sqrt(1 - pho^2));
    else
        omega = fminbnd(Gvrhofun, 0, 2);
    end
end

G = Gfun(omega);
f = omega * (D - omega * L)^-1 * b;

end

%% Iterative matrix and residual vector of adaptive choose method
function [G, f, method] = AdaptiveLE(A, b, omega)
% Compute for iterative matrix and residual vector of adaptive choose method
%
% Input:
%   A  - Coefficient matrix
%   b  - Constant column vector
%   omega - Relaxation factor
%
% Output:
%   G  - Iterative matrix
%   f  - Residual vector
%   method  - The fastest convergence method

% Author  : ZH.Yuan
% Update  : 2021/10/19 (First Version: 2021/10/19)
% Email   : 937341489@qq.com (If any suggestions or questions)

if ~exist('omega', 'var')
    omega = [];
end

[G1, f1] = JacobLE(A, b);
[G2, f2] = GSLE(A, b);
[G3, f3] = SORLE(A, b, omega);

Gvrho1 = max(abs(eig(G1)));
Gvrho2 = max(abs(eig(G2)));
Gvrho3 = max(abs(eig(G3)));

[~, minGvrho] = min([Gvrho1, Gvrho2, Gvrho3]);

switch minGvrho
    case 1
        G = G1;
        f = f1;
        method = 'Jacob';
    case 2
        G = G2;
        f = f2;
        method = 'G-S';
    case 3
        G = G3;
        f = f3;
        method = 'SOR';
end

end