clc
clearvars
close all

N = 100;
Amethod = 3;

switch Amethod
    case 1
        A = N * eye(N) + rand(N);
    case 2
        A = 2 * eye(N) + rand(N);
        A = diag(diag(A, -1), -1) + diag(diag(A, 0), 0) + diag(diag(A, 1), 1);
    case 3
        A = N * eye(N) + rand(N);
        A = A' + A;
    case 4
        M = diag(rand(N, 1));
        Z = orth(rand(N, N));
        A = Z' * M * Z;
end

Xtrue = (1 : N)';
b = A * Xtrue;
options.Display = 1;
options.ifsort = 1;
options.PlotFcns = 'on';

%% Compare and Plot
figure; hold on
% Jacob
[XstarJ, GvrhoJ, X_iterJ, errorallJ] = ileqs([A b], 'Jacob', options);
% G-S
[XstarG, GvrhoG, X_iterG, errorallG] = ileqs([A b], 'G-S', options);
% SOR
[XstarS, GvrhoS, X_iterS, errorallS] = ileqs([A b], 'SOR', options);
% Adaptive
[XstarA, GvrhoA, X_iterA, errorallA] = ileqs([A b], [], options);
legend({'Jacob', 'G-S', 'SOR', 'Adaptive'})
hold off
