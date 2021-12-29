clc
clearvars
close all

times = 200;
Nall = [10 50 100 500];

NJ = zeros(1, times);
NG = zeros(1, times);
NS = zeros(1, times);
NA = zeros(1, times);

NJall = cell(1, numel(Nall));
NGall = cell(1, numel(Nall));
NSall = cell(1, numel(Nall));
NAall = cell(1, numel(Nall));

figure(1)

for iN = 1 : numel(Nall)
    for iter = 1 : times
        N = Nall(iN);
        Amethod = 4;

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
                A = N * 5 * sqrt(N) * diag(rand(1, N)) + rand(N);
                A = A' + A;
            case 5
                M = diag(rand(N, 1));
                Z = orth(rand(N, N));
                A = Z' * M * Z;
        end

        Xtrue = (1 : N)';
        b = A * Xtrue;
        options.Display = 0;
        options.ifsort = 1;
        options.PlotFcns = 'off';

        %% Compare and Plot
        if strcmp(options.PlotFcns, 'on')
            figure; hold on
        end

        % Jacob
        [XstarJ, GvrhoJ, X_iterJ, errorallJ] = ileqs([A b], 'Jacob', options);
        % G-S
        [XstarG, GvrhoG, X_iterG, errorallG] = ileqs([A b], 'G-S', options);
        % SOR
        [XstarS, GvrhoS, X_iterS, errorallS] = ileqs([A b], 'SOR', options);
        % Adaptive
        [XstarA, GvrhoA, X_iterA, errorallA] = ileqs([A b], [], options);

        if strcmp(options.PlotFcns, 'on')
            legend({'Jacob', 'G-S', 'SOR', 'Adaptive'})
            hold off
        end

        NJ(iter) = numel(errorallJ);
        NG(iter) = numel(errorallG);
        NS(iter) = numel(errorallS);
        NA(iter) = numel(errorallA);

    end

    figure(1)
    subplot(ceil(numel(Nall) / 2), 2, iN)
    fig = boxplot([NJ', NG', NS', NA'], {'Jacob', 'G-S', 'SOR', 'Adaptive'}, 'Widths',0.3);
    set(fig, 'Linewidth', 2);
    subtitle(['Random Simulation' num2str(N) ' Times'])
    ylabel('Stop Number of Iteration')

    NJall{iN} = NJ;
    NGall{iN} = NG;
    NSall{iN} = NS;
    NAall{iN} = NA;

end

figure(1)
sgtitle('Boxplot of Numbers of Iterations')
