function crossover_time_plot()
    % ==============================================================
    % Range of drug concentrations (µg/mL)
    % ==============================================================
    D_values = 2:10:640;  

    % ==============================================================
    % DEFINE DRUG PROFILES
    % ==============================================================
    % ----- CIDAL DRUGS -----
    cidal_drugs = { ...
        struct('name','Amphotericin B','ED50',2), ...
        struct('name','Caspofungin','ED50',4) ...
    };

    % ----- STATIC DRUGS -----
    static_drugs = { ...
        struct('name','Fluconazole','ED50',64,'KT',64,'K',2), ...
        struct('name','Voriconazole','ED50',8,'KT',8,'K',0.25), ...
        struct('name','Itraconazole','ED50',0.5,'KT',0.5,'K',0.25), ...
        struct('name','Posaconazole','ED50',0.5,'KT',0.5,'K',0.25) ...
    };

    % ==============================================================
    % Simulation settings
    % ==============================================================
    Y0 = [1e4; 1e4; 0]; %(S, T, R)
    tspan = [0 200];

    % ==============================================================
    % CIDAL DRUGS
    % ==============================================================
    figure;
    hold on;
    crossover_cidal_matrix = nan(length(D_values), numel(cidal_drugs));
    cidal_labels = strings(1, numel(cidal_drugs));

    for d = 1:numel(cidal_drugs)
        drug = cidal_drugs{d};
        crossover_cidal = nan(size(D_values));

        for i = 1:length(D_values)
            D = D_values(i);
            params = get_params_cidal(D, drug.ED50);
            [t, Y] = ode45(@(t,Y) AMR_cidal(t, Y, params), tspan, Y0);
            S = Y(:,1); T = Y(:,2); R = Y(:,3);
            idx = find(R > S & R > T, 1);
            if ~isempty(idx)
                crossover_cidal(i) = t(idx);
            end
        end

        crossover_cidal_matrix(:, d) = crossover_cidal;
        cidal_labels(d) = drug.name;
        plot(D_values, crossover_cidal, 'o', 'MarkerSize', 10, ...
        'DisplayName', drug.name);
    end

    xlabel('Drug concentration (µg/mL)', 'FontSize', 18, 'FontWeight', 'bold');
    ylabel('Crossover time (hrs)', 'FontSize', 18, 'FontWeight', 'bold');
    title('Cidal Drugs: Crossover Time vs. Drug Concentration', 'FontSize', 15, 'FontWeight', 'bold');
    legend('show', 'Location', 'northeast', 'FontSize', 15, 'FontWeight', 'bold');
    set(gca, 'FontSize', 15, 'FontWeight', 'bold');
    grid on;

   
     % Fix Y-axis range from 0 to 200 hours
    ylim([0 200]);

    % Save to CSV
    T_cidal = array2table([D_values(:), crossover_cidal_matrix], ...
        'VariableNames', ['Drug_Concentration', cidal_labels]);
    writetable(T_cidal, 'crossover_times_cidal.csv');
    fprintf('Saved cidal drug results to crossover_times_cidal.csv\n');

    % ==============================================================
    % STATIC DRUGS
    % ==============================================================
   figure; 
hold on;

crossover_static_matrix = nan(length(D_values), numel(static_drugs));
static_labels = strings(1, numel(static_drugs));

% Define line styles, markers, and colors
markers = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
colors = lines(numel(static_drugs));  % MATLAB's default distinct colors

for d = 1:numel(static_drugs)
    drug = static_drugs{d};
    crossover_static = nan(size(D_values));

    for i = 1:length(D_values)
        D = D_values(i);
        params = get_params_static(D, drug.ED50, drug.K, drug.KT);
        [t, Y] = ode45(@(t,Y) AMR_static(t, Y, params), tspan, Y0);
        S = Y(:,1); T = Y(:,2); R = Y(:,3);
        idx = find(R > S & R > T, 1);
        if ~isempty(idx)
            crossover_static(i) = t(idx);
        end
    end

    crossover_static_matrix(:, d) = crossover_static;
    static_labels(d) = drug.name;


    % Plot only markers (no lines)
    plot(D_values, crossover_static, 'LineStyle', 'none', ...
         'Marker', markers{mod(d-1, length(markers)) + 1}, ...
         'Color', colors(d,:), 'MarkerSize', 10, ...
         'DisplayName', drug.name);
end

xlabel('Drug concentration (µg/mL)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Crossover time (hrs)', 'FontSize', 18, 'FontWeight', 'bold');
title('Static Drugs: Crossover Time vs. Drug Concentration', 'FontSize', 15, 'FontWeight', 'bold');
legend('show', 'Location', 'northeast', 'FontSize', 15, 'FontWeight', 'bold');
set(gca, 'FontSize', 15, 'FontWeight', 'bold');
grid on;

     % Fix Y-axis range from 0 to 200 hours
    ylim([0 200]);

    % Save to CSV
    T_static = array2table([D_values(:), crossover_static_matrix], ...
        'VariableNames', ['Drug_Concentration', static_labels]);
    writetable(T_static, 'crossover_times_static.csv');
    fprintf('Saved static drug results to crossover_times_static.csv\n');
end

% ======================================================================
% ======================== ODE FUNCTIONS ================================
% ======================================================================

function dYdt = AMR_cidal(t, Y, params)
    S = Y(1); T = Y(2); R = Y(3);
    Ptot = S + T + R;

    rS = params.rS; rT = params.rT; rR = params.rR;
    KP = params.KP; deltaI = params.deltaI;
    deltaD = params.deltaD; deltaT = params.deltaT;
    ED50 = params.ED50; D = params.D;
    kappa_TS = params.kappa_TS; kappa_ST = params.kappa_ST;
    kappa_RT = params.kappa_RT; I = params.I;

    dS = rS*S*(1 - Ptot/KP) ...
        - deltaI*I*S ...
        - deltaD*(D/(D+ED50))*S ...
        - kappa_TS*S + kappa_ST*T;

    dT = rT*T*(1 - Ptot/KP) ...
        - deltaI*I*T ...
        - deltaT*(D/(D+ED50))*T ...
        + kappa_TS*S - kappa_ST*T - kappa_RT*T;

    dR = rR*R*(1 - Ptot/KP) ...
        - deltaI*I*R ...
        + kappa_RT*T;

    dYdt = [dS; dT; dR];
end

function dYdt = AMR_static(t, Y, params)
    S = Y(1); T = Y(2); R = Y(3);
    Ptot = S + T + R;

    rS = params.rS; rT = params.rT; rR = params.rR;
    KP = params.KP; K = params.K; KT = params.KT;
    deltaI = params.deltaI; deltaD = params.deltaD;
    deltaT = params.deltaT; ED50 = params.ED50; D = params.D;
    kappa_TS = params.kappa_TS; kappa_ST = params.kappa_ST;
    kappa_RT = params.kappa_RT; I = params.I;

    dS = rS*(K/(D+K))*S*(1 - Ptot/KP) ...
        - deltaI*I*S ...
        - deltaD*(D/(D+ED50))*S ...
        - kappa_TS*S + kappa_ST*T;

    dT = rT*(KT/(D+KT))*T*(1 - Ptot/KP) ...
        - deltaI*I*T ...
        - deltaT*(D/(D+ED50))*T ...
        + kappa_TS*S - kappa_ST*T - kappa_RT*T;

    dR = rR*R*(1 - Ptot/KP) ...
        - deltaI*I*R ...
        + kappa_RT*T;

    dYdt = [dS; dT; dR];
end

% ======================================================================
% ======================== PARAMETER FUNCTIONS ==========================
% ======================================================================

function params = get_params_cidal(D, ED50)
    params.rS = 0.3;
    params.rT = 0.5;
    params.rR = 0.8;
    params.KP = 1e8;
    params.deltaI = 1e-6;
    params.deltaD = 1.0;
    params.deltaT = 0.1;
    params.ED50 = ED50;
    params.kappa_TS = 0.97;
    params.kappa_ST = 0.0053;
    params.kappa_RT = 1e-8;
    params.I = 1e6;
    params.D = D;
end

function params = get_params_static(D, ED50, K, KT)
    params.rS = 0;
    params.rT = 0.3;
    params.rR = 0.8;
    params.KP = 1e8;
    params.deltaI = 1e-6;
    params.deltaD = 0;
    params.deltaT = 0;
    params.ED50 = ED50;
    params.K = K;
    params.KT = KT;
    params.kappa_TS = 0.97;
    params.kappa_ST = 0.0053;
    params.kappa_RT = 1e-8;
    params.I = 1e7;
    params.D = D;
end
