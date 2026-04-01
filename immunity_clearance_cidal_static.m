function immunity_clearance_cidal_static()

    % ==============================================================
    % IMMUNITY RANGE
    % ==============================================================
    I_values = logspace(1,7,40);   % from 1e1 to 1e7 immune cells
    tspan = [0 200];
    Y0 = [1e6; 1e6; 0]; %(S,T,R)
    clearance_threshold = 1;

    % ==============================================================
    % DEFINE DRUG PROFILES 
    % ==============================================================

    cidal_drugs = { ...
        struct('name','Amphotericin B','ED50',2), ...
        struct('name','Caspofungin','ED50',4) ...
    };

    static_drugs = { ...
        struct('name','Fluconazole','ED50',0,'KT',64,'K',2), ...
        struct('name','Voriconazole','ED50',0,'KT',8,'K',0.25), ...
        struct('name','Itraconazole','ED50',0,'KT',0.5,'K',0.25), ...
        struct('name','Posaconazole','ED50',0,'KT',0.5,'K',0.25) ...
    };

    %% ==============================================================
    %  CIDAL DRUG PLOT
    % ==============================================================

    figure; hold on;
    colors = lines(3);   % S,T,R colors
    h = gobjects(3,1);   % legend handles

    for d = 1:numel(cidal_drugs)
        drug = cidal_drugs{d};

        clear_S = nan(size(I_values));
        clear_T = nan(size(I_values));
        clear_R = nan(size(I_values));

        for j = 1:length(I_values)
            I = I_values(j);

            params = get_params_cidal(640, drug.ED50); 
            params.I = I; % sweep immunity

            [t,Y] = ode45(@(t,Y) AMR_cidal(t,Y,params), tspan, Y0);
            S = Y(:,1); T = Y(:,2); R = Y(:,3);

            idxS = find(S < clearance_threshold, 1);
            idxT = find(T < clearance_threshold, 1);
            idxR = find(R < clearance_threshold, 1);

            if ~isempty(idxS), clear_S(j) = t(idxS); end
            if ~isempty(idxT), clear_T(j) = t(idxT); end
            if ~isempty(idxR), clear_R(j) = t(idxR); end
        end

        % Plot one subplot per drug
        subplot(1, numel(cidal_drugs), d);
        hold on;
        set(gca,'XScale','log','FontSize', 14, 'FontWeight', 'bold');

        h(1) = plot(I_values, clear_S, 'o', 'Color', colors(1,:), 'MarkerSize', 10);
        h(2) = plot(I_values, clear_T, 's', 'Color', colors(2,:), 'MarkerSize', 10);
        h(3) = plot(I_values, clear_R, 'd', 'Color', colors(3,:), 'MarkerSize', 10);

       
        xticks([1e1 1e2 1e3 1e4 1e5 1e6 1e7]);
        xticklabels({'10^1','10^2','10^3','10^4','10^5','10^6','10^7'});
        ylim([0 200]);
       
        title([drug.name], 'FontSize', 15, 'FontWeight', 'bold');
        grid on; legend show;
    end
    % Shared labels
   %han = axes('Position',[0 0 1 1],'Visible','off');
   han = axes('Visible','off');
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    %hx=xlabel(han,'Immunity (cells)','FontSize',18,'FontWeight','bold');
    %hy=ylabel(han,'Clearance time (hrs)','FontSize',18,'FontWeight','bold');
    xlabel(han,'Immunity (cells)','FontSize',18,'FontWeight','bold');
    ylabel(han,'Clearance time (hrs)','FontSize',18,'FontWeight','bold');

    % Move labels
        %hx.Units = 'normalized';
        %hx.Position(2) = -5.0;

        %hy.Units = 'normalized';
        %hy.Position(1) = -5.0;

    % Single legend
    legend(h, {'Susceptible','Tolerant','Resistant'}, ...
       'Position',[0.88 0.5 0.1 0.15],'FontSize',14);


    %% ==============================================================
    %  STATIC DRUG PLOT
    % ==============================================================

    figure;

    for d = 1:numel(static_drugs)
        drug = static_drugs{d};

        clear_S = nan(size(I_values));
        clear_T = nan(size(I_values));
        clear_R = nan(size(I_values));

        for j = 1:length(I_values)
            I = I_values(j);

            params = get_params_static(640, drug.ED50, drug.K, drug.KT);
            params.I = I;

            [t,Y] = ode45(@(t,Y) AMR_static(t,Y,params), tspan, Y0);
            S = Y(:,1); T = Y(:,2); R = Y(:,3);

            idxS = find(S < clearance_threshold, 1);
            idxT = find(T < clearance_threshold, 1);
            idxR = find(R < clearance_threshold, 1);

            if ~isempty(idxS), clear_S(j) = t(idxS); end
            if ~isempty(idxT), clear_T(j) = t(idxT); end
            if ~isempty(idxR), clear_R(j) = t(idxR); end
        end

        % Plot one subplot per drug
        h = gobjects(3,1);
        subplot(2,2,d);
        hold on; 
        set(gca,'XScale','log','FontSize', 14, 'FontWeight', 'bold');

        h(1) = plot(I_values, clear_S, 'o', 'Color', colors(1,:), 'MarkerSize', 10);
        h(2) = plot(I_values, clear_T, 's', 'Color', colors(2,:), 'MarkerSize', 10);
        h(3) = plot(I_values, clear_R, 'd', 'Color', colors(3,:), 'MarkerSize', 10);

        xticks([1e1 1e2 1e3 1e4 1e5 1e6 1e7]);
        xticklabels({'10^1','10^2','10^3','10^4','10^5','10^6','10^7'});

        ylim([0 200]);
       
        title([drug.name], 'FontSize', 14, 'FontWeight', 'bold');
        grid on; legend show;
        han = axes(gcf,'Visible','off');
        han.XLabel.Visible = 'on';
        han.YLabel.Visible = 'on';
        hx=xlabel(han,'Immunity (cells)','FontSize',18,'FontWeight','bold');
        hy=ylabel(han,'Clearance time (hrs)','FontSize',18,'FontWeight','bold');

        % Move labels
        hx.Units = 'normalized';
        hx.Position(2) = -5.0;

        hy.Units = 'normalized';
        hy.Position(1) = -5.0;

        legend(h, {'Susceptible','Tolerant','Resistant'}, ...
       'Position',[0.85 0.45 0.12 0.15],'FontSize',14);

    end

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
%CIDAL
function params = get_params_cidal(D, ED50)
    % Growth rates
    params.rS = 0.3;
    params.rT = 0.5;
    params.rR = 0.8;

    % Population carrying capacity
    params.KP = 1e8;

    % Immunity killing rate
    params.deltaI = 1e-6;

    % Drug killing rates (cidal)
    params.deltaD = 1.0;   % killing of S
    params.deltaT = 0.1;   % killing of T

    % Drug-related parameters
    params.ED50 = ED50;    % drug ED50
    params.D = 640;          % drug dose

    % Transition rates
    params.kappa_TS = 0.97;
    params.kappa_ST = 0.0053;
    params.kappa_RT = 1e-8;

end

%STATIC
function params = get_params_static(D, ED50, K, KT)
    % Growth rates for static drugs
    params.rS = 0;        % static drug blocks S growth
    params.rT = 0.3;      
    params.rR = 0.8;

    % Population carrying capacity
    params.KP = 1e8;

    % Immunity killing rate
    params.deltaI = 1e-6;

    % No killing by static drugs
    params.deltaD = 0;
    params.deltaT = 0;

    % Drug parameters
    params.ED50 = ED50;
    params.K = K;
    params.KT = KT;
    params.D = 640;

    % Transition rates
    params.kappa_TS = 0.97;
    params.kappa_ST = 0.0053;
    params.kappa_RT = 1e-8;

end
