%final simulation

function dYdt = AMR_with_Immunity(t, Y, params)
%[t, Y] = ode45(@(t,Y) AMR_with_Immunity(t, Y, params), tspan, Y0);

% Unpack state variables
S = Y(1); %susceptible
T = Y(2); %tolerant
R = Y(3); %resistant

Ptot = S + T + R;

% Unpack parameters
rS = params.rS;
rT = params.rT;
rR = params.rR;

KP = params.KP;
%K = params.K; %static drug
%KT = params.KT; %static drug

deltaI = params.deltaI;
deltaD = params.deltaD;
deltaT = params.deltaT;
ED50 = params.ED50;
D = params.D;  % drug concentration

kappa_TS = params.kappa_TS;
kappa_ST = params.kappa_ST;
kappa_RT = params.kappa_RT;

I = params.I;

% ODEs
%for cidal drugs
dS = rS * S * (1 - Ptot / KP) ...
    - deltaI * I * S ...
    - deltaD * (D / (D + ED50)) * S ...
    - kappa_TS * S ...
    + kappa_ST * T;

%for static drugs
%dS = rS * (K / (D + K)) * S * (1 - Ptot / KP) ...
    %- deltaI * I * S ...
    %- deltaD * (D / (D + ED50)) * S ...
    %- kappa_TS * S ...
    %+ kappa_ST * T;

%for cidal drugs
dT = rT * T * (1 - Ptot / KP) ...
    - deltaI * I * T ...
    - deltaT * (D / (D + ED50)) * T ...
    + kappa_TS * S ...
    - kappa_ST * T ...
    - kappa_RT * T;

%for static drugs
%dT = rT * (KT / (D + KT)) * T * (1 - Ptot / KP) ...
    %- deltaI * I * T ...
    %- deltaT * (D / (D + ED50)) * T ...
    %+ kappa_TS * S ...
    %- kappa_ST * T ...
    %- kappa_RT * T;

dR = rR * R * (1 - Ptot / KP) ...
    - deltaI * I * R ...
    + kappa_RT * T; ...


% Output derivative
dYdt = [dS; dT; dR];

end

%params.rS = 0; %static drugs
params.rS = 0.8; %cidal drugs
%params.rT = 0.3; %static drugs
params.rT = 0.5; %cidal drugs
params.rR = 0.3;
params.KP = 1e8;

params.deltaI = 1e-6;

params.deltaD = 1.0; %cidal drugs
params.deltaT = 0.1; %cidal drugs
%params.deltaD = 0; %static drugs
%params.deltaT = 0; %static drugs
params.D = 0;
params.ED50 = 1e-8;
%params.K = 0.25; %static drugs
%params.KT = 0.5; %static drugs

params.kappa_TS =0.0053;
params.kappa_ST = 0.0053;
params.kappa_RT = 1e-8;

%params.lambda = 0.5;
%params.KI = 1e7;
%params.deltaI = 0.01;
params.I = 1e6;

% Initial conditions: [S0, T0, R0, I0]
Y0 = [1e7; 1e7; 0]; 
tspan = [0 140];

[t, Y] = ode45(@(t,Y) AMR_with_Immunity(t, Y, params), tspan, Y0);

% Plot results
plot(t, Y, 'LineWidth', 2);
legend('Susceptible', 'Tolerant', 'Resistant', ...
       'FontSize', 15, 'FontWeight', 'bold');

xlabel('Time (hrs)', 'FontSize', 18, 'FontWeight', 'bold');
ylabel('Concentration (cells/mL)', 'FontSize', 18, 'FontWeight', 'bold');

%title('Cidal: Amphotericin B', ...
      %'FontSize', 18, 'FontWeight', 'bold');

set(gca, 'FontSize', 16, 'FontWeight', 'bold');


