%% AMB_Mathematical_Model_and_PID_Demo_FINAL.m

clear; clc; close all;

%% 1. PARAMETERS TRANSCRIBED FROM THE REFERENCE TEST-RIG TABLE
p.m       = 46.0;       % equivalent rotor mass [kg]
p.ki      = 129.7;      % force/current coefficient [N/A]
p.ksTable = -2.4e6;     % signed displacement stiffness in the paper [N/m]
p.gap     = 250e-6;     % radial protecting-bearing half air gap [m]
p.fsw     = 20e3;       % converter switching frequency [Hz]
p.rpmMin  = 0;
p.rpmMax  = 4000;

% The reference paper writes:
%   m*x_ddot = -ksTable*x + ki*i + Fd, with ksTable < 0.
% Therefore -ksTable is POSITIVE and destabilizing.
% V6 stores the positive magnitude of that destabilizing coefficient.
p.kNegMag = -p.ksTable; % positive destabilizing stiffness magnitude [N/m]

fprintf('\n=== SINGLE-AXIS AMB MATHEMATICAL MODEL ===\n');
fprintf('m = %.3f kg\n',p.m);
fprintf('ki = %.3f N/A\n',p.ki);
fprintf('Signed ks from paper = %.4e N/m\n',p.ksTable);
fprintf('Destabilizing magnitude used by V6 = %.4e N/m\n',p.kNegMag);
fprintf('Protecting-bearing half air gap = %.1f um\n',p.gap*1e6);
fprintf('Converter switching frequency = %.1f kHz\n',p.fsw/1e3);

%% 2. CONTINUOUS-TIME MODEL
% Paper sign convention:
%   m*x_ddot = -ksTable*x + ki*i + Fd
%
% Since ksTable = -2.4e6 N/m:
%   m*x_ddot = +2.4e6*x + 129.7*i + Fd
%
% Define state z = [x; v], where v = dx/dt.
% z_dot = A*z + B*i + E*Fd.

A = [0, 1;
     p.kNegMag/p.m, 0];
B = [0;
     p.ki/p.m];
E_F = [0;
       1/p.m];       % maps disturbance FORCE [N] to acceleration
E_acc = [0;
         1];          % maps disturbance ACCELERATION [m/s^2] to acceleration
C = [1, 0];
D = 0;

fprintf('\nContinuous-time state-space model:\n');
disp('A ='); disp(A);
disp('B ='); disp(B);
disp('E_F (force input) ='); disp(E_F);
disp('E_acc (acceleration input) ='); disp(E_acc);
disp('C ='); disp(C);

% Transfer function from incremental control current to displacement:
% X(s)/I(s) = ki / (m*s^2 - p.kNegMag)
fprintf('Plant transfer function: X(s)/I(s) = %.4g / (%.4g*s^2 - %.4g)\n', ...
        p.ki,p.m,p.kNegMag);

%% 3. OPEN-LOOP POLES: VERIFY INSTABILITY
poles = eig(A);
fprintf('Open-loop poles [rad/s]:\n');
disp(poles);
fprintf('Positive pole = %.4f rad/s; open-loop plant is unstable.\n',max(real(poles)));

% Analytic result: s = +/-sqrt(kNegMag/m)
poleMagnitude = sqrt(p.kNegMag/p.m);
fprintf('Analytic poles: +%.4f and -%.4f rad/s\n',poleMagnitude,-poleMagnitude);

%% 4. SAMPLE THE PLANT USING EXACT ZERO-ORDER HOLD
% V6 uses dt = 0.0002 s -> 5 kHz controller/simulation update rate.
% This is distinct from the 20 kHz converter switching frequency.
dt = 2e-4;
T  = 0.25;
N  = floor(T/dt)+1;
t  = (0:N-1)'*dt;

% Augmented matrix exponential gives exact ZOH matrices for constant input
% over each sample interval. Ed maps disturbance acceleration a_d = F_d/m:
% [z(k+1)] = Ad*z(k) + Bd*i(k) + Ed*d_acc(k)
M = [A, B, E_acc;
     zeros(2,4)];
Md = expm(M*dt);
Ad = Md(1:2,1:2);
Bd = Md(1:2,3);
Ed = Md(1:2,4);

fprintf('\nDiscrete-time exact-ZOH model at dt = %.4g s:\n',dt);
fprintf('Controller sample rate = %.1f kHz\n',1/dt/1e3);
disp('Ad ='); disp(Ad);
disp('Bd ='); disp(Bd);
disp('Ed (maps disturbance acceleration a_d = F_d/m) ='); disp(Ed);

%% 5. PID CONTROLLER CONNECTED TO THE PLANT
% This is the baseline PID specified in V6. V6 subsequently searches for
% better gains and also evaluates 2-state and 3-state Kalman-filter variants.
Kp = 2.0e5;       % [A/m]
Ki = 1.0e6;       % [A/(m*s)]
Kd = 500;         % [A*s/m]
tauD = 2.5e-4;    % derivative-filter time constant [s]
iLimit = 15;      % current saturation [A]
reference = 100e-6; % desired displacement [m]

% Discrete derivative low-pass filter:
% dFilt(k) = alpha*dFilt(k-1) + (1-alpha)*(y(k)-y(k-1))/dt
alpha = tauD/(tauD+dt);

% PID equation used in V6:
% i_raw = Kp*e + Ki*integral(e) + Kd*de
% e = reference - measured_position
% de = -filtered_derivative(measured_position)
%
% The sign is consistent with positive ki: positive current produces
% positive incremental force, while a displacement above the reference
% makes the error/current command negative.

x = zeros(N,1);
v = zeros(N,1);
iCmd = zeros(N,1);
y = zeros(N,1);
err = zeros(N,1);
dFilt = 0;
yPrev = 0;
integralError = 0;
currentPrev = 0;

rng(2026,'twister');
sigma = 2e-6; % nominal measurement noise standard deviation [m]
sensorNoise = sigma*randn(N,1);

for k = 1:N-1
    y(k) = x(k) + sensorNoise(k);
    err(k) = reference - y(k);

    if k > 1
        dFilt = alpha*dFilt + (1-alpha)*(y(k)-yPrev)/dt;
    end
    derivativeError = -dFilt;
    integralError = integralError + err(k)*dt;

    iRaw = Kp*err(k) + Ki*integralError + Kd*derivativeError;
    iSat = min(max(iRaw,-iLimit),iLimit);

    % Simple anti-windup correction consistent with the V6 implementation.
    if iSat ~= iRaw && sign(iRaw) == sign(err(k))
        integralError = integralError - 0.5*err(k)*dt;
    end
    iCmd(k) = iSat;

    % Example rotating-unbalance disturbance, represented as acceleration.
    % mu*r is an assumed equivalent unbalance parameter, not supplied in
    % the reference Table V. V6 uses 1e-4 kg*m and 3000 rpm by default.
    muR = 1e-4; rpm = 3000;
    omega = 2*pi*rpm/60;
    disturbanceAcceleration = muR*omega^2*sin(omega*t(k))/p.m;

    zNext = Ad*[x(k);v(k)] + Bd*iCmd(k) + Ed*disturbanceAcceleration;
    x(k+1) = zNext(1);
    v(k+1) = zNext(2);

    yPrev = y(k);
    currentPrev = iCmd(k); %#ok<NASGU>
end
y(N) = x(N) + sensorNoise(N);
err(N) = reference - y(N);

%% 6. PERFORMANCE METRICS FOR THIS DEMONSTRATION
% These are illustrative baseline-PID results from this standalone demo.
% The full V6 script computes seed-averaged metrics for all four controllers.
x_um = x*1e6;
ref_um = reference*1e6;
rmse_um = sqrt(mean((x_um-ref_um).^2));
overshoot_pct = max(0,(max(x_um)-ref_um)/ref_um*100);
rmsCurrent_A = sqrt(mean(iCmd.^2));

fprintf('\n=== BASELINE PID DEMONSTRATION (one noise seed) ===\n');
fprintf('Kp = %.4g A/m, Ki = %.4g A/(m*s), Kd = %.4g A*s/m\n',Kp,Ki,Kd);
fprintf('Current limit = +/- %.1f A\n',iLimit);
fprintf('RMSE = %.4f um\n',rmse_um);
fprintf('Overshoot = %.3f %%\n',overshoot_pct);
fprintf('RMS current = %.4f A\n',rmsCurrent_A);
fprintf('These values are from this demo run only; use V6 exported tables for paper results.\n');

%% 7. PLOTS TO SHOW THE MATHEMATICAL MODEL AND FEEDBACK IMPLEMENTATION
figure('Color','w','Name','AMB mathematical model and PID demonstration');
subplot(3,1,1);
plot(t*1e3,x_um,'LineWidth',1.3); hold on;
plot(t*1e3,ref_um*ones(size(t)),'k--','LineWidth',1.0);
grid on; ylabel('Position (\mum)');
title('Single-axis AMB: baseline PID closed-loop response');
legend('Rotor displacement','Reference','Location','best');

subplot(3,1,2);
plot(t*1e3,iCmd,'LineWidth',1.1); hold on;
yline(iLimit,'k--'); yline(-iLimit,'k--');
grid on; ylabel('Current (A)');
title('PID output after current saturation');

subplot(3,1,3);
plot(t*1e3,err*1e6,'LineWidth',1.1);
grid on; xlabel('Time (ms)'); ylabel('Error (\mum)');
title('Measured position error');

% Open-loop natural response of the unstable model, with zero input.
% Small initial displacement is used only to illustrate the unstable pole.
% Stop the illustrative open-loop trajectory when it reaches the listed
% protecting-bearing half-gap. Beyond this point, the linearized free-motion
% model is not a physically credible description of the real rig.
tOL = (0:dt:0.05)';
zOL = zeros(2,numel(tOL));
zOL(:,1) = [1e-6;0];
nOL = 1;
for k = 1:numel(tOL)-1
    zOL(:,k+1) = Ad*zOL(:,k); % zero current and zero disturbance
    nOL = k+1;
    if abs(zOL(1,k+1)) >= p.gap
        break;
    end
end
tOL = tOL(1:nOL);
zOL = zOL(:,1:nOL);
figure('Color','w','Name','Open-loop AMB instability');
plot(tOL*1e3,zOL(1,:)*1e6,'LineWidth',1.3); hold on;
yline(p.gap*1e6,'r--','Protecting-bearing half-gap');
grid on; xlabel('Time (ms)'); ylabel('Displacement (\mum)');
title('Open-loop instability (illustration stops at half-gap)');

%% 8. MODEL SCOPE AND ASSUMPTIONS
fprintf('\n=== MODEL SCOPE ===\n');
fprintf(['Single-axis linearized model; gravity and cross-axis coupling omitted.\n' ...
         'Amplifier, sensor, computation delay, saturation physics and magnetic nonlinearities\n' ...
         'are not fully modeled. Current saturation and measurement noise are included.\n' ...
         'The air gap is reported as a rig parameter but does not enter the linearized equation\n' ...
         'directly; it would matter in a nonlinear magnetic-force model or safety constraints.\n' ...
         'The reference rig speed range is 0-4000 rpm; the V6 speed sweep should be aligned\n' ...
         'to that range for a strict reproduction of the reference rig.\n\n']);
