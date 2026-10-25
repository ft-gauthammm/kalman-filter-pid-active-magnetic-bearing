%% AMB_Classical_Control_Analysis.m
% Classical-control companion analysis for the single-axis AMB V6 model.
% IMPORTANT:
% (1) This script analyses the continuous-time plant and ideal continuous PID.
% (2) It does not claim that these plots alone prove stability of the sampled,
%     saturated, derivative-filtered KF-PID implementation.
% (3) Requires Control System Toolbox for tf, rlocus, bode, nyquist, margin.
clear; clc; close all;

%% Published/reference parameters and sign convention
m  = 46.0;             % kg
ki = 129.7;             % N/A
ksSigned = -2.40e6;     % N/m, signed value reported by reference
kn = abs(ksSigned);     % positive magnitude of destabilizing stiffness
fprintf('AMB classical analysis: m=%.3f kg, ki=%.3f N/A, ks(signed)=%.4g N/m\n',m,ki,ksSigned);
fprintf('Model convention: m*xdd = kn*x + ki*i + Fd, kn=|ks|=%.4g N/m\n',kn);

%% Continuous-time plant and open-loop poles
s = tf('s');
G = ki/(m*s^2 - kn);     % displacement/current plant
pPlant = pole(G);
disp('Plant X(s)/I(s) ='); G
disp('Open-loop plant poles (rad/s):'); disp(pPlant);
fprintf('Analytical poles = +/- %.6f rad/s\n',sqrt(kn/m));
fprintf('Positive pole means the uncontrolled linearized plant is unstable.\n');

%% State-space observability
A = [0 1; kn/m 0];
B = [0; ki/m];
C = [1 0];
D = 0;
Ob = obsv(A,C);
fprintf('\n2-state observability matrix:\n'); disp(Ob);
fprintf('rank(O) = %d of %d states\n',rank(Ob),size(A,1));
fprintf('For the 3-state augmented estimator, see the note below: the disturbance state must be checked with the exact discrete A3 used in V6.\n');

%% Ideal continuous-time PID characteristic polynomial and Routh-Hurwitz
% C(s)=Kp + Ki/s + Kd*s. Negative feedback with G(s)=ki/(m*s^2-kn).
% Characteristic polynomial: m*s^3 + ki*Kd*s^2 + (ki*Kp-kn)*s + ki*Ki = 0.
fprintf('\nRouth-Hurwitz conditions for ideal PID:\n');
fprintf('Kd>0, Ki>0, ki*Kp-kn>0, and ki*Kd*(ki*Kp-kn)>m*ki*Ki.\n');
fprintf('Minimum proportional gain: kn/ki = %.4f A/m\n',kn/ki);

% Gains from V6 comparison CSV; use ideal PID here, not derivative filter/KF.
names = {'Baseline PID','Optimized PID','KF-PID (2-state)','KF-PID (3-state)'};
Kp = [200000 160000 375000 375000];
Ki = [1000000 3750000 3750000 3750000];
Kd = [500 400 640 800];
for k = 1:numel(names)
    coeff = [m, ki*Kd(k), ki*Kp(k)-kn, ki*Ki(k)];
    % Routh first-column terms for cubic: a3, a2, (a2*a1-a3*a0)/a2, a0
    routhFirst = [coeff(1), coeff(2), ...
        (coeff(2)*coeff(3)-coeff(1)*coeff(4))/coeff(2), coeff(4)];
    fprintf('\n%s\n',names{k});
    fprintf('  Ideal PID characteristic coefficients: '); fprintf(' %.6g',coeff); fprintf('\n');
    fprintf('  Routh first column:'); fprintf(' %.6g',routhFirst); fprintf('\n');
    if all(routhFirst>0)
        fprintf('  Ideal continuous PID passes cubic Routh-Hurwitz test.\n');
    else
        fprintf('  Ideal continuous PID FAILS cubic Routh-Hurwitz test.\n');
    end
end

%% Root locus and frequency response for baseline and optimized ideal PID
ctrlIndices = [1 2];  % compare conventional baseline and optimized PID
figure('Name','AMB Root Locus','Color','w');
hold on;
for k = ctrlIndices
    Cpid = Kp(k) + Ki(k)/s + Kd(k)*s;
    % Root locus varies a scalar gain multiplying the fixed ideal PID shape.
    rlocus(G*Cpid);
end
grid on; title('Root locus: scalar multiplier on fixed ideal PID shapes');
xlabel('Real axis (rad/s)'); ylabel('Imaginary axis (rad/s)');
legend(names(ctrlIndices),'Location','best');

figure('Name','AMB Bode and Margins','Color','w');
for j = 1:numel(ctrlIndices)
    k = ctrlIndices(j);
    Cpid = Kp(k) + Ki(k)/s + Kd(k)*s;
    L = G*Cpid;
    subplot(2,1,j);
    margin(L); grid on;
    title(sprintf('%s: ideal continuous-time loop margins',names{k}));
    [Gm,Pm,Wcg,Wcp] = margin(L);
    fprintf('\n%s ideal-loop margins:\n',names{k});
    fprintf('  Gain margin = %.6g (%.3f dB), phase margin = %.3f deg\n',Gm,20*log10(Gm),Pm);
    fprintf('  Gain crossover = %.6g rad/s, phase crossover = %.6g rad/s\n',Wcp,Wcg);
end

figure('Name','AMB Nyquist','Color','w');
for j = 1:numel(ctrlIndices)
    k = ctrlIndices(j);
    Cpid = Kp(k) + Ki(k)/s + Kd(k)*s;
    subplot(1,numel(ctrlIndices),j);
    nyquist(G*Cpid); grid on; axis equal;
    title([names{k} ' ideal loop Nyquist']);
end

%% Closed-loop poles for the ideal continuous PID
fprintf('\nIdeal continuous PID closed-loop poles:\n');
for k = 1:numel(names)
    Cpid = Kp(k) + Ki(k)/s + Kd(k)*s;
    T = feedback(G*Cpid,1);
    fprintf('%s poles:\n',names{k}); disp(pole(T).');
end

%% Important interpretation notes
fprintf(['\nINTERPRETATION NOTES:\n' ...
 '1) These ideal-PID root-locus/Bode/Nyquist/margin results are not the exact V6 sampled loop.\n' ...
 '2) V6 uses dt=0.0002 s, a filtered derivative for classical PID, current saturation,\n' ...
 '   anti-windup, and Kalman estimator dynamics for KF-PID configurations.\n' ...
 '3) For publication, report exact sampled-loop margins only after deriving the return ratio\n' ...
 '   consistently from the implemented closed-loop architecture and verifying it independently.\n' ...
 '4) Nyquist encirclement direction/count must be interpreted with the number of RHP poles\n' ...
 '   and the precise loop-transfer sign convention; do not hard-code a counter-clockwise count.\n']);
