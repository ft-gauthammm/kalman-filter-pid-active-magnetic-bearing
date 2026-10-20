%% AMB_PID_KF_Project_V7_AUDITED.m

clear; clc; close all;
FAST = false;

%% 1. PARAMETERS
p.m  = 46.0;        % rotor mass [kg]
p.ki = 129.7;       % force/current gain [N/A]
p.ks = +2.40e6;     % negative-stiffness coefficient, POSITIVE sign => unstable [N/m]
assert(p.ks > 0,'ks must be positive: an AMB is open-loop unstable');
p.a = p.ks/p.m;     % x_ddot = a*x + b*i + d
p.b = p.ki/p.m;

fprintf('\n=== AMB model: m=%.1f kg, ki=%.1f N/A, ks=%.3g N/m ===\n',p.m,p.ki,p.ks);
fprintf('Open-loop poles: +%.2f and -%.2f rad/s (real; one is unstable)\n\n',sqrt(p.a),sqrt(p.a));

cfg.plant = struct('m',p.m,'ki',p.ki,'ks',p.ks);
cfg.dt     = 2e-4;          % controller / simulation step [s]
cfg.ref    = 100e-6;        % 100 um displacement step [m]
cfg.sigma  = 2e-6;          % nominal sensor noise std [m]
cfg.ilim   = 15;            % current limit [A]
cfg.mur    = 1e-4;          % unbalance mu*r [kg m]
cfg.rpm    = 3000;
cfg.Tfin   = 3.0;           % validation horizon [s]
cfg.q1     = 1e-3;          % small velocity process noise used in the 3-state KF
cfg.W      = struct('RMSE',0.25,'OS',0.15,'ST',0.15,'JIT',0.20,'CUR',0.25);

% Baseline PID by textbook rules (all values lie inside the search grid):
%   Kp = 2e5 A/m, Kd for PD damping ratio ~1, Ki = Kp/Ti with Ti = 0.2 s,
%   derivative filter tau_d = Td/N with Td = Kd/Kp and N = 10.
cfg.base = struct('Kp',2.0e5,'Ki',1.0e6,'Kd',500,'tau',2.5e-4,'qi',1);

mk = @(fm,fk,fs) struct('m',p.m*fm,'ki',p.ki*fk,'ks',p.ks*fs);
if ~FAST
    cfg.Topt     = 1.5;                       % optimisation horizon [s]
    cfg.seedsOpt = [101 102];
    cfg.seedsFin = 1:20;
    cfg.seedsSw  = 1:8;
    cfg.seedsMis = 1:5;
    cfg.refine   = true;
    cfg.ens      = [mk(1,1,1), mk(1,0.9,1), mk(1,1,1.1)];   % tuning ensemble (model error!)
    cfg.grid     = struct('Kp',[1.2e5 2.0e5 3.0e5],'Ki',[3e5 1e6 3e6], ...
                          'Kd',[300 500 800],'tau',[2.5e-4 1e-3 3e-3]);
    cfg.qGrid    = struct('kf2',[1e-1 1e1 1e3],'kf3',[1e0 1e2 1e4]);
else
    cfg.Topt     = 1.0;
    cfg.seedsOpt = 101;
    cfg.seedsFin = 1:4;
    cfg.seedsSw  = 1:3;
    cfg.seedsMis = 1:2;
    cfg.refine   = false;
    cfg.ens      = mk(1,1,1);
    cfg.grid     = struct('Kp',[1.2e5 3.0e5],'Ki',[3e5 3e6],'Kd',[300 800],'tau',[2.5e-4 3e-3]);
    cfg.qGrid    = struct('kf2',[1e1 1e3],'kf3',[1e2 1e4]);
end
nE = numel(cfg.ens);

%% 2. BASELINE ON THE TUNING ENSEMBLE (normalisers for the cost)
ctrlBase = makeCtrl('pid',cfg.base,[]);
baseMet = cell(1,nE);
for j = 1:nE
    baseMet{j} = runMetrics(ctrlBase,cfg.ens(j),cfg,cfg.Topt,cfg.sigma,cfg.seedsOpt);
    if isnan(baseMet{j}.ST), baseMet{j}.ST = cfg.Topt; end
end

%% 3. OPTIMISATION (grid + local refinement)
fprintf('=== OPTIMISATION ===\n');
SR.pid = searchController('pid',cfg,baseMet);
printBest('Optimized PID',SR.pid,'pid');
SR.kf2 = searchController('kf2',cfg,baseMet);
printBest('KF-PID (2-state KF)',SR.kf2,'kf2');
SR.kf3 = searchController('kf3',cfg,baseMet);
printBest('KF-PID (3-state KF)',SR.kf3,'kf3');

names = {'Baseline PID';'Optimized PID';'KF-PID (2-state)';'KF-PID (3-state)'};
ctrls = {ctrlBase, makeCtrl('pid',SR.pid.best,[]), ...
         makeCtrl('kf2',SR.kf2.best,SR.kf2.kfLib), makeCtrl('kf3',SR.kf3.best,SR.kf3.kfLib)};
nC = numel(ctrls);

%% 4. VALIDATION (nominal plant, independent seeds, mean +/- std over seeds)
fprintf('\n=== VALIDATION: nominal plant, %.1f um noise, %d seeds, %.1f s ===\n', ...
    cfg.sigma*1e6,numel(cfg.seedsFin),cfg.Tfin);
% Frequency grid stays below Nyquist (pi/dt); 0.9*pi/dt is a safety margin.
wGrid = logspace(0.5,log10(0.9*pi/cfg.dt),160);
M = cell(1,nC); XS = cell(1,nC); US = cell(1,nC);
Sens = zeros(nC,numel(wGrid)); Comp = zeros(nC,numel(wGrid));
for i = 1:nC
    [M{i},XS{i},US{i}] = runMetrics(ctrls{i},cfg.plant,cfg,cfg.Tfin,cfg.sigma,cfg.seedsFin);
    [Sens(i,:),Comp(i,:)] = loopResp(ctrls{i},cfg.plant,cfg,wGrid);
    if M{i}.div, fprintf('WARNING: %s diverged on the nominal plant\n',names{i}); end
    fprintf('%-18s RMSE %.2f+/-%.2f um | OS %5.1f %% | Ts %s | jitter %.2f um | I_ss %.2f+/-%.2f A | Ms %.2f | compl.pk %.3f um/N\n', ...
        names{i},M{i}.RMSE,M{i}.RMSEstd,M{i}.OS,fmtST(M{i}.ST),M{i}.JIT,M{i}.CUR,M{i}.CURstd, ...
        max(Sens(i,:)),max(Comp(i,:)));
end
fprintf('\nKF-PID (2-state) vs Optimized PID: effort %+.1f %%, jitter %+.1f %%, overshoot %+.1f %%\n', ...
    100*(M{3}.CUR/M{2}.CUR-1),100*(M{3}.JIT/M{2}.JIT-1),100*(M{3}.OS/M{2}.OS-1));
fprintf('KF-PID (3-state) vs Optimized PID: effort %+.1f %%, jitter %+.1f %%, overshoot %+.1f %%\n', ...
    100*(M{4}.CUR/M{2}.CUR-1),100*(M{4}.JIT/M{2}.JIT-1),100*(M{4}.OS/M{2}.OS-1));

%% 5. SENSOR-NOISE SWEEP (gains fixed; KF R matched to the noise level, PID cannot adapt)
noiseUm = [0.5 1 2 5 10 20];
kinds = {'pid','kf2','kf3'};
sweepJit = nan(numel(noiseUm),3); sweepCur = nan(numel(noiseUm),3); sweepRm = nan(numel(noiseUm),3);
fprintf('\n=== NOISE SWEEP (steady-state jitter um / current A) ===\n');
for n = 1:numel(noiseUm)
    sg = noiseUm(n)*1e-6;
    ln = sprintf('noise %5.1f um |',noiseUm(n));
    for k = 1:3
        if k == 1
            c = ctrls{2};
        else
            pp = SR.(kinds{k}).best; pp.qi = 1;
            lib = designKF(kinds{k},SR.(kinds{k}).qBest,sg,cfg);
            c = makeCtrl(kinds{k},pp,lib);
        end
        m = runMetrics(c,cfg.plant,cfg,2.0,sg,cfg.seedsSw);
        sweepJit(n,k) = m.JIT; sweepCur(n,k) = m.CUR; sweepRm(n,k) = m.RMSE;
        if m.div, ln = [ln sprintf(' %s DIV |',kinds{k})]; %#ok<AGROW>
        else,     ln = [ln sprintf(' %s %5.2f/%5.2f |',kinds{k},m.JIT,m.CUR)]; end %#ok<AGROW>
    end
    disp(ln);
end

%% 6. MODEL-MISMATCH STUDY (controller + KF keep the NOMINAL model)
misNames = {'nominal','ki -20%','ki -10%','ki +10%','ki +20%','ks -20%','ks -10%','ks +10%','ks +20%','m -20%','m +20%'};
misF = [1 1 1; 1 .8 1; 1 .9 1; 1 1.1 1; 1 1.2 1; 1 1 .8; 1 1 .9; 1 1 1.1; 1 1 1.2; .8 1 1; 1.2 1 1];
misRm = nan(numel(misNames),3); misCur = nan(numel(misNames),3);
fprintf('\n=== MODEL MISMATCH (RMSE um / steady-state current A) ===\n');
for n = 1:numel(misNames)
    pl = mk(misF(n,1),misF(n,2),misF(n,3));
    ln = sprintf('%-8s |',misNames{n});
    for k = 1:3
        m = runMetrics(ctrls{k+1},pl,cfg,2.0,cfg.sigma,cfg.seedsMis);
        misRm(n,k) = m.RMSE; misCur(n,k) = m.CUR;
        if m.div, ln = [ln sprintf(' %s DIVERGED |',kinds{k})]; %#ok<AGROW>
        else,     ln = [ln sprintf(' %s %6.2f/%5.2f |',kinds{k},m.RMSE,m.CUR)]; end %#ok<AGROW>
    end
    disp(ln);
end

%% 7. ROTOR-SPEED STUDY (unbalance force ~ omega^2)
% Reference rig range is 0–4000 rpm; 5000 rpm would be an extrapolation.
rpmList = [1000 2000 3000 4000];
rpmJit = nan(numel(rpmList),3);
for n = 1:numel(rpmList)
    c2 = cfg; c2.rpm = rpmList(n);
    for k = 1:3
        m = runMetrics(ctrls{k+1},cfg.plant,c2,2.0,cfg.sigma,cfg.seedsMis);
        rpmJit(n,k) = m.JIT;
    end
end

%% 8. FIGURES
t = (0:size(XS{1},2)-1)*cfg.dt;
col = lines(4);
figure('Name','Step response','Color','w');
subplot(2,1,1); hold on;
for i = 1:nC, plot(t*1e3,mean(XS{i},1)*1e6,'Color',col(i,:),'LineWidth',1.3); end
plot([0 300],[1 1]*cfg.ref*1e6,'k--'); xlim([0 300]); grid on;
ylabel('displacement, seed average (\mum)'); legend([names; {'reference'}],'Location','best');
title('100 \mum step: baseline vs optimized PID vs KF-PID');
subplot(2,1,2); hold on;
for i = 1:nC, plot(t*1e3,US{i}(1,:),'Color',col(i,:)); end
xlim([0 300]); grid on; xlabel('time (ms)'); ylabel('current, 1 seed (A)');

[X,U,V,Xh,Vh,Y] = simulate(ctrls{3},cfg.plant,cfg,0.4,cfg.sigma,1); %#ok<ASGLU>
tt = (0:numel(X)-1)*cfg.dt;
figure('Name','Kalman estimation','Color','w');
subplot(2,1,1); hold on;
plot(tt*1e3,Y*1e6,'Color',[.7 .7 .7]); plot(tt*1e3,X*1e6,'LineWidth',1.3); plot(tt*1e3,Xh*1e6,'--','LineWidth',1.3);
grid on; ylabel('position (\mum)'); legend('noisy sensor','true','Kalman estimate','Location','best');
title(sprintf('Kalman filter (2-state), sensor noise %.0f \\mum',cfg.sigma*1e6));
subplot(2,1,2); hold on;
plot(tt*1e3,V*1e3,'LineWidth',1.3); plot(tt*1e3,Vh*1e3,'--','LineWidth',1.3);
grid on; xlabel('time (ms)'); ylabel('velocity (mm/s)'); legend('true','estimate','Location','best');

figure('Name','Sensitivity and compliance','Color','w');
subplot(1,2,1); hold on;
for i = 1:nC, semilogx(wGrid/2/pi,20*log10(Sens(i,:)),'Color',col(i,:),'LineWidth',1.2); end
set(gca,'XScale','log'); grid on; xlabel('frequency (Hz)'); ylabel('|S| (dB)'); title('Sensitivity'); legend(names,'Location','best');
subplot(1,2,2); hold on;
for i = 1:nC, loglog(wGrid/2/pi,Comp(i,:),'Color',col(i,:),'LineWidth',1.2); end
set(gca,'XScale','log','YScale','log'); grid on; xlabel('frequency (Hz)'); ylabel('|x/F| (\mum/N)'); title('Dynamic compliance');

figure('Name','Noise sweep','Color','w');
subplot(1,2,1); loglog(noiseUm,sweepJit,'-o','LineWidth',1.2); grid on;
xlabel('sensor noise (\mum)'); ylabel('rotor jitter (\mum rms)'); legend('Optimized PID','KF-PID 2-state','KF-PID 3-state','Location','best');
subplot(1,2,2); loglog(noiseUm,sweepCur,'-o','LineWidth',1.2); grid on;
xlabel('sensor noise (\mum)'); ylabel('steady-state current (A rms)');

figure('Name','Model mismatch','Color','w');
bar(misRm); grid on; set(gca,'XTick',1:numel(misNames),'XTickLabel',misNames); xtickangle(45);
ylabel('RMSE (\mum)  [missing bar = diverged]'); legend('Optimized PID','KF-PID 2-state','KF-PID 3-state','Location','best');
title('Robustness to plant parameter error (controller keeps nominal model)');

figure('Name','Search trade-off','Color','w'); hold on;
scatter(SR.pid.allCur,SR.pid.allJit,18,'filled'); scatter(SR.kf2.allCur,SR.kf2.allJit,18,'filled'); scatter(SR.kf3.allCur,SR.kf3.allJit,18,'filled');
grid on; xlabel('steady-state current (A rms)'); ylabel('rotor jitter (\mum rms)');
legend('PID candidates','KF-PID 2-state','KF-PID 3-state','Location','best'); title('Every evaluated candidate (ensemble average)');

figure('Name','Rotor speed','Color','w');
plot(rpmList,rpmJit,'-o','LineWidth',1.2); grid on; xlabel('rotor speed (rpm)'); ylabel('rotor jitter (\mum rms)');
legend('Optimized PID','KF-PID 2-state','KF-PID 3-state','Location','best'); title('Unbalance disturbance vs speed');

%% 9. EXPORT
Kp = zeros(nC,1); Ki = Kp; Kd = Kp; tau_ms = Kp; Qv = Kp;
for i = 1:nC
    Kp(i) = ctrls{i}.Kp; Ki(i) = ctrls{i}.Ki; Kd(i) = ctrls{i}.Kd;
    if i <= 2, tau_ms(i) = ctrls{i}.tau*1e3; Qv(i) = NaN; else, tau_ms(i) = NaN; end
end
Qv(3) = SR.kf2.qBest; Qv(4) = SR.kf3.qBest;
RMSE_um = cellfun(@(m) m.RMSE,M)'; RMSE_std = cellfun(@(m) m.RMSEstd,M)';
Overshoot_pct = cellfun(@(m) m.OS,M)'; Settling_ms = cellfun(@(m) 1e3*m.ST,M)';
Jitter_um = cellfun(@(m) m.JIT,M)'; Current_A = cellfun(@(m) m.CUR,M)'; Current_std = cellfun(@(m) m.CURstd,M)';
Ms = max(Sens,[],2); Compliance_pk = max(Comp,[],2);
comparisonTable = table(names,Kp,Ki,Kd,tau_ms,Qv,RMSE_um,RMSE_std,Overshoot_pct,Settling_ms, ...
    Jitter_um,Current_A,Current_std,Ms,Compliance_pk);
disp(comparisonTable);

noiseTable = array2table([noiseUm(:) sweepJit sweepCur],'VariableNames', ...
    {'Noise_um','Jit_PID','Jit_KF2','Jit_KF3','Cur_PID','Cur_KF2','Cur_KF3'});
misTable = array2table([misRm misCur],'VariableNames', ...
    {'RMSE_PID','RMSE_KF2','RMSE_KF3','Cur_PID','Cur_KF2','Cur_KF3'});
misTable = [table(misNames(:),'VariableNames',{'Case'}) misTable];
rpmTable = array2table([rpmList(:) rpmJit],'VariableNames',{'RPM','Jit_PID','Jit_KF2','Jit_KF3'});

writetable(comparisonTable,'AMB_V7_Comparison.csv');
writetable(noiseTable,'AMB_V7_NoiseSweep.csv');
writetable(misTable,'AMB_V7_Mismatch.csv');
writetable(rpmTable,'AMB_V7_RotorSpeed.csv');
save('AMB_V7_results.mat','cfg','SR','comparisonTable','noiseTable','misTable','rpmTable','Sens','Comp','wGrid');
fprintf('\nSaved: AMB_V7_Comparison.csv, AMB_V7_NoiseSweep.csv, AMB_V7_Mismatch.csv, AMB_V7_RotorSpeed.csv, AMB_V7_results.mat\n');
fprintf('DONE.\n');

%% ===================================================================
%  LOCAL FUNCTIONS
%  ===================================================================

function [Ad,Bd,Ed] = plantDisc(m,ki,ks,dt)
% Exact ZOH of x'=v, v'=(ks/m)x+(ki/m)u+d  ->  x+ = Ad x + Bd u + Ed d
a = ks/m; b = ki/m;
Mx = [0 1 0 0; a 0 b 1; 0 0 0 0; 0 0 0 0];
E = expm(Mx*dt);
Ad = E(1:2,1:2); Bd = E(1:2,3); Ed = E(1:2,4);
end

function [A,B] = kfModel3(m,ki,ks,dt)
% KF model with disturbance state d (acceleration, random walk): states [x v d], input u
a = ks/m; b = ki/m;
Mx = [0 1 0 0; a 0 1 b; 0 0 0 0; 0 0 0 0];
E = expm(Mx*dt);
A = E(1:3,1:3); B = E(1:3,4);
end

function L = dlqeGain(A,Q,R)
% Steady-state Kalman (update) gain by Riccati iteration; measurement = first state.
n = size(A,1); C = zeros(1,n); C(1) = 1;
P = eye(n)*1e-9; K = zeros(n,1);
for it = 1:200000
    Pp = A*P*A' + Q;
    K  = Pp*C'/(C*Pp*C' + R);
    Pn = (eye(n)-K*C)*Pp;
    if max(abs(Pn(:)-P(:))) < 1e-12*max(abs(Pn(:)))
        P = Pn; break;
    end
    P = Pn;
end
L = K;
end

function lib = designKF(kind,qList,sigma,cfg)
% Kalman gains (always designed on the NOMINAL model) for every Q in qList.
dt = cfg.dt; pl = cfg.plant;
lib = struct('A',{},'B',{},'L',{},'q',{});
for i = 1:numel(qList)
    q = qList(i);
    if strcmp(kind,'kf2')
        [A,B] = plantDisc(pl.m,pl.ki,pl.ks,dt);
        Q = diag([0.25*dt^4*q, dt^2*q]);
    else
        [A,B] = kfModel3(pl.m,pl.ki,pl.ks,dt);
        Q = diag([0.25*dt^4*cfg.q1, dt^2*cfg.q1, q*dt]);
    end
    lib(i).A = A; lib(i).B = B; lib(i).L = dlqeGain(A,Q,sigma^2); lib(i).q = q;
end
end

function c = makeCtrl(kind,pp,lib)
c.kind = kind; c.Kp = pp.Kp; c.Ki = pp.Ki; c.Kd = pp.Kd;
if strcmp(kind,'pid')
    c.tau = pp.tau;
else
    c.A = lib(pp.qi).A; c.B = lib(pp.qi).B; c.L = lib(pp.qi).L;
end
end

function [X,U,V,Xh,Vh,Y] = simulate(ctrl,plant,cfg,T,sigma,seed)
% Closed loop: noisy position y -> (PID | KF+PID) -> current u (saturated) -> AMB.
dt = cfg.dt;
[Ad,Bd,Ed] = plantDisc(plant.m,plant.ki,plant.ks,dt);
a00 = Ad(1,1); a01 = Ad(1,2); a10 = Ad(2,1); a11 = Ad(2,2);
b0 = Bd(1); b1 = Bd(2); e0 = Ed(1); e1 = Ed(2);
N = floor(T/dt)+1;
rng(seed,'twister');
nz = sigma*randn(N,1);
w = 2*pi*cfg.rpm/60;
dist = cfg.mur*w^2*sin(w*dt*(0:N-1)')/plant.m;     % unbalance acceleration
X = zeros(N,1); U = zeros(N,1); V = zeros(N,1); Xh = zeros(N,1); Vh = zeros(N,1); Y = zeros(N,1);
x = 0; v = 0; I = 0; up = 0;
isPID = strcmp(ctrl.kind,'pid');
Kp = ctrl.Kp; Ki = ctrl.Ki; Kd = ctrl.Kd;
if isPID
    alpha = ctrl.tau/(ctrl.tau+dt); dF = 0; yPrev = 0;
else
    A = ctrl.A; B = ctrl.B; L = ctrl.L; z = zeros(size(A,1),1);
end
for k = 1:N
    y = x + nz(k);
    if isPID
        if k > 1, dF = alpha*dF + (1-alpha)*(y-yPrev)/dt; end
        e = cfg.ref - y; de = -dF;
        xh = y; vh = dF;
    else
        zp = A*z + B*up;                    % predict with the previous current
        z  = zp + L*(y - zp(1));            % measurement update
        e = cfg.ref - z(1); de = -z(2);
        xh = z(1); vh = z(2);
    end
    I = I + e*dt;
    u = Kp*e + Ki*I + Kd*de;
    us = min(max(u,-cfg.ilim),cfg.ilim);
    if us ~= u && sign(u) == sign(e)        % simple back-calculation anti-windup
        I = I - 0.5*e*dt;
    end
    u = us;
    X(k) = x; V(k) = v; U(k) = u; Y(k) = y; Xh(k) = xh; Vh(k) = vh;
    xn = a00*x + a01*v + b0*u + e0*dist(k);
    v  = a10*x + a11*v + b1*u + e1*dist(k);
    x  = xn;
    up = u;
    if isPID, yPrev = y; end
end
end

function [m,Xs,Us] = runMetrics(ctrl,plant,cfg,T,sigma,seeds)
N = floor(T/cfg.dt)+1;
Xs = zeros(numel(seeds),N); Us = zeros(numel(seeds),N);
for i = 1:numel(seeds)
    [X,U] = simulate(ctrl,plant,cfg,T,sigma,seeds(i));
    Xs(i,:) = X'; Us(i,:) = U';
end
m = computeMetrics(Xs,Us,cfg);
end

function m = computeMetrics(Xs,Us,cfg)
% Overshoot and settling from the SEED-AVERAGED response; RMSE/jitter/effort per seed.
m = struct('div',false,'RMSE',NaN,'RMSEstd',NaN,'OS',NaN,'ST',NaN,'JIT',NaN,'JITstd',NaN,'CUR',NaN,'CURstd',NaN);
xu = Xs*1e6; ref = cfg.ref*1e6;
if ~all(isfinite(xu(:))) || max(abs(xu(:))) > 1e5
    m.div = true; return;
end
N = size(xu,2); i0 = floor(0.6*N)+1;            % steady-state window = last 40 %
xa = mean(xu,1);
r = sqrt(mean((xu-ref).^2,2));
m.RMSE = mean(r); m.RMSEstd = std(r);
m.OS = max(0,max(xa)-ref)/ref*100;
xs = movmean(xa,100);                            % 20 ms centred smoothing
out = find(abs(xs-ref) > 0.05*ref);              % 5 % band
if isempty(out),          m.ST = 0;
elseif out(end) >= N,     m.ST = NaN;
else,                     m.ST = out(end)*cfg.dt;
end
jj = std(xu(:,i0:end),0,2);
cc = sqrt(mean(Us(:,i0:end).^2,2));
m.JIT = mean(jj); m.JITstd = std(jj);
m.CUR = mean(cc); m.CURstd = std(cc);
end

function J = costOf(m,b,W)
if m.div, J = 1e3; return; end
if isnan(m.ST), st = 3; else, st = min(m.ST/b.ST,3); end
J = W.RMSE*m.RMSE/b.RMSE + W.OS*(m.OS+1)/(b.OS+1) + W.ST*st ...
  + W.JIT*m.JIT/b.JIT + W.CUR*m.CUR/b.CUR;
end

function [J,cur,jit] = evalCand(kind,pp,kfLib,cfg,baseMet)
ctrl = makeCtrl(kind,pp,kfLib);
nE = numel(cfg.ens);
J = 0; cur = 0; jit = 0;
for j = 1:nE
    m = runMetrics(ctrl,cfg.ens(j),cfg,cfg.Topt,cfg.sigma,cfg.seedsOpt);
    J = J + costOf(m,baseMet{j},cfg.W);
    cur = cur + m.CUR; jit = jit + m.JIT;        % NaN if diverged (not plotted)
end
J = J/nE; cur = cur/nE; jit = jit/nE;
end

function SRk = searchController(kind,cfg,baseMet)
isPID = strcmp(kind,'pid');
G = cfg.grid;
if isPID
    kfLib = []; nq = 1; tauList = G.tau; qList = NaN;
else
    qList = cfg.qGrid.(kind); nq = numel(qList);
    kfLib = designKF(kind,qList,cfg.sigma,cfg);
    tauList = 0;
end
allCost = []; allCur = []; allJit = [];
bestJ = inf; best = struct();
for kp = G.Kp
  for kI = G.Ki
    for kd = G.Kd
      for tau = tauList
        for qi = 1:nq
            pp = struct('Kp',kp,'Ki',kI,'Kd',kd,'tau',tau,'qi',qi);
            [J,cur,jit] = evalCand(kind,pp,kfLib,cfg,baseMet);
            allCost(end+1) = J; allCur(end+1) = cur; allJit(end+1) = jit; %#ok<AGROW>
            if J < bestJ, bestJ = J; best = pp; end
        end
      end
    end
  end
end
if cfg.refine
    F = [0.8 1.0 1.25];
    if isPID, tF = F; else, tF = 1; end
    b0 = best;
    for fa = F
      for fb = F
        for fc = F
          for fd = tF
            pp = b0;
            pp.Kp = b0.Kp*fa; pp.Ki = b0.Ki*fb; pp.Kd = b0.Kd*fc; pp.tau = b0.tau*fd;
            [J,cur,jit] = evalCand(kind,pp,kfLib,cfg,baseMet);
            allCost(end+1) = J; allCur(end+1) = cur; allJit(end+1) = jit; %#ok<AGROW>
            if J < bestJ, bestJ = J; best = pp; end
          end
        end
      end
    end
end
SRk.best = best; SRk.J = bestJ; SRk.nEval = numel(allCost);
SRk.allCost = allCost; SRk.allCur = allCur; SRk.allJit = allJit;
SRk.kfLib = kfLib; SRk.qList = qList;
if isPID, SRk.qBest = NaN; else, SRk.qBest = qList(best.qi); end
end

function printBest(label,SRk,kind)
b = SRk.best;
if strcmp(kind,'pid')
    fprintf('%-20s Kp=%.4g Ki=%.4g Kd=%.4g tau_d=%.3g ms | cost %.4f | %d evals\n',label,b.Kp,b.Ki,b.Kd,b.tau*1e3,SRk.J,SRk.nEval);
else
    fprintf('%-20s Kp=%.4g Ki=%.4g Kd=%.4g Q=%.3g | cost %.4f | %d evals\n',label,b.Kp,b.Ki,b.Kd,SRk.qBest,SRk.J,SRk.nEval);
end
end

function [S,Cp] = loopResp(ctrl,plant,cfg,w)
% Frequency response of the sampled loop (nominal plant).
% Controller written as xi_k = Phi*xi_{k-1} + Gam*y_k, u_k = Cc*xi_k  (reference = 0).
dt = cfg.dt;
[Ad,Bd,Ed] = plantDisc(plant.m,plant.ki,plant.ks,dt);
if strcmp(ctrl.kind,'pid')
    al  = ctrl.tau/(ctrl.tau+dt);
    Phi = [al 0 -(1-al)/dt; 0 1 0; 0 0 0];
    Gam = [(1-al)/dt; -dt; 1];
    Cc  = [-ctrl.Kd, ctrl.Ki, -ctrl.Kp];
else
    A = ctrl.A; B = ctrl.B; L = ctrl.L; n = size(A,1);
    e1 = zeros(n,1); e1(1) = 1;
    Mx = eye(n) - L*e1';
    Cz = zeros(1,n); Cz(1) = -ctrl.Kp; Cz(2) = -ctrl.Kd;
    CI = ctrl.Ki;
    Pz  = Mx*A + (Mx*B)*Cz;
    PzI = (Mx*B)*CI;
    Phi = zeros(n+1); Phi(1:n,1:n) = Pz; Phi(1:n,n+1) = PzI;
    Phi(n+1,1:n) = -dt*Pz(1,:); Phi(n+1,n+1) = 1 - dt*PzI(1);
    Gam = [L; -dt*L(1)];
    Cc  = [Cz, CI];
end
S = zeros(size(w)); Cp = zeros(size(w));
nx = numel(Gam); Ix = eye(nx); I2 = eye(2);
for i = 1:numel(w)
    z  = exp(1j*w(i)*dt);
    K  = Cc*((Ix - Phi/z)\Gam);
    Mz = z*I2 - Ad;
    G  = [1 0]*(Mz\Bd);
    Gd = [1 0]*(Mz\Ed);
    Lp = -G*K;
    Sv = 1/(1+Lp);
    S(i)  = abs(Sv);
    Cp(i) = abs(Gd*Sv)/plant.m*1e6;             % |x/F| in um/N
end
end

function s = fmtST(x)
if isnan(x), s = 'not settled'; else, s = sprintf('%.0f ms',1e3*x); end
end
