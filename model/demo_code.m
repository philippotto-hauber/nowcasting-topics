clear; close all; clc;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%- This code generates data from a multi-frequency factor model 
%- and estimates the coefficients via an EM algorithm. 
%- The model allows for series of daily, weekly, monthly and quarterly
%- series and the code is general enough so that any frequency may be
%- excluded. In such a case, the highest frequency factors are always 
%- ordered first. For details on the model and state space representation
%- see Modugno (2011). "Nowcasting inflation with high-frequency data" 
%- (https://www.ecb.europa.eu/pub/pdf/scpwps/ecbwp1324.pdf)
%- Note that the model implicitly assumes that the lower frequency
%- series are stock rather than flow variables. The latter can be implemented
%- following the (more complicated) approach in Banbura et al. (2011). 
%- "Nowcasting with Daily Data"(https://ideas.repec.org/p/red/sed012/555.html)
%-------------------------------------------------------------------------%
% generate data
%-------------------------------------------------------------------------%

% set-up
Nd = 100; % # of daily series  
Nw = 10; % # of weekly series
Nm = 10; % # of monthly series
Nq = 10; % # of quarterly series
Nr = 2; % # of factors
Np = 3; % # of lags in factor VAR
Np_eff = Np + 1; % # of lags of f in state vector (always needs to be one higher for covariance of factors in M-step!)
Nt = 3600; % # of observations (daily frequency)

% params
lam_d = [0.6 + 0.1 * randn(Nd, 1), -0.4 + 0.1 * randn(Nd, 1)]; 
lam_w = [0.6 + 0.1 * randn(Nw, 1), -0.4 + 0.1 * randn(Nw, 1)]; 
lam_m = [0.6 + 0.1 * randn(Nm, 1), -0.4 + 0.1 * randn(Nm, 1)];
lam_q = [0.6 + 0.1 * randn(Nq, 1), -0.4 + 0.1 * randn(Nq, 1)];

sig2_d = 0.2 + unifrnd(0.0, 0.2, Nd, 1);
sig2_w = 0.2 + unifrnd(0.0, 0.2, Nw, 1);
sig2_m = 0.2 + unifrnd(0.0, 0.2, Nm, 1);
sig2_q = 0.2 + unifrnd(0.0, 0.2, Nq, 1);

Phi = [0.3 * eye(Nr), -0.4 * eye(Nr), 0.1 * eye(Nr)];
Omeg = diag(unifrnd(0.8,1.2, Nr, 1));

% initialize mats for storage
y_d = NaN(Nd, Nt);
y_w = NaN(Nw, Nt);
y_m = NaN(Nm, Nt);
y_q = NaN(Nq, Nt);
F = NaN(Nr * Np, Nt); % companion form
f_w = NaN(Nr, Nt);
f_m = NaN(Nr, Nt);
f_q = NaN(Nr, Nt);

% indices for monthly and quarterly observations
N_d_w = 5; % # of days per week => working days
N_d_m = N_d_w * 4; % # of days per month
N_d_q = N_d_m * 3; % # of days per quarter

Xi_w = ones(Nt, 1); Xi_w(1:N_d_w:end) = 0; % equals 1 if start of new week
Xi_m = ones(Nt, 1); Xi_m(1:N_d_m:end) = 0; % equals 1 if start of new month
Xi_q = ones(Nt, 1); Xi_q(1:N_d_q:end) = 0; % equals 1 if start of new quarter

% loop over t
F(:, 1) = 0; 
for t = 1:Nt
    if t == 1
        F(1:Nr,t) = mvnrnd(zeros(Nr, 1), Omeg);
    else
        % daily factor 
        F(:,t) = [Phi; eye(Nr * (Np-1)) zeros(Nr * (Np-1), Nr)] * F(:, t-1) + [eye(Nr); zeros(Nr * (Np-1), Nr)] * mvnrnd(zeros(Nr, 1), Omeg)';
    end
    
    % daily data
    y_d(:, t) = lam_d * F(1:Nr, t) + sqrt(sig2_d) .* randn(Nd, 1);

    % weekly factor and data
    if Xi_w(t) == 0
        f_w(:, t) = F(1:Nr, t);
    else
        f_w(:, t) = f_w(:, t-1) + F(1:Nr, t);
    end
    y_w(:, t) = lam_w * f_w(:, t) + sqrt(sig2_w) .* randn(Nw, 1);

    % monthly factor and data
    if Xi_m(t) == 0
        f_m(:, t) = F(1:Nr, t);
    else
        f_m(:, t) = f_m(:, t-1) + F(1:Nr, t);
    end
    y_m(:, t) = lam_m * f_m(:, t) + sqrt(sig2_m) .* randn(Nm, 1);

    % quarterly factor and data
    if Xi_q(t) == 0
        f_q(:, t) = F(1:Nr, t);
    else
        f_q(:, t) = f_q(:, t-1) + F(1:Nr, t);
    end          
    y_q(:, t) = lam_q * f_q(:, t) + sqrt(sig2_q) .* randn(Nq, 1);    
end

% extract factor from companion form representation
f = F(1:Nr, :); 

% extract actually observed monthly and quarterly values
y_d_o = y_d; 
y_w_o = NaN(Nw, Nt);y_w_o(:, end:-N_d_w:1) = y_w(:, end:-N_d_w:1);
y_m_o = NaN(Nm, Nt);y_m_o(:, end:-N_d_m:1) = y_m(:, end:-N_d_m:1);
y_q_o = NaN(Nq, Nt);y_q_o(:, end:-N_d_q:1) = y_q(:, end:-N_d_q:1);

% plot factors and obs
figure; 
plot(f', 'Color', [0, 0.4470, 0.7410]); 
hold on; 
plot(f_w', 'Color', [0.8500, 0.3250, 0.0980]);
plot(f_m', 'Color', [0.9290, 0.6940, 0.1250]);
plot(f_q', 'Color', [0.4940, 0.1840, 0.5560]);

figure; 
plot(y_d', 'Color', [0, 0.4470, 0.7410]); 
hold on; 
plot(y_w', '-', 'Color', [0.8500, 0.3250, 0.0980, 0.2], 'LineWidth',0.5); % Color(4) = alpha!
plot(y_w_o', '-o', 'Color', [0.8500, 0.3250, 0.0980]);
plot(y_m', '-', 'Color', [0.9290, 0.6940, 0.1250, 0.2]);
plot(y_m_o', '-o','Color', [0.9290, 0.6940, 0.1250]);
plot(y_q', '-','Color', [0.4940, 0.1840, 0.5560, 0.2]);
plot(y_q_o', '-o','Color', [0.4940, 0.1840, 0.5560]);

%-------------------------------------------------------------------------%
% state space representation
%-------------------------------------------------------------------------%

Ns = Nr * Np_eff; % 
if Nw > 0 && Nd > 0; Ns = Ns+Nr;end 
if Nm > 0 && (Nw > 0 || Nd > 0) ; Ns = Ns+Nr;end 
if Nq > 0 && (Nm > 0 || Nw > 0 || Nd > 0); Ns = Ns+Nr;end 

% T, => time-varying
T0 = eye(Ns);
T0(Nr*Np_eff+1:Ns,1:Nr) = repmat(-eye(Nr), (Ns - Nr*Np_eff)/Nr, 1);
iT0 = T0 \ eye(Ns); 

T = NaN(Ns, Ns, Nt);
for t = 1:Nt
    Ttmp = [Phi zeros(Nr, Ns-size(Phi, 2));
            eye(Nr * (Np_eff-1)) zeros(Nr * (Np_eff-1), Ns - Nr * (Np_eff-1))];
    if Nw > 0 && Nd > 0 % if there are weekly series and its not the highest frequency. Otherwise its dynamics are governed by phi_f!
        Ttmp = [Ttmp; zeros(Nr, size(Ttmp, 1)) Xi_w(t) * eye(Nr) zeros(Nr, Ns - (size(Ttmp, 1) + Nr))];
    end
    if Nm > 0 && (Nw > 0 || Nd > 0) 
        Ttmp = [Ttmp; zeros(Nr, size(Ttmp, 1)) Xi_m(t) * eye(Nr) zeros(Nr, Ns - (size(Ttmp, 1) + Nr))];
    end
    if Nq > 0 && (Nm > 0 || Nw > 0 || Nd > 0)
        Ttmp = [Ttmp; zeros(Nr, size(Ttmp, 1)) Xi_q(t) * eye(Nr) zeros(Nr, Ns - (size(Ttmp, 1) + Nr))];
    end
    
    T(:, :, t) = iT0 * Ttmp;
end

% R, Q 
Q = zeros(Ns);
Q(1:Nr, 1:Nr) = Omeg; 
R = iT0; 

% Z, H 
Z_d = []; if Nd > 0; Z_d = [lam_d; zeros(Nw+Nm+Nq, Nr)]; end
Z_w = []; if Nw > 0; Z_w = [zeros(Nd, Nr); lam_w; zeros(Nm+Nq, Nr)]; end
Z_m = []; if Nm > 0; Z_m = [zeros(Nd+Nw, Nr); lam_m; zeros(Nq, Nr)]; end
Z_q = []; if Nq > 0; Z_q = [zeros(Nd+Nw+Nm, Nr); lam_q]; end
Z = [Z_d zeros(Nd+Nw+Nm+Nq, Nr * (Np_eff-1)) Z_w Z_m Z_q]; 

H = diag([sig2_d; sig2_w; sig2_m; sig2_q]);

%-------------------------------------------------------------------------%
% run Kalman smoother (E-step!)
%-------------------------------------------------------------------------%

s0 = zeros(Ns, 1); P0 = 10 * eye(Ns);
tic
[stT,PtT,LL] = f_KS_DK_logL([y_d_o; y_w_o; y_m_o; y_q_o],T,Z,H,R,Q,s0,P0);
toc
figure;
subplot(2,2,1)
plot(stT(1:Nr, :)', 'b')
hold on
plot(f', 'r')
title('f_d (sampled in blue)')
subplot(2,2,2)
plot(stT(Nr+1:2*Nr, :)', 'b')
hold on
plot(f_w', 'r')
title('f_w (sampled in blue)')
subplot(2,2,3)
plot(stT(2*Nr+1:3*Nr, :)', 'b')
hold on
plot(f_m', 'r')
title('f_m (sampled in blue)')
subplot(2,2,4)
plot(stT(Ns-Nr+1:end, :)', 'b')
hold on
plot(f_q', 'r')
title('f_q (sampled in blue)')

%-------------------------------------------------------------------------%
% estimate parameters (M-step)
%-------------------------------------------------------------------------%

% lam_d and sig_d
id_f = 1:Nr; id_x = 1:Nd;
[lam_d_hat, sig_d_hat] = f_sample_lam_sig(y_d_o, stT(id_f, :), PtT(id_f, id_f, :), sig2_d);

% lam_w and sig2_w
id_f = Nr * Np_eff+1:Nr*(Np_eff+1);
[lam_w_hat, sig2_w_hat] = f_sample_lam_sig(y_w_o, stT(id_f, :), PtT(id_f, id_f, :), sig2_w);

% lam_m and sig2_m
id_f = Nr * (Np_eff+1)+1:Nr * (Np_eff+2);
[lam_m_hat, sig2_m_hat] = f_sample_lam_sig(y_m_o, stT(id_f, :), PtT(id_f, id_f, :), sig2_m);

% lam_q and sig_q
id_f = Ns-Nr+1:Ns;
[lam_q_hat, sig2_q_hat] = f_sample_lam_sig(y_q_o, stT(id_f, :), PtT(id_f, id_f, :), sig2_q);

% phi_f and Omeg
id_f = 1:Nr;
id_f_lags = Nr+1:Nr * Np_eff;
Phi_hat = (stT(id_f, :)*stT(id_f_lags, :)' + sum(PtT(id_f,id_f_lags,:),3))/(stT(id_f_lags, :)*stT(id_f_lags, :)' + sum(PtT(id_f_lags,id_f_lags,:),3)) ; 
Omeg_hat = 1/Nt * ((stT(id_f, :)*stT(id_f, :)' + sum(PtT(id_f,id_f,:),3))  - Phi_hat * (stT(id_f, :)*stT(id_f_lags, :)' + sum(PtT(id_f,id_f_lags,:),3))') ;
 
% plot estimated versus true
figure;
subplot(5,2,1)
scatter(lam_d(:, 1), lam_d_hat(:, 1), 'b')
hold on;
scatter(lam_d(:, 2), lam_d_hat(:, 2), 'r')
ylim([-1 1])
xlim([-1 1])
refline(1, 0)
title('lam_d')
ylabel('estimate')
xlabel('actual')

subplot(5,2,2)
scatter(sig2_d, sig_d_hat)
ylim([0 1])
xlim([0 1])
refline(1, 0)
title('sig_d')
ylabel('estimate')
xlabel('actual')

subplot(5,2,3)
scatter(lam_w(:, 1), lam_w_hat(:, 1), 'b')
hold on;
scatter(lam_w(:, 2), lam_w_hat(:, 2), 'r')
ylim([-1 1])
xlim([-1 1])
refline(1, 0)
title('lam_w')
ylabel('estimate')
xlabel('actual')

subplot(5,2,4)
scatter(sig2_w, sig2_w_hat)
ylim([0 1])
xlim([0 1])
refline(1, 0)
title('sig_w')
ylabel('estimate')
xlabel('actual')

subplot(5,2,5)
scatter(lam_m(:, 1), lam_m_hat(:, 1), 'b')
hold on;
scatter(lam_m(:, 2), lam_m_hat(:, 2), 'r')
ylim([-1 1])
xlim([-1 1])
refline(1, 0)
title('lam_m')
ylabel('estimate')
xlabel('actual')

subplot(5,2,6)
scatter(sig2_m, sig2_m_hat)
ylim([0 1])
xlim([0 1])
refline(1, 0)
title('sig_m')
ylabel('estimate')
xlabel('actual')

subplot(5,2,7)
scatter(lam_q(:, 1), lam_q_hat(:, 1), 'b')
hold on;
scatter(lam_q(:, 2), lam_q_hat(:, 2), 'r')
ylim([-1 1])
xlim([-1 1])
refline(1, 0)
title('lam_q')
ylabel('estimate')
xlabel('actual')

subplot(5,2,8)
scatter(sig2_q, sig2_q_hat)
ylim([0 1])
xlim([0 1])
refline(1, 0)
title('sig_q')
ylabel('estimate')
xlabel('actual')

subplot(5,2,9)
scatter(Phi(:), Phi_hat(:))
ylim([-1.0 1.0])
xlim([-1.0 1.0])
refline(1, 0)
title('Phi')
ylabel('estimate')
xlabel('actual')
subplot(5,2,10)
scatter(Omeg(:), Omeg_hat(:))
ylim([-1.5 1.5])
xlim([-1.5 1.5])
refline(1, 0)
title('Omeg')
ylabel('estimate')
xlabel('actual')
