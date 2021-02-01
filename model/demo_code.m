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
Nw = 5; % # of weekly series
Nm = 10; % # of monthly series
Nq = 1; % # of quarterly series
Nr = 2; % # of factors
Np = 1; % # of lags in factor VAR
Nt = 3600; % # of observations (daily frequency)

% params
lam_d = [0.6 + 0.1 * randn(Nd, 1), -0.4 + 0.1 * randn(Nd, 1)]; 
lam_w = [0.6 + 0.1 * randn(Nw, 1), -0.4 + 0.1 * randn(Nw, 1)]; 
lam_m = [0.6 + 0.1 * randn(Nm, 1), -0.4 + 0.1 * randn(Nm, 1)];
lam_q = [0.6 + 0.1 * randn(Nq, 1), -0.4 + 0.1 * randn(Nq, 1)];

sig_d = 0.5 + unifrnd(0.0, 0.2, Nd, 1);
sig_w = 0.5 + unifrnd(0.0, 0.2, Nw, 1);
sig_m = 0.5 + unifrnd(0.0, 0.2, Nm, 1);
sig_q = 0.5 + unifrnd(0.0, 0.2, Nq, 1);

phi_f = 0.3 * eye(Nr);
Omeg = diag(unifrnd(0.8,1.2, Nr, 1));

% initialize mats for storage
y_d = NaN(Nd, Nt);
y_w = NaN(Nw, Nt);
y_m = NaN(Nm, Nt);
y_q = NaN(Nq, Nt);
f = NaN(Nr, Nt);
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
for t = 1:Nt
    if t == 1
        f(:,t) = mvnrnd(zeros(Nr, 1), Omeg);
        f_w(:, t) = f(:, t);
        f_m(:, t) = f(:, t);
        f_q(:, t) = f(:, t);
    else
        % daily factor and data
        f(:,t) = phi_f * f(:, t-1) + mvnrnd(zeros(Nr, 1), Omeg)';
        y_d(:, t) = lam_d * f(:, t) + sig_d .* randn(Nd, 1);
        
        % weekly factor and data
        if Xi_w(t) == 0
            f_w(:, t) = f(:, t);
        else
            f_w(:, t) = f_w(:, t-1) + f(:, t);
        end
        y_w(:, t) = lam_w * f_w(:, t) + sig_w .* randn(Nw, 1);
        
        % monthly factor and data
        if Xi_m(t) == 0
            f_m(:, t) = f(:, t);
        else
            f_m(:, t) = f_m(:, t-1) + f(:, t);
        end
        y_m(:, t) = lam_m * f_m(:, t) + sig_m .* randn(Nm, 1);
        
        % quarterly factor and data
        if Xi_q(t) == 0
            f_q(:, t) = f(:, t);
        else
            f_q(:, t) = f_q(:, t-1) + f(:, t);
        end          
        y_q(:, t) = lam_q * f_q(:, t) + sig_q .* randn(Nq, 1);    
    end    
end

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

Ns = 0; % # of states
if Nd > 0; Ns = Ns+Nr;end 
if Nw > 0; Ns = Ns+Nr;end 
if Nm > 0; Ns = Ns+Nr;end 
if Nq > 0; Ns = Ns+Nr;end 

% T, => time-varying
T0 = [eye(Nr), zeros(Nr, Ns-Nr);
      repmat(-eye(Nr), (Ns-Nr)/Nr, 1) eye(Ns-Nr)];
  
iT0 = T0 \ eye(Ns); 

T = NaN(Ns, Ns, Nt);
for t = 1:Nt
    Ttmp = [phi_f zeros(Nr, Ns-Nr)];
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
Z = [Z_d Z_w Z_m Z_q]; 

H = diag([sig_d; sig_w; sig_m; sig_q]);

%-------------------------------------------------------------------------%
% run Kalman filter
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

