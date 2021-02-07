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
Nd = 50; % # of daily series  
Nw = 0; % # of weekly series
Nm = 0; % # of monthly series
Nq = 4; % # of quarterly series
Nq_flow = Nq / 2; % # of quarterly flow series
Nq_stock = Nq - Nq_flow; % # of quarterly stock series
ind_q_flow = [repelem(true, Nq_flow), repelem(false, Nq_stock)];
Nr = 1; % # of factors
Np = 1; % # of lags in factor VAR
Np_eff = Np + 1; % # of lags of f in state vector (always needs to be one higher for covariance of factors in M-step!)


% load dates corresponding to Jan 1st 1991-Dec 31 2018
tmp = importdata('dates_Xi_W_19912018.csv');
Nt = size(tmp.data, 1); % # of observations (daily frequency)
offset_nonnumvars = size(tmp.textdata, 2) - size(tmp.data, 2); % number of non-numeric vars (these are ordered first!) 
Xi_w =  tmp.data(:, find(contains(tmp.textdata(1,:), 'Xi_w')) - offset_nonnumvars); % equals 0 if start of new week
Xi_m =  tmp.data(:, find(contains(tmp.textdata(1,:), 'Xi_m')) - offset_nonnumvars); % equals 0 if start of new month
Xi_q =  tmp.data(:, find(contains(tmp.textdata(1,:), 'Xi_q')) - offset_nonnumvars); % equals 0 if start of new quarter
W_qd_c =  tmp.data(:, find(contains(tmp.textdata(1,:), 'W_qd_c')) - offset_nonnumvars); % equals 0 if start of new week
W_qd_p =  tmp.data(:, find(contains(tmp.textdata(1,:), 'W_qd_p')) - offset_nonnumvars); % equals 0 if start of new week
ind_plot = tmp.data(:, find(contains(tmp.textdata(1,:), 'ind_plot')) - offset_nonnumvars);
ind_plot = ind_plot(~isnan(ind_plot));
dates_plot = tmp.data(:, find(contains(tmp.textdata(1,:), 'dates_plot')) - offset_nonnumvars);
dates_plot = dates_plot(~isnan(dates_plot));
clearvars tmp offset_nonnumvars 

% params
switch Nr
    case 1
        lam_d = 0.6 + 0.1 * randn(Nd, 1); 
        lam_w = 0.6 + 0.1 * randn(Nw, 1); 
        lam_m = 0.6 + 0.1 * randn(Nm, 1);
        lam_q = 0.6 + 0.1 * randn(Nq, 1);
    case 2
        lam_d = [0.6 + 0.1 * randn(Nd, 1), -0.4 + 0.1 * randn(Nd, 1)]; 
        lam_w = [0.6 + 0.1 * randn(Nw, 1), -0.4 + 0.1 * randn(Nw, 1)]; 
        lam_m = [0.6 + 0.1 * randn(Nm, 1), -0.4 + 0.1 * randn(Nm, 1)];
        lam_q = [0.6 + 0.1 * randn(Nq, 1), -0.4 + 0.1 * randn(Nq, 1)];
    otherwise
        error('Nr has to be smaller than or equal to 2!')
end

sig2_d = 0.2 + unifrnd(0.0, 0.2, Nd, 1);
sig2_w = 0.2 + unifrnd(0.0, 0.2, Nw, 1);
sig2_m = 0.2 + unifrnd(0.0, 0.2, Nm, 1);
sig2_q = 0.2 + unifrnd(0.0, 0.2, Nq, 1);

switch Np
    case 1
        Phi = 0.5 * eye(Nr);
    case 2
        Phi = [0.5 * eye(Nr), -0.2 * eye(Nr)];
    case 3
        Phi = [0.3 * eye(Nr), -0.4 * eye(Nr), 0.1 * eye(Nr)];
    otherwise
        error('Np has to smaller than or equal to 3!')
end
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
f_q_c = NaN(Nr, Nt); % current quarter
f_q_p = NaN(Nr, Nt); % previous quarter

% loop over t
F(:, 1) = 0; % initialize f_0, f_-1, ..., f_-p+1
for t = 1:Nt
    % iterate factors forwars
    if t == 1 % initialize model!
        F(1:Nr,t) = mvnrnd(zeros(Nr, 1), Omeg);
        f_w(:, t) = F(1:Nr, t);
        f_m(:, t) = F(1:Nr, t);
        f_q(:, t) = F(1:Nr, t);
        f_q_c(:, t) = F(1:Nr, t);
        f_q_p(:,t) = zeros(Nr, 1); 
    else
        % daily factor 
        F(:,t) = [Phi; eye(Nr * (Np-1)) zeros(Nr * (Np-1), Nr)] * F(:, t-1) + [eye(Nr); zeros(Nr * (Np-1), Nr)] * mvnrnd(zeros(Nr, 1), Omeg)';
        
        % weekly factor 
        if Xi_w(t) == 0; f_w(:, t) = F(1:Nr, t); else;f_w(:, t) = f_w(:, t-1) + F(1:Nr, t); end
        
        % monthly factor 
        if Xi_m(t) == 0; f_m(:, t) = F(1:Nr, t); else; f_m(:, t) = f_m(:, t-1) + F(1:Nr, t); end
        
        % quarterly factor 
        if Xi_q(t) == 0 
            f_q(:, t) = F(1:Nr, t); 
            f_q_c(:, t) = f_q_p(:, t-1) + F(1:Nr, t);
            f_q_p(:, t) = zeros(Nr, 1); 
        else
            f_q(:, t) = f_q(:, t-1) + F(1:Nr, t);
            f_q_c(:, t) = f_q_c(:, t-1) + W_qd_c(t) * F(1:Nr, t);
            f_q_p(:, t) = f_q_p(:, t-1) + W_qd_p(t) * F(1:Nr, t); 
        end          
    end     
    
    % generate observables in t
    y_d(:, t) = lam_d * F(1:Nr, t) + sqrt(sig2_d) .* randn(Nd, 1);
    y_w(:, t) = lam_w * f_w(:, t) + sqrt(sig2_w) .* randn(Nw, 1);
    y_m(:, t) = lam_m * f_m(:, t) + sqrt(sig2_m) .* randn(Nm, 1);
    y_q(~ind_q_flow, t) = lam_q(~ind_q_flow, :) * f_q(:, t) + sqrt(sig2_q(~ind_q_flow)) .* randn(Nq_stock, 1);
    y_q(ind_q_flow, t) = lam_q(ind_q_flow, :) * f_q_c(:, t) + sqrt(sig2_q(ind_q_flow)) .* randn(Nq_flow, 1);   
end

% extract factor from companion form representation
f = F(1:Nr, :); 

% extract actually observed monthly and quarterly values
y_d_o = y_d;
tmp = 1:Nt;
ind_w_o = tmp(Xi_w==0) - 1;
ind_w_o = ind_w_o(ind_w_o > 0);
y_w_o = NaN(Nw, Nt);y_w_o(:, tmp(Xi_w==0) - 1) = y_w(:, tmp(Xi_w==0) - 1);
ind_m_o = tmp(Xi_m==0) - 1;
ind_m_o = ind_m_o(ind_m_o > 0);
y_m_o = NaN(Nm, Nt);y_m_o(:, ind_m_o) = y_m(:, ind_m_o);
ind_q_o = tmp(Xi_q==0) - 1;
ind_q_o = ind_q_o(ind_q_o > 0);
y_q_o = NaN(Nq, Nt);y_q_o(:, ind_q_o) = y_q(:, ind_q_o);
clearvars tmp ind_w_o ind_m_o ind_q_o


% plot factors and obs
figure; 
plot(f', 'Color', [0, 0.4470, 0.7410]); 
hold on; 
plot(f_w', 'Color', [0.8500, 0.3250, 0.0980]);
plot(f_m', 'Color', [0.9290, 0.6940, 0.1250]);
plot(f_q')
xticks(ind_plot(5:5:end))
xticklabels(dates_plot(5:5:end))
title('Simulated factors')

figure; 
plot(y_d', 'Color', [0, 0.4470, 0.7410]); 
hold on; 
plot(y_w', '-', 'Color', [0.8500, 0.3250, 0.0980, 0.2], 'LineWidth',0.5); % Color(4) = alpha!
plot(y_w_o', '-o', 'Color', [0.8500, 0.3250, 0.0980]);
plot(y_m', '-', 'Color', [0.9290, 0.6940, 0.1250, 0.2]);
plot(y_m_o', '-o','Color', [0.9290, 0.6940, 0.1250]);
plot(y_q', '-','Color', [0.4940, 0.1840, 0.5560, 0.2]);
plot(y_q_o', '-o','Color', [0.4940, 0.1840, 0.5560]);
xticks(ind_plot(5:5:end))
xticklabels(dates_plot(5:5:end))
title('Simulated observations')


figure;
y_q_flow = y_q(ind_q_flow, :);
p1 = plot(y_q_flow(1, :)', '-', 'Color', [0, 0.4470, 0.7410, 0.7]); 
hold on; 
plot(y_q_flow(2:end, :)', '-', 'Color', [0, 0.4470, 0.7410, 0.7]); 
y_q_stock = y_q(~ind_q_flow, :);
p2 = plot(y_q_stock(1, :)', '-', 'Color', [0.8500, 0.3250, 0.0980, 0.7]);
hold on
plot(y_q_stock(2:end, :)', '-', 'Color', [0.8500, 0.3250, 0.0980, 0.7]);
xticks(ind_plot(5:5:end))
xticklabels(dates_plot(5:5:end))
title('Quarterly series')
legend([p1, p2], {'flows', 'stocks'})
clearvars y_q_flow y_q_stock


%-------------------------------------------------------------------------%
% run Kalman smoother (E-step!)
%-------------------------------------------------------------------------%

% collect parameters in structure
params.lam_d = lam_d;
params.lam_w = lam_w;
params.lam_m = lam_m;
params.lam_q_flow = lam_q(ind_q_flow, :);
params.lam_q_stock = lam_q(~ind_q_flow, :);
params.sig2_d = sig2_d;
params.sig2_w = sig2_w;
params.sig2_m = sig2_m;
params.sig2_q = sig2_q;
params.Xi_w = Xi_w;
params.Xi_m = Xi_m;
params.Xi_q = Xi_q;
params.W_qd_c = W_qd_c;
params.W_qd_p = W_qd_p;
params.Phi = Phi;
params.Omeg = Omeg;

% state space form
[Z, H, T, R, Q] = f_state_space_params(params);

% Kalman smoother
s0 = zeros(size(Z, 2), 1); P0 = 10 * eye(size(Z, 2));
tic
[stT,PtT,LL] = f_KS_DK_logL([y_d_o; y_w_o; y_m_o; y_q_o],T,Z,H,R,Q,s0,P0);
toc

% plot actual and sampled states
nrow_plot = double(Nd > 0) + double(Nw > 0) + double(Nm > 0) + double(Nq > 0);

[id_f, id_f_lags, id_f_d, id_f_w, id_f_m, id_f_q] = f_id_fs(Nr, Np_eff, Nd, Nw, Nm, Nq); % get positions of factors in state vector

figure;
counter = 1;
if Nd > 0
subplot(nrow_plot,1,counter)
plot(stT(id_f_d, :)', 'b')
hold on
plot(f', 'r')
title('f_d (sampled in blue)')
xticks(ind_plot(5:5:end))
xticklabels(dates_plot(5:5:end))
counter = counter + 1;
end

if Nw > 0
subplot(nrow_plot,1,counter)
plot(stT(id_f_w, :)', 'b')
hold on
plot(f_w', 'r')
title('f_w (sampled in blue)')
xticks(ind_plot(5:5:end))
xticklabels(dates_plot(5:5:end))
counter = counter + 1;
end

if Nm > 0
subplot(nrow_plot,1,counter)
plot(stT(id_f_m, :)', 'b')
hold on
plot(f_m', 'r')
title('f_m (sampled in blue)')
xticks(ind_plot(5:5:end))
xticklabels(dates_plot(5:5:end))
counter = counter + 1;
end

if Nq > 0
subplot(nrow_plot,1,counter)
plot(stT(id_f_q, :)', 'b')
hold on
plot(f_q', 'r')
title('f_q (sampled in blue)')
xticks(ind_plot(5:5:end))
xticklabels(dates_plot(5:5:end))
end

%-------------------------------------------------------------------------%
% estimate parameters (M-step)
%-------------------------------------------------------------------------%

% get positions of factors in state vector
[id_f, id_f_lags, id_f_d, id_f_w, id_f_m, id_f_q] = f_id_fs(Nr, Np_eff, Nd, Nw, Nm, Nq);

% lam_d and sig_d
if Nd > 0
    [lam_d_hat, sig_d_hat] = f_sample_lam_sig(y_d_o, stT(id_f_d, :), PtT(id_f_d, id_f_d, :), sig2_d);
end

% lam_w and sig2_w
if Nw > 0
    [lam_w_hat, sig2_w_hat] = f_sample_lam_sig(y_w_o, stT(id_f_w, :), PtT(id_f_w, id_f_w, :), sig2_w);
end

% lam_m and sig2_m
if Nm > 0
    [lam_m_hat, sig2_m_hat] = f_sample_lam_sig(y_m_o, stT(id_f_m, :), PtT(id_f_m, id_f_m, :), sig2_m);
end

% lam_q and sig_q
if Nq > 0
    [lam_q_hat, sig2_q_hat] = f_sample_lam_sig(y_q_o, stT(id_f_q, :), PtT(id_f_q, id_f_q, :), sig2_q);
end

% Phi and Omeg
Phi_hat = (stT(id_f, :)*stT(id_f_lags, :)' + sum(PtT(id_f,id_f_lags,:),3))/(stT(id_f_lags, :)*stT(id_f_lags, :)' + sum(PtT(id_f_lags,id_f_lags,:),3)) ; 
Omeg_hat = 1/Nt * ((stT(id_f, :)*stT(id_f, :)' + sum(PtT(id_f,id_f,:),3))  - Phi_hat * (stT(id_f, :)*stT(id_f_lags, :)' + sum(PtT(id_f,id_f_lags,:),3))') ;
 
% plot estimated versus true parameters
nrow_plot = 1 + double(Nd > 0) + double(Nw > 0) + double(Nm > 0) + double(Nq > 0); % additional row for Phi and Omeg!
figure;
counter = 1;
subplot(nrow_plot,2,counter)
scatter(Phi(:), Phi_hat(:))
ylim([-1.0 1.0])
xlim([-1.0 1.0])
refline(1, 0)
title('Phi')
ylabel('estimate')
xlabel('actual')
counter = counter + 1;

subplot(nrow_plot,2,counter)
scatter(Omeg(:), Omeg_hat(:))
ylim([-1.5 1.5])
xlim([-1.5 1.5])
refline(1, 0)
title('Omeg')
ylabel('estimate')
xlabel('actual')

counter = counter + 1;

if Nd > 0
    subplot(nrow_plot,2,counter)
    scatter(lam_d(:, 1), lam_d_hat(:, 1), 'b')
    if Nr == 2
        hold on;
        scatter(lam_d(:, 2), lam_d_hat(:, 2), 'r')
    end
    ylim([-1 1])
    xlim([-1 1])
    refline(1, 0)
    title('lam_d')
    ylabel('estimate')
    xlabel('actual')
    counter = counter + 1;

    subplot(nrow_plot,2,counter)
    scatter(sig2_d, sig_d_hat)
    ylim([0 1])
    xlim([0 1])
    refline(1, 0)
    title('sig_d')
    ylabel('estimate')
    xlabel('actual')
    counter = counter + 1;
end

if Nw > 0
    subplot(nrow_plot,2,counter)
    scatter(lam_w(:, 1), lam_w_hat(:, 1), 'b')
    hold on;
    scatter(lam_w(:, 2), lam_w_hat(:, 2), 'r')
    ylim([-1 1])
    xlim([-1 1])
    refline(1, 0)
    title('lam_w')
    ylabel('estimate')
    xlabel('actual')
    counter = counter + 1;

    subplot(nrow_plot,2,counter)
    scatter(sig2_w, sig2_w_hat)
    ylim([0 1])
    xlim([0 1])
    refline(1, 0)
    title('sig_w')
    ylabel('estimate')
    xlabel('actual')
    counter = counter + 1;
end

if Nm > 0
    subplot(nrow_plot,2,counter)
    scatter(lam_m(:, 1), lam_m_hat(:, 1), 'b')
    hold on;
    scatter(lam_m(:, 2), lam_m_hat(:, 2), 'r')
    ylim([-1 1])
    xlim([-1 1])
    refline(1, 0)
    title('lam_m')
    ylabel('estimate')
    xlabel('actual')
    counter = counter + 1;

    subplot(nrow_plot,2,counter)
    scatter(sig2_m, sig2_m_hat)
    ylim([0 1])
    xlim([0 1])
    refline(1, 0)
    title('sig_m')
    ylabel('estimate')
    xlabel('actual')
    counter = counter + 1;
end

if Nq > 0
    subplot(nrow_plot,2,counter)
    scatter(lam_q(ind_q_flow, 1), lam_q_hat(ind_q_flow, 1), 'bx')
    hold on
    scatter(lam_q(~ind_q_flow, 1), lam_q_hat(~ind_q_flow, 1), 'bo')
    if Nr == 2
        scatter(lam_q(ind_q_flow, 2), lam_q_hat(ind_q_flow, 2), 'rx')
        scatter(lam_q(~ind_q_flow, 2), lam_q_hat(~ind_q_flow, 2), 'ro')
    end
    ylim([-1 1])
    xlim([-1 1])
    refline(1, 0)
    title('lam_q')
    ylabel('estimate')
    xlabel('actual')
    counter = counter + 1;

    subplot(nrow_plot,2,counter)
    scatter(sig2_q(ind_q_flow), sig2_q_hat(ind_q_flow), 'bx')
    hold on
    scatter(sig2_q(~ind_q_flow), sig2_q_hat(~ind_q_flow), 'bo')
    ylim([0 1])
    xlim([0 1])
    refline(1, 0)
    title('sig_q')
    ylabel('estimate')
    xlabel('actual')
    counter = counter + 1;
end
