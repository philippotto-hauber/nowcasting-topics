clear; close all; clc;
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%- This script....


%-------------------------------------------------------------------------%
% user settings
%-------------------------------------------------------------------------%
filename = 'vint_2010_1_30.csv';
dirname = '..\data\';
Nr = 1;
Np = 6;

%-------------------------------------------------------------------------%
% load data
%-------------------------------------------------------------------------%
tmp = importdata([dirname, filename]);

% offset as there are less numeric columns!
offset_numcols = size(tmp.textdata, 2) - size(tmp.data, 2);

% back out data for estimation and forecasting
aux.ind_sample = logical(tmp.data(:, find(strcmp('ind_sample', tmp.textdata(1,:))) - offset_numcols));

% daily data
ind_y_d = find(contains(tmp.textdata(1,:), 'y_d_')) - offset_numcols;
%ind_y_d = setdiff(ind_y_d, ind_y_d([6, 9, 23])); % manually remove T05, T07, T21
ind_y_d = [1 11 22 37 45] + 4; % topics T0, T10, T21, T36, T44 => highest correlated with GDP
y_d = tmp.data(aux.ind_sample, ind_y_d)';
y_d_fore = tmp.data(~aux.ind_sample, ind_y_d)';

y_d = f_interpol(y_d); % linearly interpolate so that there are no missings!

% quarterly data
ind_y_q = find(contains(tmp.textdata(1,:), 'y_q_')) - offset_numcols;
y_q = tmp.data(aux.ind_sample, ind_y_q)';
y_q_fore = tmp.data(~aux.ind_sample, ind_y_q)';
aux.ind_q_flow = 1; 

% weights and Xi
aux.Xi_qd = tmp.data(:, find(strcmp('Xi_qd', tmp.textdata(1,:))) - offset_numcols);
aux.W_qd_p = tmp.data(:, find(strcmp('W_qd_p', tmp.textdata(1,:))) - offset_numcols);
aux.W_qd_c = tmp.data(:, find(strcmp('W_qd_c', tmp.textdata(1,:))) - offset_numcols);

% inds for back-, now- and forecasts
ind_backcast = logical(tmp.data(:, find(strcmp('ind_backcast', tmp.textdata(1,:))) - offset_numcols));
ind_nowcast = logical(tmp.data(:, find(strcmp('ind_nowcast', tmp.textdata(1,:))) - offset_numcols));
ind_forecast = logical(tmp.data(:, find(strcmp('ind_forecast1Q', tmp.textdata(1,:))) - offset_numcols));

%-------------------------------------------------------------------------%
% prepare data for estimation
%-------------------------------------------------------------------------%

% standardize
y_d_stand = (y_d - nanmean(y_d, 2)) ./ nanstd(y_d, [], 2);
y_d_fore_stand = (y_d_fore - nanmean(y_d, 2)) ./ nanstd(y_d, [], 2);
mean_gdp = nanmean(y_q);
std_gdp = nanstd(y_q);
y_q_stand = (y_q - mean_gdp) / std_gdp; 
y_q_fore_stand = (y_q_fore - mean_gdp) / std_gdp; 

% starting values
params = f_start_vals(y_d_stand, [], [], y_q_stand, aux, Nr, Np);

%-------------------------------------------------------------------------%
% EM algorithm
%-------------------------------------------------------------------------%

params_init = params; 
params = f_EMalg(y_d_stand, [], [], y_q_stand, aux, params); 

%-------------------------------------------------------------------------%
% run KF/KS to get back-, now- and forecasts
%-------------------------------------------------------------------------%
  
dat = [[y_d_stand y_d_fore_stand]; [y_q_stand y_q_fore_stand]]; 
[Z, H, T, R, Q] = f_state_space_params(params, aux, size(dat, 2));
s0 = zeros(size(T,1),1); 
P0 = 1 * eye(size(T,1)); 
[stT, ~, ~] = f_KS_DK_logL(dat,T,Z,H,R,Q,s0,P0);

figure;
plot(stT(1:Nr,:)')
title('daily factors')

figure;
plot(stT(end-2*Nr+1:end-Nr,:)')
title('f_q: cumulated daily factors')

%-------------------------------------------------------------------------%
% plot daily q-o-q growth along with actuals 
%-------------------------------------------------------------------------%

gdp_hat_stand = params.lam_q_flow * stT(end-2*Nr+1:end-Nr,:);
gdp_hat = gdp_hat_stand * std_gdp + mean_gdp;
Nt = sum(aux.ind_sample); 
Nh = length(aux.ind_sample) - Nt; % sum(aux.ind_sample == 0)
gdp_fore = NaN(1, Nt+Nh); 
gdp_fore(1, ind_backcast) = gdp_hat(1, ind_backcast); 
gdp_fore(1, ind_nowcast) = gdp_hat(1, ind_nowcast); 
gdp_fore(1, ind_forecast) = gdp_hat(1, ind_forecast); 
gdp_fore(1, end) = gdp_hat(1, end);

% actuals
acts = NaN(1, length(aux.ind_sample));
acts(1, ind_backcast) = -0.04;
acts(1, ind_nowcast) = 0.64;
acts(1, ind_forecast) = 9.1;
acts(1, end) = 2.8;

% back out dates for plot
ind_plot = find(tmp.data(:, 2) == 1 & tmp.data(:, 3) == 1 & tmp.data(:, 4) == 1);
dates_plot = tmp.data(ind_plot, 1); 


figure; 
plot([gdp_hat(:, 1: sum(aux.ind_sample)) NaN(1, Nh)]', 'b-')
hold on
p1 = plot([y_q, NaN(1, Nh)]', 'bo');
plot([NaN(1, Nt), gdp_hat(:, sum(aux.ind_sample)+1:end)]', 'r-')
p2 = plot(gdp_fore', 'ro');
p3 = plot(acts, 'kx');
xticks(ind_plot(1:5:end))
xticklabels(dates_plot(1:5:end))
legend([p1, p2, p3], {'in-sample', 'out-of-sample', 'actual'}, 'Location','SouthWest')
title('quarterly GDP growth (ann.), forecasts and actuals')




