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
Nr = 2;
Np = 3;

%-------------------------------------------------------------------------%
% load data
%-------------------------------------------------------------------------%
tmp = importdata([dirname, filename]);

% offset as there are less numeric columns!
offset_numcols = size(tmp.textdata, 2) - size(tmp.data, 2);

% back out data for estimation and forecasting
aux.ind_sample = logical(tmp.data(:, find(strcmp('"ind_sample"', tmp.textdata(1,:))) - offset_numcols));

% daily data
ind_y_d = find(contains(tmp.textdata(1,:), '"y_d_')) - offset_numcols;
y_d = tmp.data(aux.ind_sample, ind_y_d)';
y_d_fore = tmp.data(~aux.ind_sample, ind_y_d)';

% quarterly data
ind_y_q = find(contains(tmp.textdata(1,:), '"y_q_')) - offset_numcols;
y_q = tmp.data(aux.ind_sample, ind_y_q)';
y_q_fore = tmp.data(~aux.ind_sample, ind_y_q)';
aux.ind_q_flow = 1; 

% weights and Xi
aux.Xi_qd = tmp.data(:, find(strcmp('"Xi_qd"', tmp.textdata(1,:))) - offset_numcols);
aux.W_qd_p = tmp.data(:, find(strcmp('"W_qd_p"', tmp.textdata(1,:))) - offset_numcols);
aux.W_qd_c = tmp.data(:, find(strcmp('"W_qd_c"', tmp.textdata(1,:))) - offset_numcols);

% inds for back-, now- and forecasts
ind_backcast = logical(tmp.data(:, find(strcmp('"ind_backcast"', tmp.textdata(1,:))) - offset_numcols));
ind_nowcast = logical(tmp.data(:, find(strcmp('"ind_nowcast"', tmp.textdata(1,:))) - offset_numcols));

%-------------------------------------------------------------------------%
% prepare data for estimation
%-------------------------------------------------------------------------%

% back out dims ?!

% standardize
y_d_stand = (y_d - nanmean(y_d, 2)) ./ nanstd(y_d, [], 2);
mean_gdp = nanmean(y_q);
std_gdp = nanstd(y_q);
y_q_stand = (y_q - mean_gdp) / std_gdp; 

% starting values
params = f_start_vals(y_d_stand, [], [], y_q_stand, aux, Nr, Np);

%-------------------------------------------------------------------------%
% EM algorithm
%-------------------------------------------------------------------------%

params_init = params; 
params = f_EMalg(y_d, [], [], y_q, aux, params); 


%-------------------------------------------------------------------------%
% run KF/KS to get back-, now- and forecasts
%-------------------------------------------------------------------------%
