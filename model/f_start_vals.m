function params = f_start_vals(y_d, y_w, y_m, y_q, aux, Nr)

%-----------------------------------%
%- lam_d, sig2_d
%-----------------------------------%

% initial PCA estimate of daily factors => replace with EM alg a la Stock Watson
f = f_PCA(y_d, Nr);

% scale factors to unit variance
f  = f ./ std(f, [], 2);

% lam_d, sig2_d
[params.lam_d, params.sig2_d] = f_ols(y_d, f, zeros(size(f, 2), 1));

%-----------------------------------%
%- lam_w, sig2_w
%-----------------------------------%

% cumulate factors for weekly stocks
[f_w, ~] = f_cum_f(f, aux.Xi_wd, zeros(size(aux.Xi_wd)), zeros(size(aux.Xi_wd)));

% lam_w and sig2_w
[params.lam_w, params.sig2_w] = f_ols(y_w, f_w, aux.Xi_wd);

%-----------------------------------%
%- lam_m, sig2_m
%-----------------------------------%

% cumulate factors for monthly stocks and flows
[f_m, f_m_c] = f_cum_f(f, aux.Xi_md, aux.W_md_c, aux.W_md_p);

% lam_m_stock and sig2_m_stock
[params.lam_m_stock, params.sig2_m_stock] = f_ols(y_m(~aux.ind_m_flow, :), f_m, aux.Xi_md);

% lam_m_stock and sig2_m_stock
[params.lam_m_flow, params.sig2_m_flow] = f_ols(y_m(aux.ind_m_flow, :), f_m_c, aux.Xi_md);

%-----------------------------------%
%- lam_q, sig2_q
%-----------------------------------%

% cumulate factors for quarterly stocks and flows
[f_q, f_q_c] = f_cum_f(f, aux.Xi_qd, aux.W_qd_c, aux.W_qd_p);

% lam_q_stock and sig2_q_stock
[params.lam_q_stock, params.sig2_q_stock] = f_ols(y_q(~aux.ind_q_flow, :), f_q, aux.Xi_qd);

% lam_q_flow and sig2_q_flow
[params.lam_q_flow, params.sig2_q_flow] = f_ols(y_q(aux.ind_q_flow, :), f_q_c, aux.Xi_qd);
end

function F_hat = f_PCA(Y,R)
% covariance matrix of observables
SIGMA = Y*Y'/size(Y,2);

% eigenvalue and -vector decomposition
[V,D] = eig(SIGMA);

% extract eigenvalues and change order
eigenvalues_D = diag(D);
eigenvalues_D = flipud(eigenvalues_D);
D = diag(eigenvalues_D);
% change order of eigenvectors
V = f_rev_cols(V);
F_hat = V(:, 1:R)'*Y ;
F_hat = F_hat(1:R,:);
end

function [ A_reverse ] = f_rev_cols(A)
%Reverses the columns of a matrix. 
aux = zeros(size(A));
[R,C] = size(A);
for c=0:(C-1)
    aux(:,c+1) = A(:,C-c);
end
A_reverse = aux;
end

function [lam, sig2] = f_ols(y, f, Xi)
Nt = size(f, 2); 
tmp = 1:Nt;
if all(Xi == 0)
    ind_o = repelem(true, Nt);
else
    ind_o = tmp(Xi == 0) - 1;
    ind_o = ind_o(ind_o > 0);
end
x_estim = f(:, ind_o)';
y_estim = y(:, ind_o)';

b = (x_estim' * x_estim)\ x_estim' * y_estim;
lam = b'; 
sig2 = var(y_estim - x_estim * b)';
end

function [f_, f_c] = f_cum_f(f, Xi, Wc, Wp)
    [Nr, Nt] = size(f);
    f_ = NaN(Nr, Nt);
    f_c = NaN(Nr, Nt);
    f_p = zeros(Nr, 1);
    for t = 1:Nt
        if t == 1 || Xi(t) == 0
            f_(:, t) = zeros(Nr, 1);
            f_c(:, t) = f_p;
            f_p = zeros(Nr, 1);
        else
            f_(:, t) = f_(:, t-1) + f(:, t);
            f_c(:, t) = f_c(:, t-1) + Wc(t) * f(:, t);
            f_p = f_p + Wp(t) * f(:, t);
        end
    end
end

