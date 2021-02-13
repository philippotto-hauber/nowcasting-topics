function params = f_EMalg(y_d, y_w, y_m, y_q, aux, params)
% ------------------------------------------------------------------- %
% - EM algorithm a la Banbura and Modugno (2014) -------------------- %
% ------------------------------------------------------------------- %

% back out dimensions of state space system
Nd = size(params.lam_d, 1);
Nw = size(params.lam_w, 1);
Nm_flow = size(params.lam_m_flow, 1);
Nm_stock = size(params.lam_m_stock, 1);
Nm = Nm_flow + Nm_stock;
Nq_flow = size(params.lam_q_flow, 1);
Nq_stock = size(params.lam_q_stock, 1);
Nq = Nq_flow + Nq_stock;
Nr = size(params.Phi, 1);
Np_eff = size(params.Phi, 2)/Nr + 1;
Nt = size(aux.Xi_qd, 1); 

% iterations
maxiter = 50; 
LL_prev = -999999 ; 
%tic
for iter = 1 : maxiter 
    % ------------------------------------------------------------------- %
    % - E-Step ---------------------------------------------------------- %
    % ------------------------------------------------------------------- %

    % - Kalman Smoother
    % ---------------------------    
    dat = [y_d; y_w; y_m; y_q];
    [Z, H, T, R, Q] = f_state_space_params(params, aux, size(dat, 2));
    s0 = zeros(size(T,1),1); 
    P0 = 1 * eye(size(T,1)); 
    [stT,PtT,LL] = f_KS_DK_logL(dat,T,Z,H,R,Q,s0,P0) ;
    
    % - check convergence
    % ---------------------------  
    cLL = (LL - LL_prev)/(abs(LL)/2 + abs(LL_prev)/2) ; 
    if iter>1 && cLL < 1e-03; break; end    
    LL_prev = LL;
    
    % ------------------------------------------------------------------- %
    % - M-Step ---------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    
    % get positions of factors in state vector
    [id_f, id_f_lags, id_f_d, id_f_w, id_f_m_flow, id_f_m_stock, id_f_q_flow, id_f_q_stock] = f_id_fs(Nr, Np_eff, Nd, Nw, Nm_flow, Nm_stock, Nq_flow, Nq_stock);

    % lam_d and sig_d
    if Nd > 0
        [params.lam_d, params.sig_d] = f_sample_lam_sig(y_d, stT(id_f_d, :), PtT(id_f_d, id_f_d, :), params.sig2_d);
    end

    % lam_w and sig2_w
    if Nw > 0
        [params.lam_w, params.sig2_w] = f_sample_lam_sig(y_w, stT(id_f_w, :), PtT(id_f_w, id_f_w, :), params.sig2_w);
    end

    % lam_m and sig2_m
    if Nm > 0
        if Nm_flow > 0
            [params.lam_m_flow, params.sig2_m_flow] = f_sample_lam_sig(y_m(ind_m_flow,:), stT(id_f_m_flow, :), PtT(id_f_m_flow, id_f_m_flow, :), params.sig2_m(ind_m_flow));
        end
        if Nm_stock > 0
            [params.lam_m_stock, params.sig2_m_stock] = f_sample_lam_sig(y_m(~ind_m_flow,:), stT(id_f_m_stock, :), PtT(id_f_m_stock, id_f_m_stock, :), params.sig2_m(~ind_m_flow));
        end
        params.lam_m = [params.lam_m_flow; params.lam_m_stock];
        params.sig2_m = [params.sig2_m_flow; params.sig2_m_stock];
    end

    % lam_q and sig_q
    if Nq > 0
        if Nq_flow > 0         
            [params.lam_q_flow, params.sig2_q_flow] = f_sample_lam_sig(y_q(aux.ind_q_flow,:), stT(id_f_q_flow, :), PtT(id_f_q_flow, id_f_q_flow, :), params.sig2_q(aux.ind_q_flow));
        end

        if Nq_stock > 0         
            [params.lam_q_stock, params.sig2_q_stock] = f_sample_lam_sig(y_q(~ind_q_flow,:), stT(id_f_q_stock, :), PtT(id_f_q_stock, id_f_q_stock, :), params.sig2_q(~ind_q_flow));
        end

        params.lam_q = [params.lam_q_flow; params.lam_q_stock];
        params.sig2_q = [params.sig2_q_flow; params.sig2_q_stock];
    end

    % Phi and Omeg
    params.Phi = (stT(id_f, :)*stT(id_f_lags, :)' + sum(PtT(id_f,id_f_lags,:),3))/(stT(id_f_lags, :)*stT(id_f_lags, :)' + sum(PtT(id_f_lags,id_f_lags,:),3)) ; 
    params.Omeg = 1/Nt * ((stT(id_f, :)*stT(id_f, :)' + sum(PtT(id_f,id_f,:),3))  - params.Phi * (stT(id_f, :)*stT(id_f_lags, :)' + sum(PtT(id_f,id_f_lags,:),3))') ;
end

if iter == maxiter
    disp('% ------------------------------------------------ %')
    disp(['% EM algorithm failed to converge after ' num2str(iter) ' iterations'])
    disp('% ------------------------------------------------ %')
else    
    disp('% ------------------------------------------------ %')
    disp(['% EM algorithm converged after ' num2str(iter) ' iterations'])
    disp('% ------------------------------------------------ %')
end
end

