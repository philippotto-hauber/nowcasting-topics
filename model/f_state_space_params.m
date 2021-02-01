function [Z, H, T, R, Q] = f_state_space_params(params)

% back out dimensions of state space system
Nd = size(params.lam_d, 1);
Nw = size(params.lam_w, 1);
Nm = size(params.lam_m, 1);
Nq = size(params.lam_q, 1);
Nr = size(params.Phi, 1);
Np_eff = size(params.Phi, 2)/Nr + 1;
Nt = size(params.Xi_m, 1); 

% determine size of state vector Ns
Ns = Nr * Np_eff;  
if Nw > 0 && Nd > 0; Ns = Ns+Nr;end 
if Nm > 0 && (Nw > 0 || Nd > 0) ; Ns = Ns+Nr;end 
if Nq > 0 && (Nm > 0 || Nw > 0 || Nd > 0); Ns = Ns+Nr;end 

% T, => time-varying
T0 = eye(Ns);
T0(Nr*Np_eff+1:Ns,1:Nr) = repmat(-eye(Nr), (Ns - Nr*Np_eff)/Nr, 1);
iT0 = T0 \ eye(Ns); 

T = NaN(Ns, Ns, Nt);
for t = 1:Nt
    Ttmp = [params.Phi zeros(Nr, Ns-size(params.Phi, 2));
            eye(Nr * (Np_eff-1)) zeros(Nr * (Np_eff-1), Ns - Nr * (Np_eff-1))];
    if Nw > 0 && Nd > 0 % if there are weekly series and its not the highest frequency. Otherwise its dynamics are governed by phi_f!
        Ttmp = [Ttmp; zeros(Nr, size(Ttmp, 1)) params.Xi_w(t) * eye(Nr) zeros(Nr, Ns - (size(Ttmp, 1) + Nr))];
    end
    if Nm > 0 && (Nw > 0 || Nd > 0) 
        Ttmp = [Ttmp; zeros(Nr, size(Ttmp, 1)) params.Xi_m(t) * eye(Nr) zeros(Nr, Ns - (size(Ttmp, 1) + Nr))];
    end
    if Nq > 0 && (Nm > 0 || Nw > 0 || Nd > 0)
        Ttmp = [Ttmp; zeros(Nr, size(Ttmp, 1)) params.Xi_q(t) * eye(Nr) zeros(Nr, Ns - (size(Ttmp, 1) + Nr))];
    end

    T(:, :, t) = iT0 * Ttmp;
end

% Q & R 
Q = zeros(Ns);
Q(1:Nr, 1:Nr) = params.Omeg; 
R = iT0; 

% Z 
Z_d = []; if Nd > 0; Z_d = [params.lam_d; zeros(Nw+Nm+Nq, Nr)]; end
Z_w = []; if Nw > 0; Z_w = [zeros(Nd, Nr); params.lam_w; zeros(Nm+Nq, Nr)]; end
Z_m = []; if Nm > 0; Z_m = [zeros(Nd+Nw, Nr); params.lam_m; zeros(Nq, Nr)]; end
Z_q = []; if Nq > 0; Z_q = [zeros(Nd+Nw+Nm, Nr); params.lam_q]; end
Z = [Z_d zeros(Nd+Nw+Nm+Nq, Nr * (Np_eff-1)) Z_w Z_m Z_q]; 

% H
H = diag([params.sig2_d; params.sig2_w; params.sig2_m; params.sig2_q]);

