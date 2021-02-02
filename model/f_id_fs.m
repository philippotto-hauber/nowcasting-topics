function [id_f, id_f_lags, id_f_d, id_f_w, id_f_m, id_f_q] = f_id_fs(Nr, Np_eff, Nd, Nw, Nm, Nq)

id_f = 1:Nr; % note that id_f will be identical to id_f_d, id_f_w, ... always refering to the factor of the highest frequency!
id_f_lags = Nr+1:Nr * Np_eff;

ind_s = 0;
if Nd > 0
    id_f_d = ind_s+1:ind_s+Nr;
    ind_s = Nr * Np_eff;
else
    id_f_d = [];
end

if Nw > 0 && Nd > 0
   id_f_w = ind_s + 1 : ind_s+Nr;
   ind_s = ind_s + Nr; 
elseif Nw > 0 && Nd == 0 % w is first available frequency
   id_f_w = ind_s + 1 : ind_s + Nr;
   ind_s = Nr * Np_eff + 1;
else
    id_f_w = [];
end

if Nm > 0 && (Nw > 0 || Nd > 0)
   id_f_m = ind_s + 1: ind_s+Nr;
   ind_s = ind_s + Nr;
elseif Nm > 0 && Nw == 0 && Nd == 0 % m is first available series
   id_f_m = ind_s+1 : ind_s+Nr;
   ind_s = Nr * Np_eff + 1;
else
   id_f_m = [];
end

if Nq > 0
    id_f_q = ind_s+1 : ind_s+Nr;
else
    id_f_q = []; 
end

 



