%  This produces results that are very similar to the v2 version. V2 used a
%  2D integral, here I use a 1D integral/ 

% Dg: Mean Abs Dev in natural log scale due to genetic effects
% Vg: Variance     in natural log scale due to genetic effects

function [Dg, Vg] = Expected_Expression_MAD_BLN_v3(s_H, Std_sr, P_H)
Dg = nan(size(s_H));
Vg = nan(size(s_H));

for i = 1:numel(s_H)
    [Dg(i), Vg(i)] = Estimate_Dg_forone(s_H(i), Std_sr(i), P_H(i));
end

end

function [Dg, Vg] = Estimate_Dg_forone(s, Std_sr, P_H)
if isnan(Std_sr+P_H+s)
    Dg = nan;
    Vg = nan;
    return
end
if Std_sr==0
    Std_sr=eps; % to avoid numerical errors
end

e0 = 0;
E_e = Expected_Expression(s, Std_sr, e0, P_H);

Frqs = pgenotype(P_H);

Dg = integral(@(sr)dMAD_ar(sr, Std_sr, s, e0, Frqs, E_e), -Std_sr*4,+Std_sr*4);%,'RelTol',1e-3,'AbsTol',1e-4);
Vg = integral(@(sr)dVAR_ar(sr, Std_sr, s, e0, Frqs, E_e), -Std_sr*4,+Std_sr*4);%,'RelTol',1e-3,'AbsTol',1e-4);
end

function E = Expected_Expression(s, Std_sr, e0, P_H)
Frqs = pgenotype(P_H);

E = integral(@(sr)dE_expr(sr, Std_sr, s, e0, Frqs), -Std_sr*4,+Std_sr*4);%,'RelTol',1e-3,'AbsTol',1e-4);
end

function dE = dE_expr(sr, Std_sr, s, e0, Frqs)
[e_LL, e_LH, e_HL, e_HH] = dExp_sr(sr, s, e0);

% Frqs = [p_LL, p_LH, p_HL, p_HH];
tE = ...
    e_LL*Frqs(1) + ...
    e_LH*Frqs(2) + ...
    e_HL*Frqs(3) + ...
    e_HH*Frqs(4) ;

pE = p_sr(sr,  0, Std_sr);
dE = tE.* pE;
end

function [e_LL, e_LH, e_HL, e_HH] = dExp_sr(sr, s, e0)
e_LL = log(exp(e0+sr)  + exp(e0));
e_LH = log(exp(e0+sr)  + exp(e0+s));     % noise on low  expressed hap
e_HL = log(exp(e0)     + exp(e0+s+sr));  % noise on high expressed hap
e_HH = log(exp(e0+s+sr)+ exp(e0+s));
end

function dDg = dMAD_ar(sr, Std_sr,   s, e, Frqs, E_e)
[e_LL, e_LH, e_HL, e_HH] = dExp_sr(sr, s, e);

% Frqs = [p_LL, p_LH, p_HL, p_HH];
tDg = ...
    abs(e_LL-E_e)*Frqs(1) + ...
    abs(e_LH-E_e)*Frqs(2) + ...
    abs(e_HL-E_e)*Frqs(3) + ...
    abs(e_HH-E_e)*Frqs(4) ;

pE = p_sr(sr,  0, Std_sr);
dDg = tDg.* pE;
end

function dVg = dVAR_ar(sr, Std_sr,   s, e, Frqs, E_e)
[e_LL, e_LH, e_HL, e_HH] = dExp_sr(sr, s, e);

% Frqs = [p_LL, p_LH, p_HL, p_HH];
tVg = ...
    ((e_LL-E_e).^2)*Frqs(1) + ...
    ((e_LH-E_e).^2)*Frqs(2) + ...
    ((e_HL-E_e).^2)*Frqs(3) + ...
    ((e_HH-E_e).^2)*Frqs(4) ;

pE = p_sr(sr,  0, Std_sr);
dVg = tVg.* pE;
end

function p = p_sr(sr,mu,Std)
p = normpdf(sr,mu,Std);
end

function Frqs = pgenotype(P_H)
p_LL = (1-P_H).^2;    % probability of LL genotype (hom low expressed eQTL)
p_LH = P_H.*(1-P_H);  % probability of LH genotype
p_HL = P_H.*(1-P_H);  % probability of HL genotype
p_HH = P_H.^2;        % probability of HH genotype (hom high expressed eQTL)

Frqs = [p_LL, p_LH, p_HL, p_HH];
Frqs = Frqs./sum(Frqs);
end