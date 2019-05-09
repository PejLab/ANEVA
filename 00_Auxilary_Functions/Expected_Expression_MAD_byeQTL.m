function  [Dg, Vg] = Expected_Expression_MAD_byeQTL(log2_aFC_A2R, fR)
% fR: Allele freq of the ref eqtl allele
% log2_aFC_A2R: *log2* aFC of the eQTL Alt to eQTL Ref allele

% Dg: Mean Abs Dev in *natural* log scale due to genetic effects
% Vg: Variance     in *natural* log scale due to genetic effects

E_r = zeros(size(log2_aFC_A2R)); % in log2 scale
E_a = E_r + log2_aFC_A2R; % in log2 scale

E_rr = E_r + 1; % in log2 scale
E_aa = E_a + 1; % in log2 scale
E_ar = log2(2.^E_r + 2.^E_a); % in log2 scale

Frr = fR.^2; % pop frequency for RefRef eQTL genotype
Faa = (1-fR).^2;
Far = 1 - (Frr+Faa);

EE = E_rr.*Frr ...
    + E_aa.*Faa ...
    + E_ar.*Far; % expected log2 mean totaL expression

Dg =  abs(E_rr-EE).*Frr ...
    + abs(E_aa-EE).*Faa ...
    + abs(E_ar-EE).*Far; % expectd log2 mean abs. dev. from the mean in total expression

Vg =  (E_rr-EE).^2.*Frr ...
    + (E_aa-EE).^2.*Faa ...
    + (E_ar-EE).^2.*Far; % expectd log2 var from the mean in total expression

Dg = Dg .*  log(2); % in natural log scale
Vg = Vg .* (log(2)).^2; % in natural log scale
end