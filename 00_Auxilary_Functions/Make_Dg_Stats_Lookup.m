% this script calculate Dg for M0 for a range of sigma_r, so that it can be
% used as a look up table to impute from.

function Make_Dg_Stats_Lookup(OutPath)
N = 500;
sigma_r = [linspace(0, 1.501, N)];

disp('Pre-calculating Expected_Expression_MAD_BLN_v3')
for i = N:-1:1
    [Dg(i), Vg(i)]= Expected_Expression_MAD_BLN_v3(0, sigma_r(i), 0);
    
    if i/(N/20)==round(i/(N/20))
        fprintf('%d ', i);
    end
end

save(OutPath, 'sigma_r', 'Dg', 'Vg')
disp('done')
end