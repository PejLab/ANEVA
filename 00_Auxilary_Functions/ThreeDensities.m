function W= ThreeDensities(pH_given_r, pH_given_a)
W(1)=    pH_given_r .* (1-pH_given_a); % High expressed eQTL allele on the same hap as the Ref ASE and the Low  on alt.
W(3)= (1-pH_given_r).*    pH_given_a ; % Low  expressed eQTL allele on the same hap as the Ref ASE and tje High on alt.

W(2)=  1-(W(1)+W(3)); % eQTL is Hom
end