

function C = Visualize_BLNfits(X, Xc, MLE, Label, PointLabels)
hold on
K = length(MLE.P); % number of classes
isM = ~isnan(MLE.NLL);

CBOX = jet(101); CBOX = CBOX(end:-1:1,:);
%% classify
if isM
   [MAP, C] = get_MAP_Class_BLN(X, Xc, MLE);
   CI = ceil(MLE.P*100)+1;
else
     % there's no model fitted
    MAP = ones(size(X));
    C   = ones(size(X));
    
    MAP(isnan(X)) = nan;
    C(isnan(X)) = nan;
    
    CI = ceil([.50 nan nan]*100)+1;
    CBOX(:)=.5;
end

if isM
    for i = K:-1:1
        Add_Cloud_sqrt(X+Xc, MLE.P(i), MLE.Std(i), CBOX(CI(i),:));
    end
end

for i = K:-1:1
    F = find(C==i);
    if ~isempty(F)
        for j = 1:length(F)
            plot(sqrt(X(F(j))),sqrt(Xc(F(j))), 'o', 'markersize', 6.5, 'linewidth', 2.75*(1-MAP(F(j)))+.5, 'color', [.7 .7 .7], 'markerfacecolor', CBOX(CI(i),:));%, 'displayname', ['AR=' num2str(MLE.P(i)) ',\lambda =' num2str(MLE.W(i))]);
            if nargin<5 % plot PointLabels
                %text(sqrt(X(F(j))),sqrt(Xc(F(j))), num2str(F(j)), 'color', 0*[.7 .7 .7], 'horizontalalignment', 'center', 'verticalalignment', 'bottom');%, 'displayname', ['AR=' num2str(MLE.P(i)) ',\lambda =' num2str(MLE.W(i))]);
            else
                text(sqrt(X(F(j))),sqrt(Xc(F(j))), PointLabels{F(j)}, 'color', 0*[.7 .7 .7], 'fontsize', 5, 'horizontalalignment', 'center', 'verticalalignment', 'bottom');%, 'displayname', ['AR=' num2str(MLE.P(i)) ',\lambda =' num2str(MLE.W(i))]);
            end
        end
    end
end

colormap(gcf, CBOX(2:end,:));
% colorbar
box on
% legend show
xlabel('Ref. Allele count')
ylabel('Alt. Allele count')

xm = min([xlim ylim]);
xM = max([xlim ylim]);

xlim([xm xM]);
ylim([xm xM]);

tmpTicks = get(gca, 'XTick');

set(gca, 'XTick', tmpTicks);
set(gca, 'YTick', tmpTicks);

set(gca, 'XTickLabel', tmpTicks.^2);
set(gca, 'YTickLabel', tmpTicks.^2);
Pej_Add_IdentityDiag
end

function Add_Cloud_sqrt(N, P, Std, Color)
NC = 1000;
Nsqrt = sqrt(N);
sN = round((max(1, Nsqrt(randi(length(N), NC, 1))' .* (randn(NC,1)*.25+1))).^2);
sN(sN>max(N*1.2))=[];
mu = Pej_Transform_logit(P);
sR = Pej_rnd_BLN(sN, mu, Std.^2);
plot(sqrt(sR),sqrt(sN-sR), 'o', 'markersize', .5, 'color', Color, 'linewidth', .5, 'markerfacecolor', Color);
plot(sqrt([0 max(sN)*P]), sqrt([0 max(sN)*(1-P)]), '--','color', Color, 'linewidth', 1 );
end

function [TissuIdx, TissuLbl, DonorIdx, DonorLbl] = MakeDonorTissueFilt(MetaDFile, TissueXrf, InputFile)
Fin = fopen(InputFile, 'r');
ColLs(:,1) = strsplit(fgetl(Fin), char(9));
fclose(Fin);
SampleIDs = ColLs(5:end);
MetaData = Pej_Xref(SampleIDs, MetaDFile);
TissueID   = Pej_Xref(MetaData.SMTSD, TissueXrf);

for t = length(SampleIDs):-1:1
    tmpl = strsplit(SampleIDs{t}, '-');
    DonorID(t,1) = tmpl(2);
end
% Pej_Read_Table()
[TissuLbl, ~, TissuIdx] = unique(TissueID.tissue_abbrv);
[DonorLbl, ~, DonorIdx] = unique(DonorID);
end


function [GeneID, HapCounts] = GetCounts(S)
F1 = find(S(1:250)==char(9));
GeneID = S(F1(1)+1:F1(2)-1);
HapCounts = textscan(S(F1(4)+1:end), '%f|%f');
HapCounts = [HapCounts{1} HapCounts{2}];
end
