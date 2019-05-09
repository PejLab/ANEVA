function Add_Roc(Target, Score, Label, LineSpec)
F = ~isnan(Target) & ~isnan(Score);
[X, Y , ~, AUC] = perfcurve(Target(F), Score(F), true);
MyL = [Label ':' num2str(AUC*100, 3) '%'];
plot(X,Y, LineSpec, 'DisplayName', MyL, 'linewidth', 1.2);
% plot([0 1], [0 1]);
ylabel('True positive rate')
xlabel('False positive rate')
box on
set(gcf, 'position', [1 1 200 175]*1.5)

end