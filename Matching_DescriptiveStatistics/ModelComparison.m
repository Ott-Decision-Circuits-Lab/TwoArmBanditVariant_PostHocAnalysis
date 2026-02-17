Figure = figure();
set(Figure,...
    'Unit', 'inches',...
    'Position', [0.5, 0.5, 3, 3.5]);

Axes = axes(Figure);
set(Axes,...
    'Unit', 'inches',...
    'Position', [0.8, 0.5, 2, 2.5]);

YData = [ForagingRLMinNegLogDataLikelihood; QRLMinNegLogDataLikelihood]';

Plot = plot(Axes, [0, 1], YData, 'Color', 'k', 'Marker', 'o');
set(Axes,...
    'XLim', [-0.5, 1.5],...
    'XTick', [0, 1],...
    'XTickLabel', {'last memory', 'Q-RL'})
ylabel(Axes,...
       '$-ln\mathcal{L}(\hat\theta, Y)$',...
       'Interpreter', 'latex')

[h,p,ci,stats] = ttest(ForagingRLMinNegLogDataLikelihood, QRLMinNegLogDataLikelihood);
SignificanceText = text(0.5, 550, sprintf('p=%4.2g', p));
set(SignificanceText, 'HorizontalAlignment', 'center')

% for AIC
nSessions = length(Models);
ForagingAIC = -2 * ForagingRLMinNegLogDataLikelihood + 2 * 5;
QRLAIC = -2 * QRLMinNegLogDataLikelihood + 2 * 6;
YData = [ForagingAIC', QRLAIC'];

for iSession = 1:nSessions
    set(Plot(iSession), 'YData', YData(iSession, :))
end

set(Axes, 'YLim', [-1200, 0])
set(SignificanceText, 'Position', [0.5, -1100, 0])

ylabel('AIC (a.u.)', 'Interpret', 'none')

[h,p,ci,stats] = ttest(ForagingAIC, QRLAIC);
set(SignificanceText, 'String', sprintf('p=%4.2g', p))
