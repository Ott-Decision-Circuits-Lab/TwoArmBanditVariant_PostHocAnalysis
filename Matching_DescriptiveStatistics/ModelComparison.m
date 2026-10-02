Figure = figure();
set(Figure,...
    'Unit', 'inches',...
    'Position', [0.5, 0.5, 5.0, 5.2]);

Axes = axes(Figure);
set(Axes,...
    'Unit', 'inches',...
    'Position', [0.8, 1.2, 4.0, 3.8]);

hold on

YData = [BeliefState; QRL; LastC; NoMemory; Foraging; BetaDist]';

Plot = plot(Axes, [0:5], YData, 'Color', 'k', 'Marker', 'o');
set(Axes,...
    'FontSize', 12,...
    'TickDir', 'out',...
    'XLim', [-0.5, 5.5],...
    'XTick', [0:5],...
    'XTickLabel', {'Belief state', 'Q-RL', 'Last choice', 'No memory', 'Foraging', 'BetsDist'},...
    'XTickLabelRotation', 90,...
    'YLim', [0, 650])
ylabel(Axes,...
       '$-ln\mathcal{L}(\hat\theta, Y)$',...
       'Interpreter', 'latex')

[h,p,ci,stats] = ttest(BeliefState, QRL);
Text_0_1 = text(0.5, 475, sprintf('p=%4.2g', p));
set(Text_0_1, 'HorizontalAlignment', 'center')
Text_0_1_Line = plot([0, 1], [450, 450], 'k');

[h,p,ci,stats] = ttest(QRL, LastC);
Text_1_2 = text(1.5, 500, sprintf('p=%4.2g', p));
set(Text_1_2, 'HorizontalAlignment', 'center')
Text_1_2_Line = plot([1, 2], [475, 475], 'k');

[h,p,ci,stats] = ttest(LastC, NoMemory);
Text_2_3 = text(2.5, 525, sprintf('p=%4.2g', p));
set(Text_2_3, 'HorizontalAlignment', 'center')
Text_2_3_Line = plot([2, 3], [500, 500], 'k');

[h,p,ci,stats] = ttest(QRL, Foraging);
Text_1_4 = text(3.5, 575, sprintf('p=%4.2g', p));
set(Text_1_4, 'HorizontalAlignment', 'center')
Text_1_4_Line = plot([1, 4], [550, 550], 'k');

[h,p,ci,stats] = ttest(QRL, BetaDist);
Text_1_5 = text(4.5, 625, sprintf('p=%4.2g', p));
set(Text_1_5, 'HorizontalAlignment', 'center')
Text_1_5_Line = plot([1, 5], [600, 600], 'k');

% for AIC
nSessions = length(Models);
BeliefStateAIC = 2 * BeliefState + 2 * 6;
QRLAIC = 2 * QRL + 2 * 6;
LastCAIC = 2 * LastC + 2 * 5;
NoMemoryAIC = 2 * NoMemory + 2 * 4;
ForagingAIC = 2 * Foraging + 2 * 4;
BetaDistAIC = 2 * BetaDist + 2 * 6;

YData = [BeliefStateAIC', QRLAIC', LastCAIC', NoMemoryAIC', ForagingAIC', BetaDistAIC'];

for iSession = 1:nSessions
    set(Plot(iSession), 'YData', YData(iSession, :))
end

ylabel('AIC (a.u.)', 'Interpret', 'none')
set(Axes, 'YLim', [0, 1200])

[h,p,ci,stats] = ttest(BeliefStateAIC, QRLAIC);
set(Text_0_1,...
    'Position', [0.5, 975, 0],...
    'String', sprintf('p=%4.2g', p))
set(Text_0_1_Line,...
    'YData', [950, 950]);

[h,p,ci,stats] = ttest(QRLAIC, LastCAIC);
set(Text_1_2,...
    'Position', [1.5, 1000, 0],...
    'String', sprintf('p=%4.2g', p))
set(Text_1_2_Line,...
    'YData', [975, 975]);

[h,p,ci,stats] = ttest(LastCAIC, NoMemoryAIC);
set(Text_2_3,...
    'Position', [2.5, 1025, 0],...
    'String', sprintf('p=%4.2g', p))
set(Text_2_3_Line,...
    'YData', [1000, 1000]);

[h,p,ci,stats] = ttest(QRLAIC, ForagingAIC);
set(Text_1_4,...
    'Position', [3.5, 1075, 0],...
    'String', sprintf('p=%4.2g', p))
set(Text_1_4_Line,...
    'YData', [1050, 1050]);

[h,p,ci,stats] = ttest(QRLAIC, BetaDistAIC);
set(Text_1_5,...
    'Position', [4.5, 1125, 0],...
    'String', sprintf('p=%4.2g', p))
set(Text_1_5_Line,...
    'YData', [1100, 1100]);