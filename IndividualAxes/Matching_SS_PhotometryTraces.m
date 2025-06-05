function AnalysisFigure = Matching_SS_PhotometryTraces(DataObject, FigureSize, AxeSize)
% for making an example of the task structure
% figure 'unit' set as 'inch' so that we know exactly the meters to pixels
% Developed by Antonio Lee @ BCCN Berlin
% Version 1.0 ~ April 2025

if nargin < 1
    [datafile, datapath] = uigetfile(OttLabDataServerFolderPath());
    load(fullfile(datapath, datafile));
elseif ischar(DataObject) || isstring(DataObject)
    load(DataObject);
else
    DataObject = DataObject;
end

try
    TE = DataObject.LoadProcessedPhotometryFile('TE');
catch
    disp('1st input is not a PhotometryData.mat')
    return
end

%% Load related data to local variabels
SplitedString = split(DataObject.SessionName, '_');
RatID = str2double(SplitedString{1});
if isnan(RatID)
    RatID = -1;
end
RatName = num2str(RatID);
% %%The following three lines doesn't not work, as the timestamp documented
% in the SessionData may not be the same as the one being used for saving
% Date = datestr(SessionData.Info.SessionDate, 'yyyymmdd');

SessionDateTime = strcat(SplitedString{3}, '_', SplitedString{4});

%% Initiatize figure
% create figure
if nargin < 2
    FigureSize = [   2,    5,  5.8, 2.0];
end
AnalysisFigure = figure('Position', FigureSize,...
                        'NumberTitle', 'off',...
                        'Name', strcat(RatName, '_', SessionDateTime, '_Matching'),...
                        'MenuBar', 'none',...
                        'Resize', 'off',...
                        'Unit', 'inch',...
                        'Color', 'none');

set(AnalysisFigure,...
    'Position', FigureSize,...
    'unit', 'inch');

FrameAxes = axes(AnalysisFigure, 'Position', [0 0 1 1]); % spacer for correct saving dimension
set(FrameAxes,...
    'XTick', [],...
    'YTick', [],...
    'XColor', 'none',...
    'YColor', 'none',...
    'Color', 'none')

% colour palette
ColourPalette = CommonColourPalette();

%% Photometry DataObject
GreenDemodData = DataObject.LoadProcessedPhotometryFile('Demod', 'Channel', 'Green');
Time = GreenDemodData.Time;
GreenDemodData = GreenDemodData.Data;

RedDemodData = DataObject.LoadProcessedPhotometryFile('Demod', 'Channel', 'Red');
RedDemodData = RedDemodData.Data;

ExampleIdx = 243:246;
AbsoluteTimeTrialStart = TE.AbsoluteTimeTrialStart(ExampleIdx);
AbsoluteTime = AbsoluteTimeTrialStart' + repmat(Time, length(ExampleIdx), 1);
AbsoluteTime = reshape(AbsoluteTime', 1, []);

GreenTraces = reshape(GreenDemodData(ExampleIdx, :)', 1, []);

RedTraces = reshape(RedDemodData(ExampleIdx, :)', 1, []);

%% Background DA
if nargin < 3
    AxeSize = [ 0.8, 0.55, 4.05,  0.90];
end

TracesAxes = axes(AnalysisFigure,...
                  'Position', AxeSize,...
                  'Units', 'inches',...
                  'Color', 'none');
set(TracesAxes,...
    'Position', AxeSize,...
    'Units', 'inches');
hold(TracesAxes, 'on');

yyaxis(TracesAxes, 'left')

GreenDemodDataPlot = plot(TracesAxes, AbsoluteTime, movmean(GreenTraces, 20),...
                          'LineStyle', '-',...
                          'Color', 'g',...
                          'LineWidth', 1);

set(TracesAxes,...
    'TickDir', 'out',...
    'YColor', 'g',...
    'YLim', [0.44, 0.52],...
    'FontSize', 12);
xlabel(TracesAxes, 'Time (s)')
ylabel(TracesAxes, 'GRAB-DA_{3h} (V)')

yyaxis(TracesAxes, 'right')

RedDemodDataPlot = plot(TracesAxes, AbsoluteTime, movmean(RedTraces, 20),...
                        'LineStyle', '-',...
                        'Color', 'r',...
                        'LineWidth', 1);

set(TracesAxes,...
    'YColor', 'r',...
    'YLim', [0.80, 0.90]);
ylabel(TracesAxes, 'tdTomato (V)')
title(TracesAxes, 'Example florescence traces')

exportgraphics(AnalysisFigure, strcat(RatName, '_', SessionDateTime, '_Matching_PhotometryTraces.pdf'), 'ContentType', 'vector', 'Resolution', 300)
end