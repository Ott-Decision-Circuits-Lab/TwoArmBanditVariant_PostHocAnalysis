function TwoArmBanditVariant_PoolingData(varargin)
%{
First create on 20230622 by Antonio Lee for AG Ott @HU Berlin
With the file server architecture, this functions runs through the
session folders of the designated Animals after a PeriodStartDate and before a
PeriodEndDate to parse out the TrialData necessary for Analysis.m
%}

%% obtain input and check if they are of the correct format
p = inputParser;
addParameter(p, 'Animal', {}, @iscell); % Cell array of numerics, holder for rat_id
addParameter(p, 'PeriodStartDate', {}, @iscell); % Cell array of char, date in the format of yyyyMMdd
addParameter(p, 'PeriodEndDate', {}, @iscell);
%{ Not Implemented
% addParameter(p, 'FilterOut', {}, @iscell); % Cell array of char
% addParameter(p, 'SeparateByAnimal', true, @islogical); % if true, the Analysis is done per Animal; if false, the Analysis is pooled
% addParameter(p, 'SeparateBySession', false, @islogical); % if true, the Analysis is done per Session; if false, the Analysis is pooled
%}

parse(p, varargin{:});
p = p.Results;

FileServerDataFolderPath = OttLabDataServerFolderPath();

if isempty(p.Animal)
    AnimalDataFolderPath = uigetdir(FileServerDataFolderPath);
    [~, RatID, ~] = fileparts(AnimalDataFolderPath);
    p.Animal = {str2double(RatID)};
end

if isempty(p.PeriodStartDate) || isempty(p.PeriodEndDate)
    error('Period is not defined')
end

if length(p.PeriodStartDate) ~= length(p.PeriodEndDate)
    error('Miscatch in length of PeriodStartDate and PeriodEndDate')
end

for iPeriodStartDate = 1:length(p.PeriodStartDate)
    PeriodStartDate = p.PeriodStartDate{iPeriodStartDate};
    PeriodStartDate = datetime(PeriodStartDate, 'Inputformat', 'yyyyMMdd');

    PeriodEndDate = p.PeriodEndDate{iPeriodStartDate};
    PeriodEndDate = datetime(PeriodEndDate, 'Inputformat', 'yyyyMMdd');   
end

if length(p.PeriodStartDate) ~= 1 && length(p.PeriodStartDate) ~= length(p.Animal)
    error('Length of PeriodStartDate is not applicable to the length of Animals')
end

if length(p.Animal) > length(p.PeriodStartDate)
    p.PeriodStartDate = repmat(p.PeriodStartDate, [1, length(p.Animal)]);
    p.PeriodEndDate = repmat(p.PeriodEndDate, [1, length(p.Animal)]);
end

%% start looping to select SessionData as a cell array
DataHolder = {};
for iAnimal = 1:length(p.Animal)
    RatID = p.Animal{iAnimal};
    if ~isnumeric(RatID)
        disp(['The input ', num2str(RatID), ' is not a numeric. It will be ignored.'])
        continue
    end

    if ~isfolder([FileServerDataFolderPath, num2str(RatID)])
        disp(['Folder for rat_id as ', num2str(RatID), ' is not found on the file server data folder. It will be ignored.'])
        continue
    end
    
    PeriodStartDate = datetime(p.PeriodStartDate{iAnimal}, 'InputFormat', 'yyyyMMdd');
    PeriodEndDate = datetime(p.PeriodEndDate{iAnimal}, 'InputFormat', 'yyyyMMdd');
    
    BpodSessionsFolderPath = [FileServerDataFolderPath, num2str(RatID), '\bpod_session\'];
    SessionDataFolders = dir(BpodSessionsFolderPath);
    for iSession = 1:length(SessionDataFolders)
        if ~SessionDataFolders(iSession).isdir
            continue
        end

        SessionDate = '';
        try
            SessionDate = datetime(SessionDataFolders(iSession).name, 'InputFormat', 'yyyyMMdd_HHmmss');
        catch
        end

        if isempty(SessionDate)
            continue
        end
        
        if SessionDate < PeriodStartDate || SessionDate > PeriodEndDate
            continue
        end
        
        SessionDataFolderPath = [BpodSessionsFolderPath, SessionDataFolders(iSession).name, '\'];
        SessionData = dir(SessionDataFolderPath);
        SessionFilePath = '';
        for iData = 1:length(SessionData)
            if length(SessionData(iData).name) < 19
                continue
            end

            if strcmpi(SessionData(iData).name(end-18:end), [datestr(SessionDate, 'YYYYmmDD_hhMMss'), '.mat'])
                SessionFilePath = [SessionDataFolderPath, SessionData(iData).name];
                break
            end
        end

        if isempty(SessionFilePath)
            disp(['No SessionFile in folder ', SessionDataFolders(iSession).name, ' is found.'])
            continue
        end

        load(SessionFilePath); % variable name will be 'SessionData'
        if SessionData.nTrials >= 50 && SessionData.SettingsFile.GUI.CatchTrial && SessionData.SettingsFile.GUI.FeedbackDelayGrace >= 0.2 && SessionData.SettingsFile.GUI.FeedbackDelayMax == 8
            DataHolder{end+1} = SessionData;
        end
    end
    
    SelectedDataFolderPath = [FileServerDataFolderPath, num2str(RatID), '\bpod_session_pooled_analysis\',...
                              num2str(RatID), '_TwoArmBanditVariant_', p.PeriodStartDate{iAnimal}, '_', p.PeriodEndDate{iAnimal}, '\'];
    
    if ~isfolder(SelectedDataFolderPath)
        mkdir(SelectedDataFolderPath)
    end

    SelectedDataFilePath = fullfile(SelectedDataFolderPath, '\Selected_Data.mat');
    save(SelectedDataFilePath, 'DataHolder')
    
    %% Concatenate SessionData into one single file
    ConcatenatedDataHolder.Info.Subject = num2str(RatID);
    ConcatenatedDataHolder.Info.SessionDate = datetime('20000101', 'InputFormat', 'yyyyMMdd'); % temporary entry, later corrected
    ConcatenatedDataHolder.nTrials = 0;
    
    ConcatenatedTrialData.ChoiceLeft = [];
    ConcatenatedTrialData.Baited = [];
    ConcatenatedTrialData.IncorrectChoice = [];
    ConcatenatedTrialData.NoDecision = [];
    ConcatenatedTrialData.NoTrialStart = [];
    ConcatenatedTrialData.BrokeFixation = [];
    ConcatenatedTrialData.EarlyWithdrawal = [];
    ConcatenatedTrialData.StartNewTrial = [];
    ConcatenatedTrialData.SkippedFeedback = [];
    ConcatenatedTrialData.Rewarded = [];
    
    ConcatenatedTrialData.SampleTime = [];
    ConcatenatedTrialData.MoveTime = [];
    ConcatenatedTrialData.FeedbackWaitingTime = [];
    ConcatenatedTrialData.RewardProb = [];
    ConcatenatedTrialData.LightLeft = [];
    
    ConcatenatedTrialData.BlockNumber = [];
    ConcatenatedTrialData.BlockTrialNumber = [];
    ConcatenatedTrialData.RewardMagnitude = [];

    % for files before April 2023, no DrinkingTime is available
    ConcatenatedTrialData.DrinkingTime = [];
    
    ConcatenatedDataHolder.RawEvents.Trial = {};
    ConcatenatedDataHolder.TrialStartTimestamp = [];
    ConcatenatedDataHolder.TrialEndTimestamp = [];
    for iSession = 1:length(DataHolder)
        nTrials = DataHolder{iSession}.nTrials;
        ConcatenatedDataHolder.nTrials = ConcatenatedDataHolder.nTrials + nTrials;
        
        TrialData = DataHolder{iSession}.Custom.TrialData;
        
        ConcatenatedTrialData.ChoiceLeft = [ConcatenatedTrialData.ChoiceLeft, TrialData.ChoiceLeft(1:nTrials)];
        ConcatenatedTrialData.Baited = [ConcatenatedTrialData.Baited, TrialData.Baited(:, 1:nTrials)];
        ConcatenatedTrialData.IncorrectChoice = [ConcatenatedTrialData.IncorrectChoice, TrialData.IncorrectChoice(1:nTrials)];
        ConcatenatedTrialData.NoDecision = [ConcatenatedTrialData.NoDecision, TrialData.NoDecision(1:nTrials)];
        ConcatenatedTrialData.NoTrialStart = [ConcatenatedTrialData.NoTrialStart, TrialData.NoTrialStart(1:nTrials)];
        ConcatenatedTrialData.BrokeFixation = [ConcatenatedTrialData.BrokeFixation, TrialData.BrokeFixation(1:nTrials)];
        ConcatenatedTrialData.EarlyWithdrawal = [ConcatenatedTrialData.EarlyWithdrawal, TrialData.EarlyWithdrawal(1:nTrials)];
        ConcatenatedTrialData.StartNewTrial = [ConcatenatedTrialData.StartNewTrial, TrialData.StartNewTrial(1:nTrials)];
        ConcatenatedTrialData.SkippedFeedback = [ConcatenatedTrialData.SkippedFeedback, TrialData.SkippedFeedback(1:nTrials)];
        ConcatenatedTrialData.Rewarded = [ConcatenatedTrialData.Rewarded, TrialData.Rewarded(1:nTrials)];
        
        ConcatenatedTrialData.SampleTime = [ConcatenatedTrialData.SampleTime, TrialData.SampleTime(1:nTrials)];
        ConcatenatedTrialData.MoveTime = [ConcatenatedTrialData.MoveTime, TrialData.MoveTime(1:nTrials)];
        ConcatenatedTrialData.FeedbackWaitingTime = [ConcatenatedTrialData.FeedbackWaitingTime, TrialData.FeedbackWaitingTime(1:nTrials)];
        ConcatenatedTrialData.RewardProb = [ConcatenatedTrialData.RewardProb, TrialData.RewardProb(:, 1:nTrials)];
        ConcatenatedTrialData.LightLeft = [ConcatenatedTrialData.LightLeft, TrialData.LightLeft(1:nTrials)];
        
        ConcatenatedTrialData.BlockNumber = [ConcatenatedTrialData.BlockNumber, TrialData.BlockNumber(:, 1:nTrials)];
        ConcatenatedTrialData.BlockTrialNumber = [ConcatenatedTrialData.BlockTrialNumber, TrialData.BlockTrialNumber(:, 1:nTrials)];
        ConcatenatedTrialData.RewardMagnitude = [ConcatenatedTrialData.RewardMagnitude, TrialData.RewardMagnitude(:, 1:nTrials)];
        
        % for files before April 2023, no DrinkingTime is available
        try
            ConcatenatedTrialData.DrinkingTime = [ConcatenatedTrialData.DrinkingTime, TrialData.DrinkingTime(1:nTrials)];
        catch
            ConcatenatedTrialData.DrinkingTime = [ConcatenatedTrialData.DrinkingTime, nan(1, nTrials)];
        end
        
        ConcatenatedDataHolder.RawEvents.Trial = [ConcatenatedDataHolder.RawEvents.Trial, DataHolder{iSession}.RawEvents.Trial];
        ConcatenatedDataHolder.TrialStartTimestamp = [ConcatenatedDataHolder.TrialStartTimestamp, DataHolder{iSession}.TrialStartTimestamp(:,1:nTrials)];
        ConcatenatedDataHolder.TrialEndTimestamp = [ConcatenatedDataHolder.TrialEndTimestamp, DataHolder{iSession}.TrialEndTimestamp(:,1:nTrials)];
    end
    ConcatenatedDataHolder.Custom.TrialData = ConcatenatedTrialData;
    ConcatenatedDataHolder.SettingsFile = DataHolder{iSession}.SettingsFile;
    
    SessionData = ConcatenatedDataHolder;
    ConcatenatedDataFilePath = fullfile(SelectedDataFolderPath, 'Concatenated_Data.mat');
    save(ConcatenatedDataFilePath, 'SessionData')

    % FigHandle = Analysis(ConcatenatedDataFilePath);
    % saveas(FigHandle, fullfile(SelectedDataFolderPath, 'Analysis.fig'))
    % saveas(FigHandle, fullfile(SelectedDataFolderPath, 'Analysis.png'))
    
end % iAnimal

end % function