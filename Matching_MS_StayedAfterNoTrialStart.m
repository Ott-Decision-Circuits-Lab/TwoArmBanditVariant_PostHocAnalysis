

for iNoTrialStart = 5

    String = replace(num2str([0, ones(1, iNoTrialStart), 0]), ' ', '');
    
    Stayed = [];
    RewardedStayed = [];
    UnrewardedStayed = [];
    for iSession = 1:length(DataHolder)
        SessionData = DataHolder{iSession};

        nTrials = SessionData.nTrials;
        TrialData = SessionData.Custom.TrialData;
        
        ChoiceLeft = TrialData.ChoiceLeft(1:nTrials);
        Rewarded = TrialData.Rewarded(1:nTrials);
        NoTrialStart = TrialData.NoTrialStart(1:nTrials);
        BrokeFixation = TrialData.BrokeFixation(1:nTrials);
        
        NoTrialStartChar = replace(num2str(NoTrialStart == 1 & isnan(BrokeFixation)), ' ', '');
        NoChoiceChar = replace(num2str(isnan(ChoiceLeft)), ' ', '');
    
        Index1 = strfind(NoTrialStartChar, String); % may include bounded by NoDecision
        Index2 = strfind(NoChoiceChar, String); % may include NoDecision
    
        Index = intersect(Index1, Index2);

        for iIndex = Index
            Stayed = [Stayed, ChoiceLeft(iIndex) == ChoiceLeft(iIndex + iNoTrialStart + 1)];
            
            if Rewarded(iIndex) == 1
                RewardedStayed = [RewardedStayed, ChoiceLeft(iIndex) == ChoiceLeft(iIndex + iNoTrialStart + 1)];
            elseif Rewarded(iIndex) == 0
                UnrewardedStayed = [UnrewardedStayed, ChoiceLeft(iIndex) == ChoiceLeft(iIndex + iNoTrialStart + 1)];
            end
        end
    end
end
