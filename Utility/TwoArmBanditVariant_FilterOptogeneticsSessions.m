if ~exist('DataHolder', 'var')
    disp('Warning: file loaded is not a DataHolder')
    % continue
end

nSessions = length(DataHolder);
OptoSession = {};
NonOptoSession = {};
for iSession = 1:nSessions
    SessionData = DataHolder{iSession};
    
    if SessionData.nTrials < 200
        fprintf('Session number %2.0f has less than 200 trials\n', iSession)
        continue
    end
    
    if isfield(SessionData.Custom.SessionMeta, 'OptoRemarks')
        OptoSession = [OptoSession, {SessionData}];
    else
        NonOptoSession = [NonOptoSession, {SessionData}];
    end
end