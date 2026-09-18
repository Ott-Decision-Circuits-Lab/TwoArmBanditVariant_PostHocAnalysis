function String = SelectDateAndRemarksFromDataHolder(DataHolder)

nSessions = length(DataHolder);
String = {};
for iSession = 1:nSessions
    SessionData = DataHolder{iSession};

    SessionDatetime = datetime(SessionData.Custom.General.SessionDate);
    String{iSession, 1} = char(SessionDatetime, 'yyyyMMdd');
    
    try
        SessionRemarks = SessionData.Custom.SessionMeta.BehaviouralRemarks;
    catch
        SessionRemarks = '';
    end
    String{iSession, 2} = SessionRemarks;

end % for-loop

end % function

