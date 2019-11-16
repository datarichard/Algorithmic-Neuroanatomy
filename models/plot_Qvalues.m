clear
clc
%%
y_con = []; y_deg = [];
fileNames  = dir ('108*.mat');

%%
for s = 1:size(fileNames,1)
    filename = fileNames(s).name(1:19);
    subjectnum = str2num(filename(1:3));
    [Qcon,Qdeg] = Kalman9a_Qvalues(filename);
    
    % Qcon = [blocktime,block,right=1,mu,P(causal),P(best),subjecnum]
    Qcon(:,end+1) = subjectnum; Qdeg(:,end+1) = subjectnum;
    
    % get mu
    Mu_con = Qcon(:,1); Mu_con(:,2)=Qcon(:,4); Mu_con(:,3)=Qcon(:,7);
    con = sortrows(Mu_con);
    con(:,1) = fix(con(:,1));
    con(:,1) = con(:,1)+1;
    y1 = accumarray(con(:,1),con(:,2),[],@max);
    
    Mu_deg =  Qdeg(:,1); Mu_deg(:,2)=Qdeg(:,4); Mu_deg(:,3)=Qdeg(:,7);
    deg = sortrows(Mu_deg);
    deg(:,1) = fix(deg(:,1));
    deg(:,1) = deg(:,1)+1;
    y2 = accumarray(deg(:,1),deg(:,2),[],@max);
    
    % get P(causal)
    Causal_con = Qcon(:,1); Causal_con(:,2)=Qcon(:,5); Causal_con(:,3)=Qcon(:,7);
%     con = sortrows(Causal_con);
%     con(:,1) = fix(con(:,1));
%     con(:,1) = con(:,1)+1;
%     y1 = accumarray(con(:,1),con(:,2),[],@mean);
    
    Causal_deg =  Qdeg(:,1); Causal_deg(:,2)=Qdeg(:,5); Causal_deg(:,3)=Qdeg(:,7);
%     deg = sortrows(Causal_deg);
%     deg(:,1) = fix(deg(:,1));
%     deg(:,1) = deg(:,1)+1;
%     y2 = accumarray(deg(:,1),deg(:,2),[],@mean);
    
    % get P(best)
    Best_con = Qcon(:,1); Best_con(:,2)=Qcon(:,6); Best_con(:,3)=Qcon(:,7);
%     con = sortrows(Causal_con);
%     con(:,1) = fix(con(:,1));
%     con(:,1) = con(:,1)+1;
%     y1 = accumarray(con(:,1),con(:,2),[],@mean);
    
    Best_deg =  Qdeg(:,1); Best_deg(:,2)=Qdeg(:,6); Best_deg(:,3)=Qdeg(:,7);
%     deg = sortrows(Causal_deg);
%     deg(:,1) = fix(deg(:,1));
%     deg(:,1) = deg(:,1)+1;
%     y2 = accumarray(deg(:,1),deg(:,2),[],@mean);

    plot(y1); hold on; plot(y2); hold off;
    
    % add to group data
    % y_con(:,end+1) = y1(1:120); y_deg(:,end+1) = y2(1:120);
    
end
