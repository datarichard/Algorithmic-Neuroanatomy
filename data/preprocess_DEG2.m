% clear workspace
clear
clc

%% get list of current files
filenames = dir('111_degradation*');

for s = 1:length(filenames)
    
    filename = filenames(s).name;
    fid = fopen(filename);
    mydata = textscan(fid, '%f %*s %s','Delimiter','\t');
    time = mydata{1,1};
    event = mydata{1,2};
    clear fid mydata

    % adjust for trigger time:
    idx = strfind(event,'''pressed'': False, ''port'': 0, ''key'': 4, ''time'': 0');
    triggertime = find(not(cellfun('isempty',idx)));
    triggertime = time(triggertime(1));
    time = time - triggertime;
    clear idx

    %% Get the trial start and end times, and press times

    % get the trial start times:
    idx = strfind(event,'end rating B'); %returns a cell of 1 for every true event
    idxRow = find(not(cellfun('isempty',idx))); %returns an array of every row = 1 
    endRatingB = time(idxRow);
    clear idx idxRow

    idx = strfind(event,'begin trial');
    idxRow = find(not(cellfun('isempty',idx)));
    beginTrial = time(idxRow);
    trialStarts = [beginTrial;endRatingB];
    trialStarts = trialStarts(1:end-1);
    clear idx idxRow beginTrial

    % get the trial end times:
    idx = strfind(event,'start rating A');
    idxRow = find(not(cellfun('isempty',idx)));
    trialEnds = time(idxRow);
    clear idx idxRow

    %% get times of events of interest

    % find every left button press
    idx = strfind(event,'''pressed'': True, ''port'': 0, ''key'': 0, ''time'': 0'); 
    leftPressesIdx = find(not(cellfun('isempty',idx))); 
    leftPresses = time(leftPressesIdx);

    for b = 1:6
        LeftPressesBlock{b}=leftPresses(find(leftPresses>trialStarts(b)&leftPresses<trialEnds(b)));
    end
    
    % find every right button press
    idx = strfind(event,'''pressed'': True, ''port'': 0, ''key'': 1, ''time'': 0'); 
    rightPressesIdx = find(not(cellfun('isempty',idx))); 
    rightPresses = time(rightPressesIdx);
    
    for b = 1:6
        RightPressesBlock{b}=rightPresses(find(rightPresses>trialStarts(b)&rightPresses<trialEnds(b)));
    end

    % get the time of every earned left reward 'earn L':
    idx = strfind(event,'earn L'); %returns a cell of 1s for every true event
    idxRow = find(not(cellfun('isempty',idx))); %returns an array for every 1 row index
    earnLeftTimes = time(idxRow); %creates an array of times (from time array)

    for b = 1:6
        EarnLeftBlock{b}=earnLeftTimes(find(earnLeftTimes>trialStarts(b)&earnLeftTimes<trialEnds(b)));
    end

    % get the time of every earned left reward 'earn L':
    idx = strfind(event,'free L'); %returns a cell of 1s for every true event
    idxRow = find(not(cellfun('isempty',idx))); %returns an array for every 1 row index
    freeLeftTimes = time(idxRow); %creates an array of times (from time array)

    for b = 1:6
        FreeLeftBlock{b}=freeLeftTimes(find(freeLeftTimes>trialStarts(b)&freeLeftTimes<trialEnds(b)));
    end

    % get the time of every earned right reward 'earn R':
    idx = strfind(event,'earn R'); %returns a cell of 1s for every true event
    idxRow = find(not(cellfun('isempty',idx))); %returns an array for every 1 row index
    earnRightTimes = time(idxRow); %creates an array of times (from time array)

    for b = 1:6
        EarnRightBlock{b}=earnRightTimes(find(earnRightTimes>trialStarts(b)&earnRightTimes<trialEnds(b)));
    end

    % get the time of every free right reward 'free R':
    idx = strfind(event,'free R'); %returns a cell of 1s for every true event
    idxRow = find(not(cellfun('isempty',idx))); %returns an array for every 1 row index
    freeRightTimes = time(idxRow); %creates an array of times (from time array)

    for b = 1:6
        FreeRightBlock{b}=freeRightTimes(find(freeRightTimes>trialStarts(b)&freeRightTimes<trialEnds(b)));
    end

    %% create tables of events for each block (Table)
    
    o3 = zeros(0,3);
    o4 = zeros(0,3);
    
    for b = 1:6

        % table of 1) time, 2) A1/A2 3) Wait 4) O1/O2
        a1=LeftPressesBlock{b};  a1(:,2)=1; a1(:,4)=0;
        a2=RightPressesBlock{b}; a2(:,2)=1; a2(:,4)=0;
        o1=EarnLeftBlock{b};     o1(:,4)=1; 
        o2=EarnRightBlock{b};    o2(:,4)=1; 
        
        if isempty(FreeLeftBlock{b})
            o3 = zeros(0,4);
        else
            o3=FreeLeftBlock{b}; o3(:,4)=1; 
        end
        
        if isempty(FreeRightBlock{b})
            o4 = zeros(0,4);
        else
            o4=FreeRightBlock{b};o4(:,4)=1; 
        end
        
        % find wait seconds
        responselist = sortrows([LeftPressesBlock{b};RightPressesBlock{b}],1);
        responselist = fix(responselist);
        responselist = unique(responselist);
        
        starttime = fix(trialStarts(b)); endtime = round(trialEnds(b));
        wt = [];
        while starttime < endtime
            starttime = starttime + 1;
            if ~ismember(starttime,responselist)
                waitsec = starttime;
                wt(end+1,1) = waitsec;
            end           
        end
        
        wt(:,3) = 1; wt(:,4) = 0;
        
        
        LeftTable{b} =[a1;o1;o3;wt];    % 1)time 2)A1 3)wait 4)O1
        RightTable{b} =[a2;o2;o4;wt];   % 1)time 2)A2 3)wait 4)O2
        LeftTable{b} = sortrows(LeftTable{b},1);
        RightTable{b} = sortrows(RightTable{b},1);
        
    end
    
    % find start and end of all non-trial periods (startTimes,finishTimes)
    Nontrialtimes = [0;trialEnds]; % triggertime + start of rating period
    endtime = 2.602*350;
    finishTimes = [trialStarts;endtime];
    Nontrialdurations = finishTimes - Nontrialtimes;
    
    %% save data
    subject = filename(1:3);
    savefile = strcat(subject,'preprocessed');
    save(savefile,'Nontrialtimes','Nontrialdurations','LeftPressesBlock',...
        'RightPressesBlock','EarnLeftBlock','EarnRightBlock',...
        'FreeLeftBlock','FreeRightBlock','LeftTable','RightTable','trialStarts',...
        'trialEnds');
    %% stemplot events

    for b = 1:6

        if b==1||b==4||b==5
            degpresses = LeftTable{b}(LeftTable{b}(:,2)==1,1)-trialStarts(b);
            degrewards = LeftTable{b}(LeftTable{b}(:,4)==1,1)-trialStarts(b);
            conpresses = RightTable{b}(RightTable{b}(:,2)==1,1)-trialStarts(b);
            conrewards = RightTable{b}(RightTable{b}(:,4)==1,1)-trialStarts(b);
        else
            degpresses = RightTable{b}(RightTable{b}(:,2)==1,1)-trialStarts(b);
            degrewards = RightTable{b}(RightTable{b}(:,4)==1,1)-trialStarts(b);
            conpresses = LeftTable{b}(LeftTable{b}(:,2)==1,1)-trialStarts(b);
            conrewards = LeftTable{b}(LeftTable{b}(:,4)==1,1)-trialStarts(b);
        end
        
        y1 = ones(size(degpresses,1),1);
        y2 = ones(size(degrewards,1),1);
        y3 = ones(size(conpresses,1),1);
        y4 = ones(size(conrewards,1),1);
        
        figure(1)
        subplot(6,1,b)
        stem(degrewards,y2*1.2,'Marker','none','Color','r','LineWidth',1.25)
        hold on
        stem(degpresses,y1,'Marker','none')              
        hold off
        xlim([0 130])
        ylim([0 2])
        
        figure(2)
        subplot(6,1,b)
        stem(conrewards,y4*1.2,'Marker','none','Color','r','LineWidth',1.25)
        hold on
        stem(conpresses,y3,'Marker','none')        
        hold off
        xlim([0 130])
        ylim([0 2])
        
    end
    
end


        
     
     




