function [Q] = Kalman_converge_dP(filename)
% previously v8a

% get free parameters
subj = filename(1:3);
paramfile = ['fmin',subj,'.mat'];
load(paramfile,'param');
slope = param(1);   % sigmoid decay parameter 0 to -inf (e.g., -10)
tangent = param(2); % sigmoid decay parameter 0 to 100 (e.g., 0.5)
v = param(3);       % observation uncertainty 0.001 to 1 (e.g., 0.5)
d = param(4);       % process uncertainty 0 to 1 (e.g., 0.5)
%tau = param(5);     % exploration/exploitation 0 to 10 (e.g., 2)

load(filename); % containing LeftTable, RightTable

Q = []; % table for Q-values created below

for b = 1:6 % for each block...

    for r = 1:2 % for left and right actions separately...

        %% ... create a table of events and event traces.
        if r == 1
            table = LeftTable{b};
            alttable = RightTable{b};
        else
            table = RightTable{b};
            alttable = LeftTable{b};
        end
        
        % remove duplicate entries
        table = unique(table,'rows');
        alttable = unique(alttable,'rows');
        
        % Get the response times.
        responseTable = table((table(:,2)==1),1); responseTable(:,2) = 1;
        waitTable = alttable((alttable(:,2)==1),1); waitTable(:,3) = 1;
        rewardTable = table((table(:,3)==1),1); rewardTable(:,4) = 1; 
        
        % Recreate table of responses, waits and rewards
        responseTable(:,3:4) = 0; waitTable(:,4) = 0;
        table = sortrows([responseTable;waitTable;rewardTable],1);

        % 2. Calculate the trace of each reward:
        for i = 1:size(rewardTable,1)

            % 2.1 Get current reward time
            current = rewardTable(i,1);

            % 2.2 Get previous response time
            last = responseTable(find(responseTable(:,1)<current,1,'last'));

            % 2.3 Calculate reverse reward trace
            if ~isempty(last) % if a previous response exists...

                % 2.3.1 Get the delay time 
                delay = current - last;

                % Calculate trace using decreasing function
                trace = sigmf(delay,[slope,tangent]); % sigmoidal function
                %trace = exp(slope.*(delay)); % exponential function

                % And add trace to table(:,4) when action is present
                responserow = find(table(:,1) == last);
                table(responserow,4) = table(responserow,4)+trace;
                
                % subtract trace from table(:,4) at reward time
                rewardrow = find(table(:,1) == current);
                table(rewardrow,4) = table(rewardrow,4) - trace;
                            
            end            
        end
        
        % Calculate context trace
        table(:,3) = 1;

        %% kalman filter loop 
        % starting values
        mu = [0,0]; % cue expected values (beliefs)
        C = eye(2); % cue covariance matrix (uncertainty matrix)
        variance = [1,1]; 

        for i = 1:size(table,1)

            c = table(i,2:3)';                  % set cue values (transposed)
            t = table(i,4);                     % set outcome values

            alpha = (v + c'*C*c);               % calculate learning rate
            PE = (t - mu(i,:)*c)*C*c;           % calculate prediction-error
            delta_mu = alpha^-1*PE;             % calculate delta mu
            new_mu = mu(i,:) + delta_mu';       % update mu
            mu(end+1,:) = new_mu;               % store mu

            cmatrix = C*(c*c')*C;               % intermediate step to calculate C
            delta_C = -alpha^-1*cmatrix;        % calculate delta_C
            C = C + delta_C+eye(2)*d;           % update C (reduce uncertainty)
            variance(end+1,1:4) = [C(1,1),C(2,2),C(1,2),C(2,1)]; % store belief uncertainty
            
        end
        
        %% Delta P 
        if r == 1
            earned = size(EarnLeftBlock{b},1);
            free = size(FreeLeftBlock{b},1);
            responses = size(LeftPressesBlock{b},1);
%             waits = size(RightPressesBlock{b},1);
        else
            earned = size(EarnRightBlock{b},1);
            free = size(FreeRightBlock{b},1);
            responses = size(RightPressesBlock{b},1);
%             waits = size(LeftPressesBlock{b},1);
        end
        
        waits = 120 - size(unique(round([LeftPressesBlock{b};RightPressesBlock{b}])),1);
        
        posP = earned./responses; negP = free./waits;
        
%         if negP < 1
%             dP = posP - negP;
%         else
%             waits = unique(round([LeftPressesBlock{b};RightPressesBlock{b}]));
%             dP = posP - free./waits;
%         end

        dP = posP - negP;

        
        %% find actions and calculate values
        rows = find(table(:,2)==1); % find the row index of responses
        mu1 = mu(rows,1); mu2 = mu(rows,2); mu_diff = mu1 - mu2;
        var1 = variance(rows,1); var2 = variance(rows,2);
        var_diff = var1 + var2; sigma = sqrt(var_diff);
        Paction = normcdf(0,mu_diff,sigma,'upper');
        
        % find trial start time
        starttime = Nontrialtimes(b) + Nontrialdurations(b);
        % create action value table q = [time, qleft, qright, response]
        if isempty(rows)
            q{r} = zeros(1,12);
        else
            q{r}(:,1) = table(rows,1)-starttime;  % times adusted for start
            q{r}(:,2) = b;              % block
            q{r}(:,3) = r;              % response
            q{r}(:,4) = mu1;              
            q{r}(:,5) = mu2;              
            q{r}(:,6) = mu_diff;
            q{r}(:,7) = var1;
            q{r}(:,8) = var2;
            q{r}(:,9) = sigma;
            q{r}(:,10) = Paction;
            q{r}(:,11) = posP; % positive contingency
            q{r}(:,12) = dP;
        end
        
    end
    
    q{3} = [q{1};q{2}]; 
    q{3} = sortrows(q{3},1);
    Q(end+1:end+length(q{3}),1:12) = q{3};
    clear q
    
end

%% clear zero time values
Q(Q(:,1)==0,:) = [];

%% identify degraded responses ( == 1)
LeftDegradedBlocks = Q(:,2)==1|Q(:,2)==4|Q(:,2)==5; % blocks 1, 4 and 5 = 1
LeftResponses = Q(:,3) == 1;                          % left responses = 1
LeftDegradedResponses = LeftDegradedBlocks.*LeftResponses; % degraded = 1

RightDegradedBlocks = 1 - LeftDegradedBlocks;
RightResponses = Q(:,3) == 2;
RightDegradedResponses = RightDegradedBlocks.*RightResponses;

DegradedResponses = LeftDegradedResponses + RightDegradedResponses;
Q(:,end+1) = DegradedResponses;

%% plot results
% x = Q(:,1); y1 = Q(:,2); y2 = Q(:,3);
% y1 = y1 - mean(y1); y2 = y2 - mean(y2);
% plot(x,y1); hold on; plot(x,y2,'r'); hold off

end

