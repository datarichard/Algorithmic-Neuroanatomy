function [Qcon,Qdeg] = Kalman9a_Qvalues(filename)
% based on Kalman_sim2

% get free parameters
subj = filename(1:3);
paramfile = ['fmin',subj,'.mat'];
load(paramfile,'param');
slope = -inf;
tangent = param(1);
v = param(2);       % observation uncertainty 
d = param(3);       % process uncertainty 

load(filename); % containing LeftTable, RightTable

Q = []; % table for Q-values created below

for block = 1:6
    %% create table of inputs and outputs
    % table: [time, a1, a2, w, o1, o2]
    table = LeftTable{block};
    table(:,4:5) = table(:,3:4); table(:,3) = 0; table(:,6) = 0;
    temp = RightTable{block};
    temp(:,3:5) = temp(:,2:4);
    temp(:,2) = 0; temp(:,6) = temp(:,5); temp(:,5) = 0;
    table = sortrows([table;temp]); table = unique(table,'rows');
    table(:,7:8) = table(:,5:6); % add columns 7 and 8 for reward traces (outputs)

    %% reverse trace reward times
    %slope = -1000; tangent = 1.05;

    for action = 1:2

        if action == 1; col1=7;col2=2; else col1=8;col2=3; end

        rewards = table(table(:,col1)==1,1); responses = table(table(:,col2)==1,1);
        for r = 1:length(rewards)
            reward = rewards(r);
            response = responses(find(responses<reward,1,'last'));
            if ~isempty(response)
                delay = reward - response;
                trace = sigmf(delay,[slope,tangent]); % sigmoid function
                %trace = exp(slope.*(delay)); % exponential with negative slope
                responserow = find(table(:,1)==response);
                table(responserow,col1) = table(responserow,col1) + trace;
                rewardrow = find(table(:,1)==reward);
                table(rewardrow,col1) = table(rewardrow,col1) - trace;
            end
        end

    end


    %% define inputs and outputs
    ntrial = size(table,1);
    inputs = table(:,2:4); inputs(:,3) = 1; % [a1, a2, b]
    outputs = table(:,7:8);                 % [o1, o2]
    
    % create transition matrix
    M = zeros(3,2)/3; C = eye(3)*1;
    matrix = reshape(M,[1 6]); vmatrix = reshape(C,[1 9]);
    
    for trial = 1:size(table,1)
        
        %% Kalman filter
        i = inputs(trial,:)';           % inputs are a column vector
        o = outputs(trial,:);           % outputs are a row vector

        alpha = (v + i'*C*i);           % total uncertainty (outcome + expected)
        PE = C*i*(o' - M'*i)';          % prediction-error
        delta_mu = PE*1./alpha;         % calculate delta mu
        M = M + delta_mu;               % update mu
        matrix(end+1,:) = reshape(M,[1,6]);

        cmatrix = C*(i*i')*C;           % adjust uncertainty by observed outcome?
        delta_C = 1/alpha*cmatrix;      % uncertainty reduction by observation
        C = C - delta_C+eye(3)*d;       % update C (reduce uncertainty)
        vmatrix(end+1,:) = reshape(C,[1 9]);

        % M is a 3 by 2 matrix:
        % 
        % A1-O1   A1-O2   
        %    
        % A2-O1   A2-O2   
        % 
        % B-O1    B-O2    
        %
        % transitions: 
        % col1    col2    col3    col4    col5    col6    
        % A1-O1   A2-O1   B-O1    A1-O2   A2-O2   B-O2 
    end
    
    %% find actions and calculate values
    % q = [time, left, right, mu1, mu2]
    q = [table(:,1:3),matrix(1:end-1,1),matrix(1:end-1,5)];
    
    % q = [time, left, right, mu1, mu2, P(left), P(right)]
    var1 = vmatrix(1:end-1,1);
    var2 = vmatrix(1:end-1,5);
    q(:,6) = normcdf(0,q(:,4),sqrt(var1),'upper'); % P(left>0)
    q(:,7) = normcdf(0,q(:,5),sqrt(var2),'upper'); % P(right>0)
    
    % q = [time, left, right, mu1, mu2, P(left), P(right), P(left=best)]
    mu_diff = q(:,4) - q(:,5); 
    var_sum = var1+var2; sigma = sqrt(var_sum);
    q(:,8) = normcdf(0,mu_diff,sigma,'upper'); % P(left=best)
    
    % q = [time, left, right, mu1, mu2, P(left), P(right), P(best), block]
    q(:,9) = block; % add block
    q(:,1) = q(:,1)-q(1,1); % adjust event times
  
    Q = [Q;q];
    
end % end of block

%% %%%%%%%% Sort degraded and contingent values %%%%%%%%%%%%%%%%%%%%%%%
Q = Q(Q(:,2)+Q(:,3)>0,:); % select response rows

% find degraded left response rows
LeftDegradedBlocks = Q(:,end)==1|Q(:,end)==4|Q(:,end)==5; % index left = degraded
LeftResponses = Q(:,2) == 1;
LeftDegradedRows = and(LeftDegradedBlocks,LeftResponses);

% find degraded right response rows
RightDegradedBlocks = ~LeftDegradedBlocks;
RightResponses = Q(:,3) == 1;
RightDegradedRows = and(RightDegradedBlocks,RightResponses);

% combine and concatenate degraded rows
DegradedRows = or(RightDegradedRows,LeftDegradedRows);
Qdeg = Q(DegradedRows,:);

Qdeg(:,4) = Qdeg(:,4).*Qdeg(:,2); 
Qdeg(:,5) = Qdeg(:,5).*Qdeg(:,3);
Qdeg(:,4) = Qdeg(:,4) + Qdeg(:,5);

Qdeg(:,6) = Qdeg(:,6).*Qdeg(:,2); 
Qdeg(:,7) = Qdeg(:,7).*Qdeg(:,3);
Qdeg(:,5) = Qdeg(:,6) + Qdeg(:,7);

% Q_deg = [blocktime,block,right=1,mu,P(causal),P(best)]
Qdeg(:,6) = Qdeg(:,8); Qdeg(:,2) = Qdeg(:,9);
Qdeg(:,9) = []; Qdeg(:,8) = []; Qdeg(:,7) = [];

% find contingent right response rows
RightContingentRows = and(LeftDegradedBlocks,RightResponses);

% find contingent left response rows
LeftContingentRows = and(RightDegradedBlocks,LeftResponses);

% combine and concatenate contingent rows
ContingentRows = or(LeftContingentRows,RightContingentRows);
Qcon = Q(ContingentRows,:);

Qcon(:,4) = Qcon(:,4).*Qcon(:,2); 
Qcon(:,5) = Qcon(:,5).*Qcon(:,3);
Qcon(:,4) = Qcon(:,4) + Qcon(:,5);

Qcon(:,6) = Qcon(:,6).*Qcon(:,2); 
Qcon(:,7) = Qcon(:,7).*Qcon(:,3);
Qcon(:,5) = Qcon(:,6) + Qcon(:,7);

% Q_con = [blocktime,block,right=1,mu,P(causal),P(best)]
Qcon(:,6) = Qcon(:,8); Qcon(:,2) = Qcon(:,9);
Qcon(:,9) = []; Qcon(:,8) = []; Qcon(:,7) = [];

end

