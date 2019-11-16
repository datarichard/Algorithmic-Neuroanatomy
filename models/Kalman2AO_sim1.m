clear
clc

% set free parameters
ntrial = 120;                 % number of trials (10000)
contingency = [0.2 0.2 0.2];  % AO contingencies (0.5 0.5 0.5)
v = 5;                        % learning rate (100)
tau = 10;

% define inputs and outputs
inputs = zeros(ntrial,3);       % [a1, a2, b]
outputs = zeros(ntrial,2);      % [o1, o2]
rewards = zeros(ntrial,3);

% define priors
M = zeros(3,2); % belief matrix
C = eye(3);     % covariance (transition) matrix

% store values
beliefs = reshape(M,[1,6]); covariance = [];
deltamu = [];
deltacov = [];
delta_P = zeros(ntrial,2);
pResp = zeros(ntrial,2);
best = [];

% M is a 3 by 2 matrix:     C is a 3 by 3 matrix:
% 
%   A1-O1   A1-O2           A1.A1   A1.A2   A1.B
%    
%   A2-O1   A2-O2           A2.A1   A2.A2   A2.B
% 
%   B-O1    B-O2            B.A1    B.A2    B.B

% beliefs & deltas: 
% col1    col2    col3    col4    col5    col6    
% A1-O1   A2-O1   B-O1    A1-O2   A2-O2   B-O2  

for trial = 1:ntrial
    
    %% action selection
    maxM = max(M'); A1 = maxM(1); A2 = maxM(2); B = maxM(3);
    pA1 = exp(tau*A1)./(exp(tau*A1) + exp(tau*A2) + exp(tau*B));
    pA2 = exp(tau*A2)./(exp(tau*A1) + exp(tau*A2) + exp(tau*B));
    pB = exp(tau*B)./(exp(tau*A1) + exp(tau*A2) + exp(tau*B));
    best(trial,:) = [pA1,pA2,pB];
    
    action = rand();
     
    if action < pA1
        inputs(trial,:) = [1,0,1];  % action 1       
    elseif action > pA1 && action < (pA1+pA2)
        inputs(trial,:) = [0,1,1];  % action 2
    else
        inputs(trial,:) = [0,0,1];  % background (wait)
    end
        
    % calculate reward vector
    if inputs(trial,1) == 1
        rewards(trial,:) = [1,0,0].*contingency;
    elseif inputs(trial,2) == 1
        rewards(trial,:) = [0,1,0].*contingency;
    else
        rewards(trial,:) = [0,0,1].*contingency;
    end
    
    r=rand();
    rewards(trial,:)=r<rewards(trial,:);
     
    % create output state vector
    if any(rewards(trial,:))                    % If action is rewarded
        if rewards(trial,2) == 1
            outputs(trial,:) = [0,1];           % output = O2
        else
            outputs(trial,:) = [1,0];           % output = O1
        end
    else
        outputs(trial,:) = [0,0];               % output = no reward
    end
    
    %% Kalman filter
    i = inputs(trial,:)';           % inputs are a column vector
    o = outputs(trial,:);           % outputs are a row vector
    
    alpha = (v + i'*C*i);           % total uncertainty (outcome + expected)
    PE = C*i*(o' - M'*i)';          % uncertainty * summed prediction-error
    delta_mu = PE*1./alpha;         % calculate delta mu
    M = M + delta_mu;               % update mu
    
    cmatrix = C*(i*i')*C;           % covariance^2 by inputs
    delta_C = 1/alpha*cmatrix;      % adjust covariance by uncertainty
    C = C - delta_C;                % update C (reduce uncertainty)
    
    % store values
    beliefs(end+1,:) = reshape(M,[1,6]);
    covariance(end+1,:) = reshape(C,[1,9]);
    deltamu(end+1,:) = reshape(delta_mu,[1,6]);
    deltacov(end+1,:) = reshape(delta_C,[1,9]);
    
    % calculate delta P
    a = sum(inputs(:,1).*outputs(:,1))/sum(inputs(:,1)); % A = P(outcome|cue)
    absent = find(inputs(1:trial,1)==0 & inputs(1:trial,2)==0);
    c = sum(outputs(absent,1))/length(absent);          % C = P(outcome|~cue)
    delta_P(trial,1) = a - c;
    
    % calculate delta P
    a = sum(inputs(:,2).*outputs(:,2))/sum(inputs(:,2)); % A = P(outcome|cue)
    absent = find(inputs(1:trial,1)==0 & inputs(1:trial,2)==0);
    c = sum(outputs(absent,2))/length(absent);          % C = P(outcome|~cue)
    delta_P(trial,2) = a - c;
    
end

%% plot action values
figure(1)
subplot(3,1,1)
plot(beliefs(:,1)); hold on; plot(beliefs(:,5)); 
legend('A1 belief','A2 belief'); ylim([-.5 .5]); xlim([0 ntrial])
hold off; 

subplot(3,1,2) 
plot(delta_P); ylim([-1 1]); xlim([0 ntrial])
legend('dP(A1O1)','dP(A2O2)')

subplot(3,1,3)
span = 11;
x1 = smooth(inputs(:,1),span);
x2 = smooth(inputs(:,2),span);
diff = x2 - x1;
plot([x1,x2]); hold on
plot(diff,'k:'); hold off


% stem(deltamu(:,1));hold on
% stem(deltamu(:,3));
% plot(deltacov(:,3),':k'); hold off
% xlim([0 ntrial])
% legend('deltaAO','deltaXO','Cov')

% figure(2)
% subplot(2,1,1)
% stem(inputs(:,1),'Marker','none','LineWidth',2);
% ylim([-1 2]); xlim([0 ntrial])
% hold off
% 
% subplot(2,1,2)
% stem(outputs(:,1),'Color','r','Marker','none','LineWidth',2); 
% ylim([-1 2]); xlim([0 ntrial])
% hold off

a = sum(inputs(:,1:2),1)./ntrial;
a(3) = 1-(a(1)+a(2))

% %% function outputs
% model_values = transitions(2:end,1);
% model_values(:,2) = transitions(2:end,5);
% model_delta(:,1) = transitions(2:end,1)-delta_P(:,1);
% model_delta(:,2) = transitions(2:end,5)-delta_P(:,2);
% resp = inputs(:,1:2);

