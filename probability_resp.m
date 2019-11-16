
function [Pdeg,Pcon] = probability_resp(filename)

load(filename,'trialStarts','LeftPressesBlock','RightPressesBlock');

Fdeg = zeros(150,6);
Fcon = zeros(150,6);

for b = 1:6
    if b ==1||b==4||b==5
        degPresses = LeftPressesBlock{b};
        conPresses = RightPressesBlock{b};
    else
        degPresses = RightPressesBlock{b};
        conPresses = LeftPressesBlock{b};        
    end
    
    degPresses = degPresses - trialStarts(b);
    degPresses = unique(fix(degPresses));
    degPresses(degPresses == 0) = 1;
    Fdeg(degPresses, b) = 1;

    conPresses = conPresses - trialStarts(b);
    conPresses = unique(fix(conPresses));
    conPresses(conPresses == 0) = 1;
    Fcon(conPresses, b) = 1; 
    
end

Fdeg = Fdeg(1:120,:);
Fcon = Fcon(1:120,:);

if contains(filename, 'deg')
    denominator = 65; % Average time per block in this condition
else
    denominator = 50; % Average time per block in this condition
end

Pdeg = mean(sum(Fdeg) / denominator);
Pcon = mean(sum(Fcon) / denominator);

end