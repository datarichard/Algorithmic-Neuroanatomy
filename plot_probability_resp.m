degfiles  = dir ('data/*deg_preprocessed*');
sigfiles = dir ('data/*sig_preprocessed*');

deg_means = get_means(degfiles);
sig_means = get_means(sigfiles);

X = categorical({'Unsignalled','Signalled'});
X = reordercats(X,{'Unsignalled','Signalled'}); 
vals = [deg_means; sig_means];
bar(X, vals); 
ylim([0, 1])
% title()
legend('Con', 'Deg')

function [means] = get_means(filenames)
degraded = zeros(length(filenames),1); 
conting = zeros(length(filenames),1);

for i = 1:size(filenames,1)
    fullfilename = fullfile('data', filenames(i).name);
    [Pdeg,Pcon] = probability_resp(fullfilename);
    degraded(i) = Pdeg;
    conting(i) = Pcon;
end

means = [mean(conting), mean(degraded)];

end


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
    denominator = 65;
else
    denominator = 50;
end

Pdeg = mean(sum(Fdeg) / denominator);
Pcon = mean(sum(Fcon) / denominator);

end