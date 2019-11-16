
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