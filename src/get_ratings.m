
function [means] = get_ratings(filenames)

means = zeros(length(filenames), 2);

for i = 1:length(filenames)
    fullfilename = fullfile('data', filenames(i).name);
    load(fullfilename, 'ratings');
 
    rt = zeros(6,2);  
    % non-degraded        % degraded
    rt(1,1) = ratings(1,2); rt(1,2) = ratings(1,1);
    rt(2,1) = ratings(2,1); rt(2,2) = ratings(2,2);
    rt(3,1) = ratings(3,1); rt(3,2) = ratings(3,2);
    rt(4,1) = ratings(4,2); rt(4,2) = ratings(4,1);
    rt(5,1) = ratings(5,2); rt(5,2) = ratings(5,1);
    rt(6,1) = ratings(6,1); rt(6,2) = ratings(6,2);
    means(i, :) = mean(rt);
    
end

end