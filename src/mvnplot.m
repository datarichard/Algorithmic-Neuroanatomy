function E = mvnplot(Mu,Sigma,h)
    % Displays bivariate normal plot and returns expected value
    % Written by R.W.Morris 7th April 2016
    % e.g.,
    % h = figure(1);
    % Mu = [0 0.5];
    % Sigma = [.1 .05; .05 .05];
    % E = mvnplot(Mu,Sigma,h);
    
    x1 = -1:.05:1; x2 = -1:.05:1;
    [X1,X2] = meshgrid(x1,x2);
    F = mvnpdf([X1(:) X2(:)],Mu,Sigma);
    F = reshape(F,length(x2),length(x1));
    if ~exist('h','var')
        h = 1;
    end
    figure(h)
    surf(x1,x2,F); view(0,90)
    caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
    axis([-1 1 -1 1 0 inf])
    xlabel('Ax'); ylabel('Cx'); zlabel('Probability Density');
    E = max(max(F));
end