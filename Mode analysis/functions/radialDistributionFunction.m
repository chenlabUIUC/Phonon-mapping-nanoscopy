function rdf = radialDistributionFunction(points, binWidth, maxDistance, plotting)
    if nargin<4
        plotting = 0;
    end

    % Points: n*2 array of point coordinates
    % binWidth: width of the bins for the histogram
    % maxDistance: maximum distance to consider for the RDF
    
    % Number of points
    n = size(points, 1);
    
    % Calculate pairwise distances
    distances = pdist(points);
    
    % Create bins for the histogram
    edges = 0:binWidth:maxDistance;
    binCenters = edges(1:end-1) + binWidth / 2;
    
    % Histogram of pairwise distances
    histCounts = histcounts(distances, edges);
    
    % Normalize the histogram to obtain RDF
    area = pi * ((edges(2:end)).^2 - (edges(1:end-1)).^2);
    idealGasDensity = n / (pi * max(max(points(:,1)) - min(points(:,1)), max(points(:,2)) - min(points(:,2)))^2);
    rdf = histCounts ./ (area * idealGasDensity * n);
    
    % Plot the RDF
    % figure;
    if plotting
        plot(binCenters, rdf, 'LineWidth', 2);
        xlabel('Distance');
        ylabel('g(r)');
        title('Radial Distribution Function');
    end
end