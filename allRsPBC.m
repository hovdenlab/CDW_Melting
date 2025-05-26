function rs = allRsPBC(xs, ys, r_c, Lx, Ly)
    % allRsPBC computes the pairwise distances between points (xs, ys)
    % under periodic boundary conditions and returns those distances
    % which are less than the cutoff distance r_c.
    %
    % Inputs:
    %   xs - A vector of x coordinates of the points
    %   ys - A vector of y coordinates of the points
    %   r_c - The cutoff distance
    %   Lx - The width of the periodic domain
    %   Ly - The height of the periodic domain
    %
    % Output:
    %   rs - A vector of distances less than r_c

    % Number of points
    numPoints = numel(xs);

    % Find pairwise differences in x and y using meshgrid (accounting for PBC)
    [x1, x2] = meshgrid(xs, xs);
    [y1, y2] = meshgrid(ys, ys);
    
    % Calculate pairwise distance components considering PBC using minimum image convention
    dx = min(abs(x1 - x2), Lx - abs(x1 - x2));
    dy = min(abs(y1 - y2), Ly - abs(y1 - y2));

    % Calculate pairwise distances
    distances = sqrt(dx.^2 + dy.^2);

    % Exclude distances greater than the cutoff and distances of each point to itself
    distances = distances(triu(true(numPoints), 1) & distances < r_c);

    % Reshape the matrix into a vector of distances less than r_c
    rs = distances(:);
end