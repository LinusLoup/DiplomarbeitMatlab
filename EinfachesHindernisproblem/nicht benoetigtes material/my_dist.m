function d = my_dist(points,dist_points)
%MY_DIST evaluates a vector of distances. Each component k is the distance
%between the k-th column-vector of the matrix "points" to the points of the
%matrix "dist_points".
%
%If points is a (n,1)-matrix, my_dist will evaluate just the distance
%between the vector "points" and the given set of points "dist_points".
%
    n = size(dist_points,2);
    j = size(points,2);
    all_length = zeros(1,n);
    d = zeros(1,j);
    
    for l = 1 : j
        for k = 1 : n
            all_length(k) = norm(points(:,l)-dist_points(:,k));
        end
        d(l) = min(all_length);
        all_length = zeros(1,n);
    end

end

