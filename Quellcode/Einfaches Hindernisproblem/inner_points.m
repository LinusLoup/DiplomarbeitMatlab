function non_bound_index = inner_points(adjacency)
%DIRICHLET_BOUNDARY assembles the inner points, that is the non boundary
%points, by using the adjacency matrix given by assemb.

% evaluating the inner points out of the adjacency matrix:
[m,n] = size(adjacency);
non_bound_index = find(prod(adjacency == zeros(m,n)));

end