function u_out = dirichlet_boundary(u_in,adjacency,bounds)
%DIRICHLET_BOUNDARY assembles the boundary-values of the Dirichlet-Boundary
%   Value Problem. The vector u is the solution of the non-Boundary points
%   and adjacency is a Matrix, in which the adjacency of the Boundary
%   points is stored. The vector bounds contents the boundary-values of
%   each boundary-point.

[m,n] = size(adjacency);

%% Berechnung der Adjazenzmatrix f�r die inneren Punkte:
non_bound_ind = find(prod(adjacency == zeros(m,n)));
l = length(non_bound_ind);
non_bound_adj = zeros(l,n);
non_bound_adj(:,non_bound_ind) = eye(l,l);

%% Berechnung der L�sung durch "Aufbl�hen" u_in und bounds Vektoren:
b = adjacency' * bounds;
uu = non_bound_adj' * u_in;

u_out = uu + b;

end