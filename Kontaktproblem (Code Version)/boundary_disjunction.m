function [gamma_D,gamma_N,gamma_C] = boundary_disjunction(boundary_indices,points)
%BOUNDARY_DISJUNCTION Summary of this function goes here
%   Detailed explanation goes here

bound_points = points(:,boundary_indices);

gamma_C = boundary_indices(find(bound_points(2,:)>=0 & bound_points(2,:)<=1));
gamma_D = boundary_indices(find(bound_points(2,:)<0));
gamma_N = boundary_indices(find(bound_points(2,:)>1));

end

