function Nplus_set = Nplus(contact_set,inner_points,nodes)
%NPLUS computes the set N^+ of the inner non-contact nodes, given the set
%of contact nodes N^0, inner points and the set of all nodes.

np = length(nodes);
index = setdiff(1:np,contact_set);
Nplus_set = intersect(index,inner_points);

end