function [neighbours,flag] = neighbourhood(edges_or_point,triangles,option)
%EDGE_NEIGHBOURHOOD computes the indices of the triangles, which are the neighbours for the given edges or point. The vector flag shows, if an edge has got only one neighbour (flag = 0) or two (flag = 1).

n = length(edges_or_point);
flag = zeros(1,n);

switch lower(option)
    case {'edge','edges'}
        neighbours = zeros(2,n);

        for i = 1:n
            [~,neighbours(:,i)] = find(triangles==edges_or_point(i));
    
            if neighbours(1,i) ~= neighbours(2,i)
                flag(i) = 1;
            end
        end
        
    case {'point'}
        [~,neighbours] = find(triangles(1:3,:)==edges_or_point);
end