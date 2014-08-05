function [neighbours,flag] = neighbourhood(edges_or_point,triangles,option)
%EDGE_NEIGHBOURHOOD berechnet die Dreiecksindizes, die zu mitgegebenen
%Kanten edges, f�r ein Set aus Dreicken triangles, die Nachbarn sind.
%
%Der Vektor flag besteht hierbei aus Einsen und Nullen, um anzugeben,
%welche Kanten nur einen Nachbar haben.
%0: Es gibt nur einen Nachbarn; 1: Es gibt zwei Nachbarn.

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