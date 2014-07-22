function [neighbours,flag] = edge_neighbourhood(edges,triangles)
%EDGE_NEIGHBOURHOOD berechnet die Dreiecksindizes, die zu mitgegebenen
%Kanten edges, für ein Set aus Dreicken triangles, die Nachbarn sind.
%
%Der Vektor flag besteht hierbei aus Einsen und Nullen, um anzugeben,
%welche Kanten nur einen Nachbar haben.
%0: Es gibt nur einen Nachbarn; 1: Es gibt zwei Nachbarn.

n = length(edges);
neighbours = zeros(2,n);
flag = zeros(1,n);
my_triangle = [triangles;zeros(1,size(triangles,2))];

for i = 1:n
    help_index = find(my_triangle==edges(i));
    neighbours(:,i) = ceil(help_index/4);
    
    if neighbours(1,i) ~= neighbours(2,i)
        flag(i) = 1;
    end
end

end

