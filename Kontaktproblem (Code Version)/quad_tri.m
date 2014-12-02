function [wi,gauss,ansatz_values] = quad_tri(points,ansatz_fun,num_of_nodes)
%QUAD_TRI computes the weights, Gauss-points and functionvalues of the shape functions, which are needed for the quadrature. The matrix points stores the x-,y-values of the nodes of the triangle, ansatz_fun a vector of shape functions and num_of_nodes the number of nodes for the quadrature formula.

% x- and y-values of the points:
x = points(1,:);
y = points(2,:);

switch num_of_nodes
        % determination of the weights and the local points xi,eta:
        case 3
            wi = [1/6,1/6,1/6];
            local_poi = [1/2,1/2,0;0,1/2,1/2];
            xi = local_poi(1,:);
            eta = local_poi(2,:);
            
        case 7
            wh = [(155-sqrt(15))/2400,(155+sqrt(15))/2400];
            
            wi = [9/80,wh(1),wh(1),wh(1),wh(2),wh(2),wh(2)];
            local_poi = [1/3,(6-sqrt(15))/21,(9+2*sqrt(15))/21,...
                (6-sqrt(15))/21,(6+sqrt(15))/21,(9-2*sqrt(15))/21,...
                (6+sqrt(15))/21;1/3,(6-sqrt(15))/21,(6-sqrt(15))/21,...
                (9+2*sqrt(15))/21,(6+sqrt(15))/21,(6+sqrt(15))/21,...
                (9-2*sqrt(15))/21];
            xi = local_poi(1,:);
            eta = local_poi(2,:);
            
        otherwise
            error('keine Quadraturformel bekannt');
            
end  

% evaluation of the Gauss-points by using xi and eta with the points x and y:
gauss(1,:) = x(1)+(x(2)-x(1))*xi+(x(3)-x(1))*eta;
gauss(2,:) = y(1)+(y(2)-y(1))*xi+(y(3)-y(1))*eta;

% computing the functionvalues of the shape functions in the local coordinates: 
ansatz_values = ansatz_fun(xi,eta);

end