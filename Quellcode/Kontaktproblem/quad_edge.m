function [wi,gauss,J,ansatz_values] = quad_edge(points,ansatz_fun,num_of_nodes)
%QUAD_EDGE computes the weights, Gauss-points and functionvalues of the 
%shape functions, which are needed for the quadrature over an edge of a 
%triangle. The matrix points stores the x-,y-values of the nodes of the 
%edge, ansatz_fun a vector of shape functions and num_of_nodes the number 
%of nodes for the quadrature formula.

% x- and y-values of the points:
x = points(1,:);
y = points(2,:);

switch num_of_nodes
    % determination of the weights and the local points xi,eta:
    case 1
        wi = 2;
        local_poi = 0;
            
    case 2
        wi = [1,1];
        local_poi = [-sqrt(1/3),sqrt(1/3)];
        
    case 3 
        wi = [5/9,8/9,5/9];
        local_poi = [-sqrt(3/5),0,sqrt(3/5)];
        
    case 4
        wi = [18-sqrt(30),18+sqrt(30),18+sqrt(30),18-sqrt(30)]/36;
        local_poi = [-sqrt(3/7+2/7*sqrt(6/5)),-sqrt(3/7-2/7*sqrt(6/5)),...
            sqrt(3/7-2/7*sqrt(6/5)),sqrt(3/7+2/7*sqrt(6/5))];
            
    case 5
        wi = [(322-13*sqrt(70))/900,(322+13*sqrt(70))/900,128/225,...
            (322+13*sqrt(70))/900,(322-13*sqrt(70))/900];
        local_poi = [-1/3*sqrt(5+2*sqrt(10/7)),-1/3*sqrt(5-2*sqrt(10/7)),...
            0,1/3*sqrt(5-2*sqrt(10/7)),1/3*sqrt(5+2*sqrt(10/7))];
        
    otherwise
        error('keine Quadraturformel bekannt');
            
end  

% evaluation of the Gauss-points by using xi and eta with the points x and y:
l = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
parameter = l/2.*local_poi+l/2;
gauss(1,:) = x(1) + parameter.*(x(2)-x(1));
gauss(2,:) = y(1) + parameter.*(y(2)-y(1));

% evaluating the jacobian:
J = l/2;

% computing the functionvalues of the shape functions in the local coordinates: 
ansatz_values = ansatz_fun(local_poi);

end