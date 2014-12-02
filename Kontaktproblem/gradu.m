function grad_u = gradu(points,zvalues)
%GRADU determines for given (x,y,z)-values the gradient of the function u=a*x+b*y+c*z, that is grad u=(a,b).

% x- and y-values of the points:
x = points(1,:);
y = points(2,:);

% the Jacobian:
J = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));

% evaluating grad_u:
grad_u = 1/J * ([y(3)-y(1), y(1)-y(2); x(1)-x(3), x(2)-x(1)]*...
    [zvalues(2)-zvalues(1);zvalues(3)-zvalues(1)])';

end