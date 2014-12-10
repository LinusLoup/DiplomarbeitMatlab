function [S,f_local,J] = local_mat(points,uS_local,option)
%LOCAL_MAT computes the Jacobian J, a local stiffness matrix S and if 
%option='bubble' a local vector f_local, which will be needed for a(u_S,*) 
%in rho_S. LOCAL_MAT expects the nodes (points) of the local triangle, the 
%z-values uS_local of this nodes and an option, which can be 'linear' or 
%'bubble'.

% initialization of f_local:
f_local = zeros(3,1);

% x- and y-values of the nodes:
x = points(1,:);
y = points(2,:);

% the Jacobian J and other factors for the affine transformation from the
% reference element onto a arbitrary triangle:
J = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));
a = 1/J * ((x(3)-x(1))^2 + (y(3)-y(1))^2);
b = -1/J * ((y(2)-y(1))*(y(3)-y(1)) + (x(2)-x(1))*(x(3)-x(1)));
c = 1/J * ((x(2)-x(1))^2 + (y(2)-y(1))^2);

switch lower(option)
    case {'linear'}
    % local stiffness matrix for linear shape functions:
    S = [ a/2+b+c/2, -a/2-b/2, -b/2-c/2;
        -a/2-b/2,    a/2,      b/2;
        -b/2-c/2,    b/2,      c/2     ];   
    case {'bubble'}
    % local stiffness matrix for quadratic shape functions:
    S_quad = 4/3* [ a+b+c, -b-c, b;
                    -b-c,  a+b+c, -a-b;
                      b,    -a-b, a+b+c ];
    S = diag(diag(S_quad));
    
    % gradient of u_S:
    gradu_S = gradu(points,uS_local);
    % transformation onto the reference element:
    gradu_S = [x(2)-x(1), y(2)-y(1); x(3)-x(1), y(3)-y(1)]*gradu_S';
    % f_local:
    f_local = -(a*gradu_S(1)+b*gradu_S(2))*[0; 2/3; -2/3]...
        -(b*gradu_S(1)+c*gradu_S(2))*[-2/3; 2/3; 0];
    case {'quadratic'}
    % local stiffness matrix for quadratic shape functions:
    S = 4/3* [ a+b+c, -b-c, b;
               -b-c,  a+b+c, -a-b;
                 b,    -a-b, a+b+c ];        
    otherwise
        error('Unknown option');
end

end

