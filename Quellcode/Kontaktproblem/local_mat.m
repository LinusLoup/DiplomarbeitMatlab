function [S,f_local,J] = local_mat(points,lambda,mu,uS_local_x,uS_local_y,option)
%LOCAL_MAT computes the Jacobian J, a local stiffness matrix S and if 
%option='bubble' a local vector f_local, which will be needed for a(u_S,*) 
%in rho_S. LOCAL_MAT expects the nodes (points) of the local triangle, the 
%z-values uS_local of this nodes and an option, which can be 'linear' or 
%'bubble'.

% initialization of f_local:
f_local = zeros(6,1);

% x- and y-values of the nodes:
x = points(1,:);
y = points(2,:);

% the Jacobian J and other factors for the affine transformation from the 
% reference element onto a arbitrary triangle:
J = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));
a = 1/J * (mu*(x(3)-x(1))^2 + (2*mu+lambda)*(y(3)-y(1))^2);
b = -1/J * (mu*(x(3)-x(1))*(x(2)-x(1)) + (2*mu+lambda)*(y(3)-y(1))*(y(2)-y(1)));
c = -1/J * (mu*(y(3)-y(1))*(y(2)-y(1)) + (2*mu+lambda)*(x(3)-x(1))*(x(2)-x(1)));
d = 1/J * (mu*(x(3)-x(1))*(y(2)-y(1)) + lambda*(y(3)-y(1))*(x(2)-x(1)));
e = -1/J * (mu+lambda)*(x(3)-x(1))*(y(3)-y(1));
f = 1/J * (mu*(x(2)-x(1))^2 + (2*mu+lambda)*(y(2)-y(1))^2);
g = 1/J * (mu*(y(3)-y(1))^2 + (2*mu+lambda)*(x(3)-x(1))^2);
h = 1/J * (mu*(x(2)-x(1))*(y(3)-y(1)) + lambda*(x(3)-x(1))*(y(2)-y(1)));
i = -1/J * (mu+lambda)*(x(2)-x(1))*(y(2)-y(1));
j = 1/J * (mu*(y(2)-y(1))^2 + (2*mu+lambda)*(x(2)-x(1))^2);

switch lower(option)
    case {'linear'}
    % local stiffness matrix for linear shape functions:
    S = 1/2 * [ a+2*b+f, d+e+h+i, -a-b, -e-h, -b-f, -d-i;
                d+e+h+i, 2*c+g+j, -d-e, -c-g, -h-i, -c-j;
                -a-b,   -d-e,       a,   e,   b,    d;
                -e-h,   -c-g,       e,   g,   h,    c;
                -b-f,   -h-i,       b,   h,   f,    i;
                -d-i,   -c-j,       d,   c,   i,    j];
    
    case {'bubble'}
    % local stiffness matrix for quadratic shape functions:
    S = 4/3* diag([a+b+f, c+g+j, a+b+f, c+g+j, a+b+f, c+g+j]);
    
    % gradient of u_S:
    gradu_S_x = gradu(points,uS_local_x);
    gradu_S_y = gradu(points,uS_local_y);
    % transformation onto the reference element:
    gradu_S_x = [x(2)-x(1), y(2)-y(1); x(3)-x(1), y(3)-y(1)]*gradu_S_x';
    gradu_S_y = [x(2)-x(1), y(2)-y(1); x(3)-x(1), y(3)-y(1)]*gradu_S_y';
    gradu_S = [gradu_S_x; gradu_S_y];
    % f_local:
    f_local = -(([a, b, e, d]*gradu_S).*[0; 0; 2/3; 0; -2/3; 0] + ...
        ([b, f, h, i]*gradu_S).*[-2/3; 0; 2/3; 0; 0; 0]...
        +([e, h, g, c]*gradu_S).*[0; 0; 0; 2/3; 0; -2/3] + ...
        ([d, i, c, j]*gradu_S).*[0; -2/3; 0; 2/3; 0; 0]);
        
    otherwise
        error('Unknown option');
end

end

