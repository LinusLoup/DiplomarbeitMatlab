function x = projected_jacobi(A,b,g,x0,tol)
%PROJECTED_JACOBI evaluates the solution of a linear varational inequality
%(Ax*-b)^T(x-x*)>=0 s.t. x>=g by using the projected Jacobi method.
%
%x0 is the start-solution and tol the error between an iteration x_i and
%the (projected) next one x_{i+1}.


%% Initialisierung
n = length(b);

if nargin == 3
    x0 = zeros(n,1);
    tol = 1e-10;
end

if nargin == 4
    tol = 1e-10;
end

xn = x0;
xn_new = zeros(n,1);


%% Schleife des Jacobiverfahrens
while 1
    for i = 1:n
        xn_new(i) = (b(i)-A(i,:)*xn+A(i,i)*xn(i))/A(i,i);
        
        if xn_new(i) < g(i)
            xn_new(i) = g(i);
        end
    end
    
    err = norm(xn-xn_new);
    
    if err < tol
        break;
    end
    
    xn = xn_new;
end

x = xn_new;

end

