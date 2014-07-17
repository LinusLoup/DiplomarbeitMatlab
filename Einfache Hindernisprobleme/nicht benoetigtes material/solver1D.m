n = 100;

L = 1;

x = linspace(-L,L,n)';
h = L/(n-1);

a = 0.3;
b = 1;
g = a-b*x.^2;
%g = -(3*x.^2-1).^2+1;
u = zeros(n,1);
u_quad = u;
u_jacobi = u;
f = -0.1*zeros(n, 1);

A = zeros(n,n);

for i=1:n
    for j=1:n
        if abs(j-i) == 1
            A(i,j) = -1/h;
        elseif j==i
            A(i,j) = 2/h;
        else
            A(i,j) = 0;
        end
    end
end

H = zeros(2,n);
H(1,1) = 1;
H(2,n) = 1;
u_quad = quadprog(A,-f,[],[],H,zeros(2,1),g,[]);

A_new = A(2:n-1,2:n-1);
f_new = f(2:n-1);
uu = projected_jacobi(A_new,f_new,g(2:n-1),ones(n-2,1),1e-18);
u_jacobi(2:n-1) = uu;

% Analytische Lösung:
xhat = 1-sqrt(1-a/b);
alpha = (a-b*xhat^2)/(1-xhat);
u_ana = zeros(n,1);
for i = 1:n
   if x(i) >= -xhat && x(i) <= xhat
       u_ana(i) = g(i);
   elseif x(i) < -xhat
       u_ana(i) = alpha*(1+x(i));
   else
       u_ana(i) = alpha*(1-x(i));
   end
end


subplot(3,1,1);
plot(x,g,x,u_quad);axis([-1,1,-0.2,1.2]);
subplot(3,1,2); 
plot(x,g-u_ana);axis([-1,1,-1,1]);
subplot(3,1,3);
plot(x,g,x,u_ana);axis([-1,1,-0.2,1.2]);