function [u,u_derivative,lambda,energy_exact,u_H1seminorm_squared,R,f,phi]=u_exact_square_constant_obstacle(xy)

Lx=2;
Ly=2;

x=xy(:,1)-(Lx/2);
y=xy(:,2)-(Ly/2);

distance=sqrt(x.^2+y.^2);

R=0.7;
index=sort(find((distance<R)));

u=(distance.^2-R^2).^2;
u(index)=0*ones(size(index));

f=-8*(2*distance.^2-R^2);
f(index)=-8*R^2*(1-(distance(index).^2-R^2));

u_x=4*x.*(distance.^2-R^2);
u_y=4*y.*(distance.^2-R^2);
u_x(index)=zeros(size(index));
u_y(index)=zeros(size(index));

%f(index)=0*f(index);
%f=0*f;


lambda=zeros(size(x));
lambda(index)=-f(index);

I1=2/7+2/5-1.96/5-1.96/9+(2*0.49^2/3);
I2=2/5+2/9-1.96/3+0.49^2;
I3=2*pi*(0.7^8/8-0.98*0.7^6/6+0.49^2*0.7^4/4);
I4=2*pi*(0.7^6/6-0.98*0.7^4/4+0.49^2*0.7^2/2);

energy_exact=96*I1-24*I3-15.84*I2+3.92*I4;  
u_H1seminorm_squared=64*I1-16*I3;      

energy_exact=19.49622366;
u_H1seminorm_squared=15.26106850;
    
u_derivative=[u_x u_y];

phi=0*u;

P=(24/35)-((56/45)*R^2)+((2/3)*R^4);
D=((28/45)-((4/3)*R^2)+(R^4))*(R^2);
T=(2/3)*pi*(R^8);

energy_exact=(96*P)-(32*D)+T;
u_H1seminorm_squared=(64*P)-(2*T);

end
