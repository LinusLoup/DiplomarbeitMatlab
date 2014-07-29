function [u,u_derivative,lambda,energy_exact,u_H1seminorm_squared,R]=u_exact_ring_spherical_obstacle(xy,f,phi,Lx,Ly,rho)

x=xy(:,1);
y=xy(:,2);
clear xy

x=x-1;
y=y-1;

distance=sqrt(x.^2+y.^2);

if (f/phi >4)
    
    %podminka pro mezni psi
    %f*rho^2*(sin(psi))^2*(1-2*log(rho*sin(psi)))-4*rho*(1-cos(psi))+(4*rho*(sin(psi))^2*log(rho*sin(psi))/cos(psi))-f+4*phi=0 
    fun=@(psi)(f*rho^2*(sin(psi)).^2).*(1-2*log(rho*sin(psi)))-4*rho*(1-cos(psi))+((4*rho*(sin(psi)).^2).*log(rho*sin(psi))/cos(psi))-f+4*phi;
    psi=fsolve(fun,atan(1/rho)/2,optimset('Display','off'));      
    
    A=(f*rho^2*(sin(psi))^2-f+4*phi-4*rho*(1-cos(psi)))/(4*log(rho*sin(psi)));
  
    R=rho*sin(psi);
    index_in_contact=sort(find((distance<R)));
    index_off_contact=sort(find((distance>=R)));

    u=zeros(size(distance));
    u(index_off_contact)=A*log(distance(index_off_contact))+(1/4)*f*(1-distance(index_off_contact).^2);
    u(index_in_contact)= phi-rho+sqrt(rho^2-distance(index_in_contact).^2);

     u_x=zeros(size(distance));
     u_y=zeros(size(distance));
     u_x(index_off_contact)=(A*x(index_off_contact))./distance(index_off_contact).^2-f*x(index_off_contact)/2;
     u_y(index_off_contact)=(A*y(index_off_contact))./distance(index_off_contact).^2-f*y(index_off_contact)/2;

     Bvector=sqrt(rho^2-distance(index_in_contact).^2);

     u_x(index_in_contact)=-x(index_in_contact)./Bvector;
     u_y(index_in_contact)=-y(index_in_contact)./Bvector;

     lambda=zeros(size(x));
     lambda(index_in_contact)=-f+(2*Bvector+(distance(index_in_contact).^2)./Bvector)./(Bvector.^2);
    
     psi_max=atan(1/rho);
     
     psi_save=psi;
     psi_testing=(psi_max/1000):(psi_max/1000):psi_max;
     %psi_testing=0.15:(psi_max/1000):0.25;
     psi_testing=[psi_testing psi_save];
     energy_testing=zeros(size(psi_testing));
     
     for i=1:numel(psi_testing)
            psi=psi_testing(i);    
            %seminorma
            I1(i)=-2*pi*A^2*log(rho*sin(psi))-A*f*pi*(1-rho^2*(sin(psi))^2)+(1/8)*pi*f^2*(1-rho^4*(sin(psi))^4);
            I2(i)=-pi*rho^2*((sin(psi))^2+log((cos(psi))^2));
            u_H1seminorm_squared(i)=I1(i)+I2(i);
    
            %energie
            I3(i)=(pi/2)*A*f*(rho^2*(sin(psi))^2-2*rho^2*(sin(psi))^2*log(rho*sin(psi))-1)+...
               (1/4)*pi*f^2*(1-rho^2*(sin(psi))^2)-(1/8)*pi*f^2*(1-rho^4*(sin(psi))^4);
            I4(i)=pi*f*rho^2*(phi-rho)*(sin(psi))^2+(1/3)*2*pi*f*rho^3*(1-sqrt((1-(sin(psi))^2)^3));
                       
            %energy_testing(i)=u_H1seminorm_squared;
            %energy_testing(i)=-(I3+I4);
     end
     
     energy_testing=(1/2)*u_H1seminorm_squared-(I3+I4);
     
     psi=psi(end);
     energy_exact=energy_testing(end);
     u_H1seminorm_squared=u_H1seminorm_squared(end); 
     
     
     psi_testing(end)=[];
     energy_testing(end)=[];
     I1(end)=[];
     I2(end)=[];
     I3(end)=[];
     I4(end)=[];
     
     energy_minimum_numerical=min(energy_testing);
     
%      figure(1000)
%      hold on
%      plot(psi_testing,I1/2+I2/2-I3-I4);
%      plot(psi_testing,fun(psi_testing));
%      plot(psi,fun(psi),'ro');
%      plot(psi,energy_exact,'ro');
%      grid on
%      hold off
     

     
     
    %energy testing
    
    
 
%     
%     energy_exact=(-pi*f^2/16)*(1-R^4+R^4/log(R))-pi*f*R^2*(phi-f/4)/2/log(R)-pi*(phi-f/4)^2/log(R);
%     
%     A=(4*phi+f*R^2-f)/4/log(R);
%     u_H1seminorm_squared=(pi*f^2/8)*(1-R^4)-A*pi*f*(1-R^2)-2*A^2*pi*log(R);
    
else
    R=0;
    
    u=(f/4)*(-distance.^2+1);  
    
    u_x=-f*(x-(Lx/2))/2;
    u_y=-f*(y-(Ly/2))/2;
    
    lambda=zeros(size(x));
    
    energy_exact=-pi*f^2/16;
    
    u_H1seminorm_squared=pi*(f^2)/8;
    
end

u_derivative=[u_x u_y];

 


end
