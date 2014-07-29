function [u,u_derivative,lambda,energy_exact,u_H1seminorm_squared,R]=u_exact_ring_constant_obstacle(xy,f,phi,Lx,Ly)

x=xy(:,1);
y=xy(:,2);
distance=sqrt(sum((xy-[ones(size(xy,1),1)*Lx/2 ones(size(xy,1),1)*Ly/2]).^2,2));

if (f/phi >=4)
    R=fsolve(@(x)(x.^2).*(-2*log(x)+1)-1+4*phi/f,0.5,optimset('Display','off'));
    index=sort(find((distance<R)));
    
    u=(f/4)*(-distance.^2+1)+(4*phi+f*R^2-f)/4/log(R)*log(distance);  
    u(index)=phi*ones(size(index));
    
    u_x=(4*phi+f*R^2-f)*(x-(Lx/2))/4/log(R)./distance.^2-f*(x-(Lx/2))/2;
    u_y=(4*phi+f*R^2-f)*(y-(Ly/2))/4/log(R)./distance.^2-f*(y-(Ly/2))/2;
    u_x(index)=zeros(size(index));
    u_y(index)=zeros(size(index));
    
    lambda=zeros(size(x));
    lambda(index)=-f*ones(size(index));
    
    energy_exact=(-pi*f^2/16)*(1-R^4+R^4/log(R))-pi*f*R^2*(phi-f/4)/2/log(R)-pi*(phi-f/4)^2/log(R);
    
    A=(4*phi+f*R^2-f)/4/log(R);
    u_H1seminorm_squared=(pi*f^2/8)*(1-R^4)-A*pi*f*(1-R^2)-2*A^2*pi*log(R);
    
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

 

% if (r_contact<=0)
%     energy_exact=-f^2/24;    
%     
%     u=-f*x.^2/2+f*x/2; 
%     u_derivative=(1/2)*f*(1-2*x);
% 
%     
% else
%     energy_exact=f*phi*((4/3)*sqrt(2*phi/f)-1);
%     
%     index_left=sort(find((x-1/2)<-r_contact));
%     index_right=numel(x)+1-numel(index_left):numel(x);
%     index_middle= index_left(end)+1:index_right(1)-1;
% 
%     u_left=-(f/2)*x.^2+((phi(1)+(f/2)*(1/2-r_contact)^2)/(1/2-r_contact))*x;
% 
%     u(index_left)=u_left(index_left);
%     u(index_middle)=phi;
%     u(index_right)=u_left(sort(index_left,'descend'));
% 
%     u_derivative_left=-f*x-sqrt(2*phi(1)*f);
% 
%     u_derivative(index_left)=u_derivative_left(index_left);
%     u_derivative(index_middle)=0;
%     u_derivative(index_right)=-u_derivative_left(sort(index_left,'descend'));       
% end



end
