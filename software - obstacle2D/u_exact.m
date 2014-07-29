function [u,u_gradient,u_laplacian,mu,distance]=u_exact(x,phi,Lx,Ly,radius,mu_aplitude)
%edge 1 - down, 2 - right, 3 - top, 4 - left
%f returns -u_laplacian + mu, where mu only adds to inside the prescribed radius 

epsilonX=0.3*Lx;
epsilonY=0.3*Ly;
amplitude_radial=0.5;
amplitude_boundary=1;

distance=sqrt(sum((x-[ones(size(x,1),1)*Lx/2 ones(size(x,1),1)*Ly/2]).^2,2));

u1=u_boundary(x(:,2),epsilonY);
u2=u_boundary(Lx-x(:,1),epsilonX);
u3=u_boundary(Ly-x(:,2),epsilonY);
u4=u_boundary(x(:,1),epsilonX);
uRadial=u_radial_function(distance,radius,phi,amplitude_radial);

u1_gradient=[zeros(size(x,1),1) u_derivative_boundary(x(:,2),epsilonY)];
u2_gradient=[-u_derivative_boundary(Lx-x(:,1),epsilonX) zeros(size(x,1),1)];
u3_gradient=[zeros(size(x,1),1) -u_derivative_boundary(Ly-x(:,2),epsilonY)];
u4_gradient=[u_derivative_boundary(x(:,1),epsilonX) zeros(size(x,1),1)]; 
uRadial_gradient=u_radial_function_gradient(distance,x,radius,Lx,Ly,amplitude_radial);

u1_laplacian=u_secondderivative_boundary(x(:,2),epsilonY);
u2_laplacian=u_secondderivative_boundary(Lx-x(:,1),epsilonX);
u3_laplacian=u_secondderivative_boundary(Ly-x(:,2),epsilonY);
u4_laplacian=u_secondderivative_boundary(x(:,1),epsilonX);
uRadial_laplacian=u_radial_function_laplacian(distance,radius,amplitude_radial); 

%uRadial=ones(size(uRadial));
%uRadial_gradient=zeros(size(uRadial_gradient));

u=((((u1.*u2).*u3).*u4).*uRadial);
%u=((((u1.*u2).*u3).*u4));
%u=uRadial;


u_gradient= u1_gradient.*kron([1 1],((u2.*u3).*u4).*uRadial) ... 
           +u2_gradient.*kron([1 1],((u1.*u3).*u4).*uRadial) ...
           +u3_gradient.*kron([1 1],((u1.*u2).*u4).*uRadial) ...
           +u4_gradient.*kron([1 1],((u1.*u2).*u3).*uRadial) ...
           +uRadial_gradient.*kron([1 1],((u1.*u2).*u3).*u4);      
%u_gradient=uRadial_gradient;   

u_laplacian= 2*(sum(u1_gradient.*uRadial_gradient,2)).*((u2.*u3).*u4) ...
           +2*(sum(u2_gradient.*uRadial_gradient,2)).*((u1.*u3).*u4) ...
           +2*(sum(u3_gradient.*uRadial_gradient,2)).*((u1.*u2).*u4) ...
           +2*(sum(u4_gradient.*uRadial_gradient,2)).*((u1.*u2).*u3) ...
           +u1_laplacian.*(((u2.*u3).*u4).*uRadial) ... 
           +u2_laplacian.*(((u1.*u3).*u4).*uRadial) ...
           +u3_laplacian.*(((u1.*u2).*u4).*uRadial) ...
           +u4_laplacian.*(((u1.*u2).*u3).*uRadial) ...
           +uRadial_laplacian.*(((u1.*u2).*u3).*u4);
       
mu=mu_inside_radius(distance,radius,mu_aplitude);
           
%u_laplacian=uRadial_laplacian



%*smoother(x(:,2),epsilonX);

    function value=u_boundary(d,epsilon)
    %d means the distance from the boundary
        value=amplitude_boundary*ones(size(d));
        index=find(d<=epsilon);
        value2=(-1/epsilon^2)*d.^2+(2/epsilon)*d;
        value(index)=amplitude_boundary*value2(index);
    end

    function value=u_derivative_boundary(d,epsilon)
    %d means the distance from the boundary
        value=zeros(size(d));
        index=find(d<=epsilon);
        value2=-2*(d-epsilon)/(epsilon*epsilon);
        value(index)=amplitude_boundary*value2(index);
    end

    function value=u_secondderivative_boundary(d,epsilon)
    %d means the distance from the boundary  
        value=zeros(size(d));
        index=find(d<=epsilon);
        value2=-2*ones(size(d))/(epsilon*epsilon);
        value(index)=amplitude_boundary*value2(index);
    end   
        
    function value=u_radial_function(distance,radius,phi,amplitude_radial)
    %distance means the distance from the center Lx/2 Ly/2
    
        value=phi*ones(size(distance));
        index=find(distance>=radius);
        value2=amplitude_radial*(distance-radius).^2+phi;
        value(index)=value2(index);
    end

    function value=u_radial_function_gradient(distance,x,radius,Lx,Ly,amplitude_radial)
    %distance means the distance from the center Lx/2 Ly/2
    
        value=zeros(size(distance),2);
        index=find(distance>=radius);
        value2=amplitude_radial*[(2*(distance-radius).*(x(:,1)-Lx/2))./distance ...
                (2*(distance-radius).*(x(:,2)-Ly/2))./distance];
        value(index,:)=value2(index,:);
    end

    function value=u_radial_function_laplacian(distance,radius,amplitude_radial)
    %distance means the distance from the center Lx/2 Ly/2
    
        value=zeros(size(distance));
        index=find(distance>=radius);
        value2=amplitude_radial*4*(distance-radius)./distance;
        value(index)=value2(index);
    end
    
    function value=mu_inside_radius(distance,radius,mu)
    %distance means the distance from the center Lx/2 Ly/2
    
        value=zeros(size(distance));
        index=find(distance<radius);
        value(index)=ones(numel(index),1)*mu;
    end
    
end
