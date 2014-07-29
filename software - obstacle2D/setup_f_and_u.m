function [u,uGradient,f,lambda,phi,u_H1seminorm_squared,energy_exact,r_contact]=setup_f_and_u(benchmark,coordinates,f_amplitude,phi_amplitude,Lx,Ly)
switch benchmark
    case 'exact_ring_constant_obstacle'
        f=f_amplitude*ones(size(coordinates,1),1);  
           
        phi=phi_amplitude*ones(size(coordinates,1),1); 
        [u,uGradient,lambda,energy_exact,u_H1seminorm_squared,r_contact]=u_exact_ring_constant_obstacle(coordinates,f_amplitude,phi_amplitude,Lx,Ly);    
                  
    case 'exact_ring_spherical_obstacle'
        rho=1.2;    
        
        f=f_amplitude*ones(size(coordinates,1),1);  
           
        xy=coordinates;
        distance=sqrt(sum((xy-[ones(size(xy,1),1)*Lx/2 ones(size(xy,1),1)*Ly/2]).^2,2));
        phi=(phi_amplitude-rho)*ones(size(coordinates,1),1); 
        
        index=find(distance<=rho);    %above spherical obstacle   
        phi(index)=phi(index)+sqrt(rho^2-distance(index).^2);     
        
        [u,uGradient,lambda,energy_exact,u_H1seminorm_squared,r_contact]=u_exact_ring_spherical_obstacle(coordinates,f_amplitude,phi_amplitude,Lx,Ly,rho);    
               
    case 'exact_square_constant_obstacle'                  
        [u,uGradient,lambda,energy_exact,u_H1seminorm_squared,r_contact,f,phi]=u_exact_square_constant_obstacle(coordinates);
        
    case 'exact_square_no_obstacle'
        fSquare=@(x)(f_amplitude*2)*[Lx*x(:,1)-x(:,1).*x(:,1)+Ly*x(:,2)-x(:,2).*x(:,2)]; 
        exactSolutionSquare=@(x)(f_amplitude)*(Lx*x(:,1)-x(:,1).*x(:,1)).*(Ly*x(:,2)-x(:,2).*x(:,2));
        exactGradientSquare=@(x)(f_amplitude)*[(Lx-2*x(:,1)).*(Ly*x(:,2)-x(:,2).*x(:,2)) ...
                                          (Ly-2*x(:,2)).*(Lx*x(:,1)-x(:,1).*x(:,1))];    
        
        u=exactSolutionSquare(coordinates);
        uGradient=exactGradientSquare(coordinates);
        f=fSquare(coordinates);
        lambda=zeros(size(coordinates,1),1);      
   
        energy_exact=-(f_amplitude^2*Lx^3*Ly^3*(Lx^2+Ly^2))/180;
        r_contact=0;
        phi=phi_amplitude*ones(size(coordinates,1),1); 
        
        u_H1seminorm_squared=(f_amplitude^2*Lx^3*Ly^3*(Lx^2+Ly^2))/90;        
end

end
