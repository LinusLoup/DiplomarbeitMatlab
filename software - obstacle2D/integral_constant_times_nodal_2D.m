function [integral_value, integral_density]=integral_constant_times_nodal_2D(sigma,u,elements,areas) 
        
        u_average=sum(u(elements),2)/size(elements,2);
        integral_density=(sigma.*u_average);
        integral_value=integral_density'*areas;
        
        

        
end