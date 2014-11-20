function exact_solution_obstacle_plot

figure(1);

[x,y] = meshgrid(-1.5:0.04:1.5,-1.5:0.04:1.5);
z = zeros(size(x));
index = find(x.^2+y.^2>=1);
z(index) = 1/2*(x(index).^2+y(index).^2)-1/2*log(x(index).^2+y(index).^2)-1/2;
surf(x,y,z);

figure(2);

[x,y] = meshgrid(-2:0.04:2,-2:0.04:2);
    function z = my_fun(xn,yn)
        [m,n] = size(xn);
        z = zeros(m,n);
        r = sqrt(xn.^2+yn.^2);
        r_new = 2*r-1/2;
        
        index_gamma1 = find(r>=1/4 & r<3/4 & yn>=0);
        
        z(index_gamma1) = (r(index_gamma1)).^(2/3)...
            .*sin(2/3*atan2(yn(index_gamma1),xn(index_gamma1)))...
            .*(-6*r_new(index_gamma1).^5+15*r_new(index_gamma1).^4 ...
            -10*r_new(index_gamma1).^3+1);
        
        index_gamma1 = find(r>=1/4 & r<3/4 & yn<0);
        
        z(index_gamma1) = (r(index_gamma1)).^(2/3)...
            .*sin(2/3*(atan2(yn(index_gamma1),xn(index_gamma1))+2*pi))...
            .*(-6*r_new(index_gamma1).^5+15*r_new(index_gamma1).^4 ...
            -10*r_new(index_gamma1).^3+1);
            
        index_gamma2 = find(r<1/4 & yn>=0);
        z(index_gamma2) = (r(index_gamma2)).^(2/3)...
            .*sin(2/3*atan2(yn(index_gamma2),xn(index_gamma2)));
        
        index_gamma2 = find(r<1/4 & yn<0);
        z(index_gamma2) = (r(index_gamma2)).^(2/3)...
            .*sin(2/3*(atan2(yn(index_gamma2),xn(index_gamma2))+2*pi));
        
        index3 = find(xn>=0 & yn<=0);
        z(index3) = 0;
    end
my_fun(0,-1)
u = my_fun(x,y);
surf(x,y,u);

end