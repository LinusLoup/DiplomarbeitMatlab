%% Analytische Lösung (?) ohne Last (f=0) auf einem Kreis:
% 
% a = 0.3;
% b = 1;
% xhat = 1-sqrt(1-a/b);
% alpha = (a-b*xhat^2)/(1-xhat);
% u_ana = zeros(np,1);
% for i = 1:np
%     r = sqrt(p(1,i)^2+p(2,i)^2);
%    if r >= -xhat && r <= xhat
%        u_ana(i) = a-b*r^2;
%    elseif r < -xhat
%        u_ana(i) = alpha*(1+r);
%    else
%        u_ana(i) = alpha*(1-r);
%    end
% end
% 
% subplot(3,2,[5,6]); pdeplot(p,e,t,'zdata',u_ana);
% 
% fval
% 1/2*u_quadprog'*A*u_quadprog - f'*u_quadprog
% 1/2*u_jacobi'*A*u_jacobi - f'*u_jacobi
% 1/2*u_ana'*A*u_ana