function [osc1_val,osc1_vec] = osc1(N0plus_set,obstacle_values,nodes,triangles,uS_values)
%OSC1 evaluates the obstacle oscillation osc_1(u_S,phi). The general condition in this case is that psi is affine.

% Initializing:
osc1_vec = zeros(length(nodes),1);

% Evaluation of the sum over the points of N^{0+}:
for i = 1:length(N0plus_set)
    % the help value int_h for the integral of the gradient:en:
    int_h = 0;
    % evaluating the support of the shape funktions and the indix of the local shape function:
    [~,w_p] = find(triangles(1:3,:)==N0plus_set(i));
    
    % computation of the single norms (integrals):
    for j = 1:length(w_p)
        % nodes of the triangle, Jacobian and factors for the affine trafo:
        p_index = triangles(1:3,w_p(j));
        mypoi = nodes(:,p_index);
        x = mypoi(1,:);
        y = mypoi(2,:);
        J = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));
        a = 1/J * ((x(3)-x(1))^2 + (y(3)-y(1))^2);
        b = -1/J * ((y(2)-y(1))*(y(3)-y(1)) + (x(2)-x(1))*(x(3)-x(1)));
        c = 1/J * ((x(2)-x(1))^2 + (y(2)-y(1))^2);
        
        % evaluate the gradient of psi and u_S:
        grad_psi = gradu(mypoi,obstacle_values(p_index));
        grad_u = gradu(mypoi,uS_values(p_index));
        
        % transformation of the gradients onto the reference element:
        grad_psi = [x(2)-x(1), y(2)-y(1); x(3)-x(1), y(3)-y(1)]*grad_psi';
        grad_u = [x(2)-x(1), y(2)-y(1); x(3)-x(1), y(3)-y(1)]*grad_u';
            
            % to integrate the gradients, we just multiplicate the area of the reference element with the high (the functionvalue), because the gradients are constant:
            int_h = int_h + 1/2*(a*grad_psi(1)^2+2*b*grad_psi(1)...
                *grad_psi(2)+c*grad_psi(2)^2-2*(a*grad_psi(1)*grad_u(1)...
                +b*(grad_psi(1)*grad_u(2)+grad_psi(2)*grad_u(1))...
                +c*grad_psi(2)*grad_u(1))+a*grad_u(1)^2+2*b*(grad_u(1)...
                *grad_u(2))+c*grad_u(2)^2);
    end
    
    osc1_vec(N0plus_set(i)) = int_h;
end

% square root of the sum of the integrals:
osc1_val = sqrt(sum(osc1_vec));

end