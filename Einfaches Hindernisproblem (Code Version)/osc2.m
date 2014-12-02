function [osc2_val,osc2_vec] = osc2(Nplusplus_set,N0minus_set, nodes,triangles,edges,f_load)
%OSC2 evaluates the oscillation osc_2(u_S,psi,f).

% initialization:
np = length(nodes);
osc2_vec = zeros(np,1);
N_set = 1:np;
N_sum2 = setdiff(N_set,union(Nplusplus_set,N0minus_set));

% evaluations of h_P for the two sums of osc_2:
    function h_P = eval_hP(set,global_nodes,global_edges)
        % initialization:
        h_P = zeros(size(set));
        
        for l = 1:length(h_P)
        % computing of the edges, which contain p:
        [~,E_p] = find(global_edges(3:4,:)==set(l));
    
        % calculation of the lengths of E in E_p:
        points_1 = global_nodes(:,global_edges(3,E_p));
        points_2 = global_nodes(:,global_edges(4,E_p));
        points_diff = points_1 - points_2;
        h_P(l) = max(sqrt(points_diff(1,:).^2 + points_diff(2,:).^2));
        end
    end
h_P1 = eval_hP(Nplusplus_set,nodes,edges);
h_P2 = eval_hP(N_sum2,nodes,edges);

% evaluating of the sums:
for i = 1:length(Nplusplus_set)
    int_h = 0;
    % support of phi_p, that is the indices of the triangles:
    [~,w_p] = find(triangles(1:3,:)==Nplusplus_set(i));
    
    % computing of the mean-values f_p:
    area_wp = 0;
    intfP_h = 0;
    for k = 1:length(w_p)
        % local points of the considered triangle and the Jacobian:
        mypoi = nodes(:,triangles(1:3,w_p(k)));
        x = mypoi(1,:);
        y = mypoi(2,:);
        J = (x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1));
        % evaluating the Gauss-points for the quadrature:
        [wi,gauss,~] = quad_tri(mypoi,@(x,y) 0,7);
        % quadrature:
        intfP_h = intfP_h + J * sum(wi.*f_load(gauss(1,:),gauss(2,:)).^2);
        % computing of |w_p|:
        area_wp = area_wp + J/2;
    end    
    f_P = 1/area_wp*intfP_h;
    
    % calculation of the norm:
    for j = 1:length(w_p)
        % local points of the considered triangle and the Jacobian:
        mypoi = nodes(:,triangles(1:3,w_p(j)));
        x = mypoi(1,:);
        y = mypoi(2,:);
        J = (x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1));
        % evaluating the Gauss-points for the quadrature:
        [wi,gauss,~] = quad_tri(mypoi,@(x,y) 0,7);
        % quadrature:
        int_h = int_h + J * sum(wi.*(f_load(gauss(1,:),gauss(2,:))...
            -f_P).^2);
    end
    
    osc2_vec(Nplusplus_set(i)) = osc2_vec(Nplusplus_set(i))+h_P1(i)^2*int_h;
end

for i = 1:length(N_sum2)
    int_h = 0;
    % support of phi_p, that is the indices of the triangles:
    [~,w_p] = find(triangles(1:3,:)==N_sum2(i));
    
    % calculation of the norm:
    for j = 1:length(w_p)
        % local points of the considered triangle and the Jacobian:
        mypoi = nodes(:,triangles(1:3,w_p(j)));
        x = mypoi(1,:);
        y = mypoi(2,:);
        J = (x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1));
        % evaluating the Gauss-points for the quadrature:
        [wi,gauss,~] = quad_tri(mypoi,@(x,y) 0,7);
        % quadrature:
        int_h = int_h + J * sum(wi.*f_load(gauss(1,:),gauss(2,:)).^2);
    end
    
    osc2_vec(N_sum2(i)) = osc2_vec(N_sum2(i))+h_P2(i)^2*int_h;
end

osc2_val = sqrt(sum(osc2_vec));
end