function rho_p = eval_rho_p(nodes,triangles,midpoints,midtri,u_S,eps_V,fun)
%EVAL_RHO_P evaluates the local contribution rho_P of rho_S, given the nodes, triangles, midpoints of the edges (of the triangles) and die midpoints-triangle-ordering (midtri). Also it is given the solution auf die local defect problem eps_V and the function fun, which is part of the integral in rho_S. u_S are, as always, the Galerkin solution.

% Initializing:
np = size(nodes,2);
rho_p = zeros(np,1);

% Hat- and Bubble-functions:
phi_P = @(xi,eta) [1-xi-eta; xi; eta];
phi_E = @(xi,eta) [4*xi.*(1-xi-eta); 4*xi.*eta; 4*eta.*(1-xi-eta)];

% Evaluation of the weights and values of the function for the surface integral:
[wi,~,phi_E_values] = quad_tri([0,1,0;0,0,1],phi_E,7);
[~,~,phi_P_values] = quad_tri([0,1,0;0,0,1],phi_P,7);
        
% local weights and nodes for the Gauss quadrature for the line integral:
nodes_local = [-sqrt(1/3),sqrt(1/3)];
phiPE_local = @(xi) (xi+1)/2.*(-xi.^2+1);
local_phiPE_int = sum(phiPE_local(nodes_local));

for k = 1:np
    % Ordering triangle <-> reference function and support of phi_P:
    [phi_p_local,w_p] = find(triangles(1:3,:)==k);
    % indices of the edges in the considered triangle:
    E_index = midtri(:,w_p);
    
    % Evaluation of the surface integral of rho_p over w_p:
    for j = 1:length(w_p)
        % points of the considered triangle with Jacobi determinant:
        mypoi = nodes(:,triangles(1:3,w_p(j)));
        x = mypoi(1,:);
        y = mypoi(2,:);
        J = (x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1));
        % Gauss-points and local contributions of eps_V:
        eps_V_local = eps_V(E_index(:,j));
        [~,gauss,~] = quad_tri(mypoi,@(x,y) 0,7);
        % the functionvalues of the considered local hatfunction phi_P:
        phi_Pl_values = phi_P_values(phi_p_local(j),:);
        % evaluation of the first integral of rho_p:
        rho_p(k) = rho_p(k) + J * sum(wi.*fun(gauss(1,:),gauss(2,:)).*...
            (eps_V_local'*phi_E_values).*(phi_Pl_values));
    end
    
    % Evaluation of the line integral over the edges E\in E_p:
        % Evaluaton of the set E_p of the edges, which contain the point P:
        [~,E_p] = find(midpoints(3:4,:)==k);
        % Evaluation of the normal-fluxes j_E for all edges in E_p:
        j_E = normal_flux(E_p,nodes,triangles,midpoints,midtri,u_S);

        for i = 1:length(E_p)    
            % evaluating the points of the edges:
            edge_poi_ind = midpoints(3:4,E_p(i));
            edge_poi = nodes(:,edge_poi_ind);
            
            % affin transformation of the local integral on the global edge by multiplication of the jacobian:
            laenge = norm(edge_poi(:,1)-edge_poi(:,2));
            global_phiPE_int = local_phiPE_int*1/2*laenge;
            
            % determination of the functionvalues of eps_V(x_E):
            epsV_loc = eps_V(E_p(i));
            
            % adding the line integral to the local contribution rho_p:
            rho_p(k) = rho_p(k)+j_E(i)*epsV_loc*global_phiPE_int;
        end
end

end