function [N0plus_set,N0minus_set] = N0plusminus(N0_set,nodes, triangles,edges,edge_triangles,f_load,obstacle,uS_values)
%N0PLUSMINUS evaluates the sets N^{0+} and N^{0-} for the calculation of the oscillationterms, given the inner contact nodes N^0_set, the nodes of the mesh, the triangle ordering, the edges and the midpoints-triangle ordering, the loadfunction, the obstacle and vector of u_S values.

% parameter for the approximate of zero:
zero = 0.0001;

% initializing of the sets N^{0+} and N^{0-}:
N0plus_set = zeros(size(N0_set));
N0minus_set = zeros(size(N0_set));

% calculating the points of the reference element to compare the functionvalues later:
[xi,eta] = meshgrid(0:0.01:1,0:0.01:1);
index_out_of_bounds = find(xi+eta>1);
xi(index_out_of_bounds) = 0;
eta(index_out_of_bounds) = 0;
[m_mesh,n_mesh] = size(xi);
np = m_mesh * n_mesh;
index_within = setdiff(1:np,index_out_of_bounds);

% hat functions with functionvalues:
phi_P = @(xi,eta) [1-xi-eta; xi; eta];
phi_P_values = phi_P(xi,eta);

% evaluation of the indices of the N^{0+} nodes:
for i = 1:length(N0_set)
    % support of the shape function and index of the local shape function:
    [phi_p_local,w_p] = find(triangles(1:3,:)==N0_set(i));
    flag_plus = 0;      % verification value for the condition of N0+
    flag_minus = 0;     % verification value for the condition of N0-
    
    for j = 1:length(w_p)
        % vertices of the triangle and functionvalues of u_S:
        p_index = triangles(1:3,w_p(j));
        uS_local = uS_values(p_index);
        mypoi = nodes(:,p_index);
        x = mypoi(1,:);
        y = mypoi(2,:);
        
        % transformation from the local onto the global triangle:
        xval = x(1)+(x(2)-x(1)).*xi+(x(3)-x(1)).*eta;
        yval = y(1)+(y(2)-y(1)).*xi+(y(3)-y(1)).*eta;

        % evaluating the functionvalues of u_S-psi:
        z = uS_local(1)*phi_P_values(1:m_mesh,:)+uS_local(2)*...
            phi_P_values(1+m_mesh:2*m_mesh,:)+uS_local(3)*...
            phi_P_values(1+2*m_mesh:3*m_mesh,:)-obstacle(xval,yval);
        
        % computing the functionvalues of the load f_load:
        f = f_load(xval,yval);
        
        % verificating, if u_S-psi>0 and if u_S=psi & f<=0:
        switch phi_p_local(j)
            case 1
                new_within = setdiff(index_within,1);
            case 2
                new_within = setdiff(index_within,sub2ind([m_mesh,...
                    n_mesh],1,n_mesh));
            case 3
                new_within = setdiff(index_within,m_mesh);
        end
        
        % If any z-value is less than zero, the point cannot be an element of N^{0+}, that is, why flag_plus = 1. If all f_load-values and all z-values are less than zero, p could be in N^{0-}:
        if any(z(new_within) <= zero)
            flag_plus = 1;
            
            z_nodes = uS_local-obstacle(x,y)';
            if (all(abs(z_nodes) <= zero) && all(f(index_within)<=0))
                flag_minus = flag_minus + 1;
            end
        end
    end
    
    % verification if there is no contact except of the point p:
    if flag_plus == 0
        N0plus_set(i) = N0_set(i);
        continue;
    end
    
    % verification of j_E <= 0 for all E in E_p: evaluating of the set of edges E_p, which contains the point p:
    [~,E_p] = find(edges(3:4,:)==N0_set(i));
    % evaluating the normal fluxes of all edges of E_p:
    j_E = normal_flux(E_p,nodes,triangles,edges,edge_triangles,uS_values);
    % verificating, if there is full contact and f<=0 and j_E<=0 for all E in E_p:
    if (flag_minus == length(w_p) && all(j_E <= zero))        
        N0minus_set(i) = N0_set(i);
    end
end

% elimination of the zeros:
N0plus_set = setdiff(N0plus_set,0);
N0minus_set = setdiff(N0minus_set,0);

end