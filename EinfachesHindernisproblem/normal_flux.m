function j_E = normal_flux(E_p,nodes,triangles,edges,edge_triangles,uS_values)
%NORMAL_FLUX computes the normal flux for all given edges E \in E_p. The other given objects are: The nodes of the triangulation, the indices of the triangles-nodes ordering, the midpoints/edges-matrix edges, the edges-triangle-ordering edge_triangle and the functionvalues of the Galerkin solution uS_values.

% Initializing:
j_E = zeros(size(E_p));
    
% Evaluate the neigbours of all edges in the edges-set E_p:
[neighbours,flag] = neighbourhood(E_p,edge_triangles,'edges');
        
% Computation of the single normal fluxes:
for k = 1:length(E_p)     
    % Evaluation of the edge-points:
    edge_poi_ind = edges(3:4,E_p(k));
    edge_poi = nodes(:,edge_poi_ind);
            
    % Calculation of the gradient of u_S on T_1 and T_2:
    neigh_tri_ind = neighbours(:,k);
    neigh_tri = triangles(1:3,neigh_tri_ind);
    p_T = nodes(:,neigh_tri);
    uS_T = uS_values(neigh_tri);
            
    p_T1 = p_T(:,1:3);
    p_T2 = p_T(:,4:6);
    uS_T1 = uS_T(:,1);
            
    if flag(k) == 0
        uS_T2 = zeros(3,1);
    else
        uS_T2 = uS_T(:,2);
    end
    
    graduS_T = [gradu(p_T1,uS_T1);gradu(p_T2,uS_T2)];
    
    % the normalvector from T_1 to T_2:
    connect = edge_poi(:,2)-edge_poi(:,1);
    orth_connect = [-connect(2);connect(1)];
    n = 1/norm(orth_connect)*orth_connect;
            
    % Verification, if n shows from T_1 to T_2:
    p3_T1 = setdiff(p_T1',edge_poi','rows')';
    edge_test = (p3_T1-edge_poi(:,1))';
            
    if edge_test*n > 0
        n = -n;
    end
    
    normal = graduS_T*n;
    j_E(k) = normal(2)-normal(1);
end

end