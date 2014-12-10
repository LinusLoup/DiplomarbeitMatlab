function j_E = normal_flux(E_p,given_flag,nodes,triangles,edges, edge_triangles,uS_values,lambda,mu)
%NORMAL_FLUX computes the stress normal flux for all given edges E in E_p. 
%The other given objects are: The nodes of the triangulation, the indices 
%of the triangles-nodes ordering, the midpoints/edges-matrix edges, the 
%edges-triangle-ordering edge_triangle and the functionvalues of the 
%Galerkin solution uS_values. Also given are the Lame-constants.

% Initializing:
j_E = zeros(2,size(E_p,2));
    
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
    p_T1 = p_T(:,1:3);
    p_T2 = p_T(:,4:6);
    
        % The calculated u_S values onto T_1 and T_2, where the first three 
        % rows are x-coordinates and the last three the y-coordinates:
    uS_T1 = uS_values([2*neigh_tri(:,1)-1,2*neigh_tri(:,1)]);
            
    if flag(k) == 0
        uS_T2 = zeros(6,1);
    else
        uS_T2 = uS_values([2*neigh_tri(:,2)-1,2*neigh_tri(:,2)]);
    end
    
    graduS_T_x = [gradu(p_T1,uS_T1(1:3));gradu(p_T2,uS_T2(1:3))];
    graduS_T_y = [gradu(p_T1,uS_T1(4:6));gradu(p_T2,uS_T2(4:6))];
    
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
    
    strain1 = [graduS_T_x(1,1), 1/2*(graduS_T_x(1,2)+graduS_T_y(1,1));
              1/2*(graduS_T_x(1,2)+graduS_T_y(1,1)), graduS_T_y(1,2)];
    strain2 = [graduS_T_x(2,1), 1/2*(graduS_T_x(2,2)+graduS_T_y(2,1));
              1/2*(graduS_T_x(2,2)+graduS_T_y(2,1)), graduS_T_y(2,2)];
    
    if (given_flag(k) > 0)
        stress1 = 2*mu(1)*strain1 + lambda(1)*trace(strain1)*eye(2);
        stress2 = 2*mu(1)*strain2 + lambda(1)*trace(strain2)*eye(2);
    else
        stress1 = 2*mu(2)*strain1 + lambda(2)*trace(strain1)*eye(2);
        stress2 = 2*mu(2)*strain2 + lambda(2)*trace(strain2)*eye(2);
    end
    
    j_E(:,k) = stress2*n-stress1*n;
end

end