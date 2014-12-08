function [u_S,points,edges,triangles,midtri,midpoints,rhoS_plot,IQ_plot, J_error,recursion_depth, degree_of_freedom,time_vec] = adaptive_refinement_solution (points,edges,triangles,lambda,mu,solution,vol_load_fun_x, vol_load_fun_y,surf_load_fun_x,surf_load_fun_y,obstacle,geo_data, J_exact,max_error,para_rho,max_points,max_recursion)
%ADAPTIVE_REFINEMENT_SOLUTION uses the adaptive refinement strategy shown
%in chapter 4 and evaluates the solution on a adaptive refined mesh.

% initializing all local values:
u_S = solution;
J_u = J_exact;
rhoS_plot = zeros(max_recursion,1);
IQ_plot = zeros(max_recursion,1);  
J_error = zeros(max_recursion,1);  
recursion_depth = 1;
time_vec = zeros(max_recursion,1);
degree_of_freedom = zeros(max_recursion,1);

while 1
    tic
    % further initializations:
    np = size(points,2);
    degree_of_freedom(recursion_depth) = np;

    % evaluate the boundary indices and conditions:
    [~,~,H,~]=assemb(geo_data.mygeomb,points,edges);
    [~,boundary_ind] = find(H(:,1:np));
    [gamma_D,gamma_N,gamma_C] = boundary_disjunction(boundary_ind,points);
    
    num_of_dirichlet_points = length(gamma_D);
    adjacency_D = sparse(2*num_of_dirichlet_points,2*np);
    for i = 1:num_of_dirichlet_points
        adjacency_D(2*i-1:2*i,2*gamma_D(i)-1:2*gamma_D(i)) = eye(2);
    end
    dirichlet_bound = sparse(2*(num_of_dirichlet_points),1);
    
    num_of_contact_points = length(gamma_C);
    B = zeros(num_of_contact_points,2*np);
        % for simplicity just n = (0,1):
    for j = 1:num_of_contact_points
        if points(2,gamma_C(j))>0
            B(j,2*gamma_C(j)-1:2*gamma_C(j)) = [0,-1];
        else
            B(j,2*gamma_C(j)-1:2*gamma_C(j)) = [0,1];
        end
    end
    
    % assembling of the matrix/vector data and the evaluation of the Dirichlet boundary data: 
    [A,f] = assemble(points,triangles,vol_load_fun_x,vol_load_fun_y, surf_load_fun_x,surf_load_fun_y,gamma_N,lambda,mu,7,'linear');
    
    % evaluation of the midpoints:
    [midpoints,midtri] = midpoints_of_triangle(triangles,points);
    nmp = size(midpoints,2);

    % computation of the functionvalues of the obstacle onto the mesh nodes and midpoints:
    z_obs_prob = obstacle(points(1,gamma_C),points(2,gamma_C));
    
    contact_midpoints = zeros(5,num_of_contact_points-1);
    contact_count = 1;
    for k = 1:nmp
        if (ismember(midpoints(3,k),gamma_C) && ismember(midpoints(4,k),gamma_C))
           contact_midpoints(1:4,contact_count) = midpoints(:,k);
           contact_midpoints(5,contact_count) = k;
           contact_count = contact_count+1;
        end   
    end
    z_obs_midpoints = obstacle(contact_midpoints(1,:),contact_midpoints(2,:));
    
    if size(z_obs_prob,1) == 1
        z_obs_prob = z_obs_prob';
        z_obs_midpoints = z_obs_midpoints';
    end

    % solution of the variational inequality with active-set/inner-points-method:
    u_S = sparse(u_S);
    z_obs_prob = sparse(z_obs_prob);
    opts = optimset('Algorithm','interior-point-convex','LargeScale', 'on','Display','off');
    [u_S,J_uS] = quadprog(A,-f,B,z_obs_prob,adjacency_D,dirichlet_bound, [],[],u_S,opts);
    
    % plot of the deformed cylinder:
    point_shift = zeros(2,np);
    point_shift(1:2*np) = u_S;
    new_points = points+point_shift;
    figure(recursion_depth+1)
    pdeplot(new_points,edges,triangles)
    
    % computing the functionvalues of u_S onto the midpoints:
    u_S_mid = zeros(num_of_contact_points-1,1);
    
    for j = 1 : num_of_contact_points-1
        index = contact_midpoints([3,4],j);
         
        if contact_midpoints(2,j)>0
            u_S_mid(j) = [(u_S(2*index(1)-1)+u_S(2*index(2)-1))/2, (u_S(2*index(1))+u_S(2*index(2)))/2]*[0; -1/2];
        else
            u_S_mid(j) = [(u_S(2*index(1)-1)+u_S(2*index(2)-1))/2, (u_S(2*index(1))+u_S(2*index(2)))/2]*[0; 1/2];
        end
    end
    
    % the solution of the local defect problem by assembling the matrix with the bubble functions and using the equations in (4.10) and (4.11):
    [A_Q,rhoS_phiE] = assemble(points,triangles,vol_load_fun_x,vol_load_fun_y, surf_load_fun_x,surf_load_fun_y,gamma_N,lambda,mu,7,'bubble', u_S);
    %[eps_V,rho_E,d_E,~] = defect_problem_solution(points,triangles, midtri,rhoS_phiE,u_S_mid,z_obs_midpoints);
    B_Q = zeros(num_of_contact_points-1,2*nmp);
        % for simplicity just n = (0,1):
    for j = 1:num_of_contact_points-1
        if contact_midpoints(2,j)>0
            B_Q(j,2*contact_midpoints(5,j)-1:2*contact_midpoints(5,j)) = [0,-1/2];
        else
            B_Q(j,2*contact_midpoints(5,j)-1:2*contact_midpoints(5,j)) = [0,1/2];
        end
    end
    [eps_V,~] = quadprog(A_Q,-rhoS_phiE,B_Q,z_obs_midpoints+u_S_mid, [],[],[],[],[],opts);
    
    % the hierarchical error estimator: evaluation of rho_S(eps_V) with the equation beyond (4.68):
    rhoS_glob = eps_V'*rhoS_phiE;
    rhoS_plot(recursion_depth) = rhoS_glob;

    % computation of -I_Q(eps_V) with (4.6):
    IQ_plot(recursion_depth) = -1/2*(eps_V.^2)'*diag(A_Q) + rhoS_glob;
    
    % calculation of the error of -I(e):
    J_error(recursion_depth) = J_uS-J_u;

    % evaluating the local contributions of rho_S with Lemma 4.14:
%    rho_p = eval_rho_p(points,triangles,edges,midpoints,midtri,u_S, eps_V,vol_load_fun_x, vol_load_fun_y,surf_load_fun_x,surf_load_fun_y, gamma_N,lambda,mu);
    
    % calculating the indices of the triangles, which have to be refined:
%    refine_triangle = find_triangle_refinement(rho_p,rhoS_glob,triangles,para_rho,'symmetric');
    
    % refinement of the mesh:
    [p_h,e_h,t_h,uS_h] = refinemesh(geo_data.mygeomg,points,edges, triangles,u_S);%,refine_triangle);
    
    % termination criterion: if the number of the nodes is too large, no triangle will be refined or the hierarchical error estimator is small enough:
    if ( rhoS_glob < max_error || recursion_depth == max_recursion || length(p_h) > max_points)
        u_S = point_shift;
        fprintf('%s %f.\n','Die Rekursionstiefe ist ', recursion_depth);
        break;
    else
        points = p_h;
        edges = e_h;
        triangles = t_h;
        u_S = uS_h;
        recursion_depth = recursion_depth + 1;
    end
    time_vec(recursion_depth) = toc;
end

% elimination of the zeros in the vectors of the error/-estimator:
rhoS_plot = rhoS_plot(1:recursion_depth);
IQ_plot = IQ_plot(1:recursion_depth);
J_error = J_error(1:recursion_depth);
degree_of_freedom = degree_of_freedom(1:recursion_depth);
end