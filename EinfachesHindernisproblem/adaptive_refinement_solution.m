function [u_S,points,edges,triangles,midtri,midpoints,rhoS_plot,IQ_plot, J_error,osc_term,recursion_depth] = adaptive_refinement_solution (points,edges,triangles,solution,load_fun,geo_data,J_exact, max_error,para_rho,para_osc,max_points,max_recursion)
%ADAPTIVE_REFINEMENT_SOLUTION uses the adaptive refinement strategy shown
%in chapter 4 and evaluates the solution on a adaptive refined mesh.

% initializing all local values:
u_S = solution;
J_u = J_exact;
rhoS_plot = zeros(max_recursion,1);
IQ_plot = zeros(max_recursion,1);  
J_error = zeros(max_recursion,1);  
osc_term = zeros(max_recursion,1); 
nodes_vec = zeros(max_recursion,1); 
recursion_depth = 1;               

while 1
    % further initializations:
    np = size(points,2);
    nodes_vec(recursion_depth) = np;
    
    % evaluation of the midpoints:
    [midpoints,midtri] = midpoints_of_triangle(triangles,points);
    nmp = size(midpoints,2);

    % computation of the functionvalues of the obstacle onto the mesh nodes and midpoints:
    obstacle = @(x,y) -zeros(size(x))';
    z_obs_prob = obstacle(points(1,:),points(2,:));
    z_obs_midpoints = obstacle(midpoints(1,:),midpoints(2,:));

    % assembling of the matrix/vector data and the evaluation of the Dirichlet boundary data: 
    [A,f] = assemble(points,triangles,load_fun,7,'linear');
    [~,~,H,R]=assemb(geo_data.mysquareb,points,edges);

    % solution of the variational inequality with active-set/inner-points-method:
    u_S = sparse(u_S);
    z_obs_prob = sparse(z_obs_prob);
    opts = optimset('Algorithm','interior-point-convex','LargeScale', 'on','Display','off');
    [u_S,J_uS] = quadprog(A,-f,[],[],H,R,z_obs_prob,[],u_S,opts);

    % computing the functionvalues of u_S onto the midpoints:
    u_S_mid = zeros(nmp,1);
    
    for j = 1 : nmp
        index = midpoints([3,4],j);
        u_S_mid(j) = (u_S(index(1))+u_S(index(2)))/2; 
    end
    
    % the solution of the local defect problem by assembling the matrix with the bubble functions and using the equations in (4.10) and (4.11):
    [A_Q,rhoS_phiE] = assemble(points,triangles,load_fun,7,'bubble',u_S);
    [eps_V,rho_E,d_E,~] = defect_problem_solution(points,triangles, midtri,rhoS_phiE,u_S_mid,z_obs_midpoints);

    % the hierarchical error estimator: evaluation of rho_S(eps_V) with the equation beyond (4.68):
    rhoS_glob = eps_V'*rhoS_phiE;
    rhoS_plot(recursion_depth) = rhoS_glob;
    
    % computation of -I_Q(eps_V) with (4.6):
    IQ_plot(recursion_depth) = -1/2*(eps_V.^2)'*diag(A_Q) + rhoS_glob;
    
    % calculation of the error of -I(e):
    J_error(recursion_depth) = J_uS-J_u;

    % evaluating the local contributions of rho_S with Lemma 4.14:
    rho_p = eval_rho_p(points,triangles,midpoints,midtri,u_S,eps_V, load_fun);

    % determination of the sets N0, N0+, N+, N++, N0- for the oscillationterms:
    inner_points_omega = inner_points(H);
    N0_set = N0(u_S,inner_points_omega,z_obs_prob);
    Nplus_set = Nplus(N0_set,inner_points_omega,points);
    [N0plus_set,N0minus_set] = N0plusminus(N0_set,points, triangles,midpoints,midtri,load_fun,obstacle,u_S);
    Nplusplus_set = Nplusplus(Nplus_set,midpoints,rho_E,d_E);
    
    % evaluation of the oscillationterms:
    [osc1_term,osc1_local] = osc1(N0plus_set,z_obs_prob,points, triangles,u_S);
    [osc2_term,osc2_local] = osc2(Nplusplus_set,N0minus_set,points, triangles,midpoints,load_fun);
    osc_local = osc1_local + osc2_local;
    osc_term(recursion_depth) = sqrt(osc1_term^2 + osc2_term^2);

    % calculating the indices of the triangles, which have to be refined:
    refine_triangle = find_triangle_refinement(rho_p,rhoS_glob, osc_local,osc_term(recursion_depth),triangles,para_rho, para_osc,'symmetric');
   
    % refinement of the mesh:
    [p_h,e_h,t_h,uS_h] = refinemesh(geo_data.mysquareg,points,edges, triangles,u_S,refine_triangle);
    
    % termination criterion: if the number of the nodes is too large, no triangle will be refined or the hierarchical error estimator is small enough:
    if (isempty(refine_triangle) || rhoS_glob < max_error || recursion_depth == max_recursion || length(p_h) > max_points)
        length(p_h);
        fprintf('%s %f.\n','Die Rekursionstiefe ist ', recursion_depth);
        break;
    else
        points = p_h;
        edges = e_h;
        triangles = t_h;
        u_S = uS_h;
        recursion_depth = recursion_depth + 1;
    end
end

% elimination of the zeros in the vectors of the error/-estimator:
rhoS_plot = rhoS_plot(1:recursion_depth);
IQ_plot = IQ_plot(1:recursion_depth);
J_error = J_error(1:recursion_depth);
osc_term = osc_term(1:recursion_depth);

end