clear all
close all

addpath('.\library_vectorization');

benchmark='exact_ring_constant_obstacle';   
%'exact_ring_constant_obstacle', 'exact_ring_spherical_obstacle',
%'exact_square_constant_obstacle', 'exact_square_no_obstacle'

iterations_majorant_all=10;  %number of majorant minimization iterations
order_f=0;   %0 or 1 (integrations of rhs f, 0 - piecewise constant f, 1 - piecewise nodal f)

min_level_refinement=2; 
max_level_refinement=5; %levels of mesh refinement 

f_amplitude=-10; phi_amplitude=-1;  %does not work for 'exact_square_constant_obstacle' benchmark

draw=011010111;
draw=11000010;
%draw format:            bit 1 - exact solution, multiplier and flux
%                        bit 2 - discrete_solution and fluxes
%                        bit 3 - loading
%                        bit 4 - exact_error distribution (also recomputed
%                                                          on finer meshes)
%                        bit 5 - multipliers: nodal, rescaled constant
%                        bit 6 - majorant distribution (all three terms separately + total value)
%                        bit 7 - computed flux and multiplier from majorant 
%                        bit 8 - convergence plot (interestiong for more refinements) 
%                        bit 9 - separated figures for publication

Dirichlet_boundary_value_correction='no';  %only for testing of ring  benchmarks
Lx=2;   % beam length in x direction
Ly=2;   % beam length in y direction  
C_Omega=sqrt(2)/pi;    %0.225;
C_Omega=2;    %0.225;


for i=1:numel(iterations_majorant_all)
    
iterations_majorant=iterations_majorant_all(i);

for level=min_level_refinement:max_level_refinement
    close all
    mesh_setup;    %created elementes and coordinates, not edges
    
    if strfind(benchmark,'exact_ring')
        fprintf('difference of areas = %d:\n', (size(elementsRectangular_circumscribed,1)-size(elementsRectangular,1))*hx*hy);    
    end

    fprintf('Rectangulation consists of:\n');
    fprintf('%d elements, ',size(elementsRectangular,1));
    fprintf('%d nodes, ',size(coordinates,1));
    fprintf('%d edges, \n',size(edge2nodes,1));
    fprintf('%d nodes in x direction, ',nx);
    fprintf('%f step size in x direction, \n',hx);
    fprintf('%d nodes in y direction, ',ny);
    fprintf('%f step size in y direction. \n',hy);
    %fprintf('=> %d DOF \n',size(edge2nodes,1));

    %nodal values f and u
    [u,uGradient,f_nodal,dummy,phi,u_H1seminorm_squared,energy_exact,r_contact]=setup_f_and_u(benchmark,coordinates,f_amplitude,phi_amplitude,Lx,Ly);
    [u_finer,dummy, f_nodal_finer,dummy,dummy,u_H1seminorm_squared_finer]=setup_f_and_u(benchmark,coordinates_finer,f_amplitude,phi_amplitude,Lx,Ly);
    [u_finest,dummy, f_nodal_finest,dummy,dummy,u_H1seminorm_squared_finest]=setup_f_and_u(benchmark,coordinates_finest,f_amplitude,phi_amplitude,Lx,Ly);


    %dirichlet value u 
    uDirichlet=zeros(size(u));
    if  isempty(strfind(benchmark,'exact_ring'))
        uDirichlet(nodes_dirichlet)=setup_f_and_u(benchmark,coordinates(nodes_dirichlet,:),f_amplitude,phi_amplitude,Lx,Ly);
    end

    %midpoint values f and u
    [u_midpointRectangular,dummy,f_midpointRectangular,lambda_midpointRectangular,dummy]=setup_f_and_u(benchmark,elementsRectangular2midpoint,f_amplitude,phi_amplitude,Lx,Ly);
    [u_midpointRectangular_finer,dummy,f_midpointRectangular_finer,lambda_midpointRectangular_finer,dummy]=setup_f_and_u(benchmark,elementsRectangular2midpoint_finer,f_amplitude,phi_amplitude,Lx,Ly);
    [u_midpointRectangular_finest,dummy,f_midpointRectangular_finest,lambda_midpointRectangular_finest,dummy]=setup_f_and_u(benchmark,elementsRectangular2midpoint_finest,f_amplitude,phi_amplitude,Lx,Ly);
    [u_midpointTriangular,dummy,f_midpointTriangular,lambda_midpointTriangular,dummy]=setup_f_and_u(benchmark,elementsTriangular2midpoint,f_amplitude,phi_amplitude,Lx,Ly);
    uGradientX=uGradient(:,1);
    uGradientY=uGradient(:,2);

    %assembly of FEM matrices and areas
    [KTriangular,areasTriangular]=stifness_matrixP1_2D(elementsTriangular,coordinates);
    MTriangular=mass_matrixP1_2D(elementsTriangular,areasTriangular);
    [KRectangular,MRectangular, KRectangular_3Dmatrix, MRectangular_3Dmatrix]=assembly_FEM_matrices_rectangle_bilinear(elementsRectangular,hx,hy);
    [KRectangular_finer, dummy, KRectangular_finer_3Dmatrix, dummy]=assembly_FEM_matrices_rectangle_bilinear(elementsRectangular_finer,hx_finer,hy_finer);
    [KRectangular_finest, dummy, KRectangular_finest_3Dmatrix, dummy]=assembly_FEM_matrices_rectangle_bilinear(elementsRectangular_finest,hx_finest,hy_finest);

    areasRectangular=hx*hy*ones(size(elementsRectangular,1),1);

    %assembly of rhs
    switch order_f;
    case 0,  
            f_averageRectangular=f_midpointRectangular;
            f_averageRectangular_finer=f_midpointRectangular_finer;
            f_averageRectangular_finest=f_midpointRectangular_finest;
            f_averageTriangular=f_midpointTriangular;
            
            bTriangular=sparse(elementsTriangular,ones(size(elementsTriangular)),(1/6)*hx*hy*kron(f_averageTriangular,ones(1,3))); 
            [bRectangular, bRectangular_3Dmatrix]=linearForm_fQ0_vQ1(f_averageRectangular,elementsRectangular,hx,hy); 
            [bRectangular_finer, dummy]=linearForm_fQ0_vQ1(f_averageRectangular_finer,elementsRectangular_finer,hx_finer,hy_finer); 
            [bRectangular_finest, dummy]=linearForm_fQ0_vQ1(f_averageRectangular_finest,elementsRectangular_finest,hx_finest,hy_finest); 
                        
    case 1, 
            f_averageRectangular=evaluate_elements_average(elementsRectangular,f_nodal);
            f_averageTriangular=evaluate_elements_average(elementsTriangular,f_nodal);
            
            bTriangular=MTriangular*f_nodal;
            %bRectangular2=MRectangular*f_nodal;  
            
            [bRectangular, bRectangular_3Dmatrix]=linearForm_fQ1_vQ1(f_nodal,elementsRectangular,hx,hy); %using nodal value of f     
            [bRectangular_finer, dummy]=linearForm_fQ1_vQ1(f_nodal_finer,elementsRectangular_finer,hx_finer,hy_finer); %using nodal value of f
            [bRectangular_finest, dummy]=linearForm_fQ1_vQ1(f_nodal_finest,elementsRectangular_finest,hx_finest,hy_finest); %using nodal value of f    
            
    case 2,
            f_averageRectangular=evaluate_elements_average(elementsRectangular,f_nodal);
            f_averageTriangular=evaluate_elements_average(elementsTriangular,f_nodal);
            
            elementsRectangularExtended=[elementsRectangular,(1+max(max(elementsRectangular)):size(elementsRectangular,1)+max(max(elementsRectangular)))'];
            f_bubbleRectangular=f_midpointRectangular-sum(f_nodal(elementsRectangular),2)/4;
            f_nodalRectangularExtended=[f_nodal; f_bubbleRectangular];
            [bRectangular, bRectangular_3Dmatrix]=linearForm_fQ1plusbubble_vQ1(f_nodalRectangularExtended,elementsRectangularExtended,hx,hy);               
            [bRectangular2, bRectangular_3Dmatrix2]=linearForm_fQ1_vQ1(f_nodal,elementsRectangular,hx,hy); 
            bTriangular=MTriangular*f_nodal;
    end

    options = optimset('Algorithm','interior-point-convex','TolX',1e-12,'TolFun',1e-12,'TolCon',1e-12,'Display','off');
    %options = optimset('Algorithm','interior-point-convex','Display','off');
    %options = optimset('Algorithm','interior-point','Display','off');


    bTriangular_correction=-KTriangular*uDirichlet;    
    bRectangular_correction=-KRectangular*uDirichlet;  

    beq=MRectangular*phi;

    [vTriSol,dummy,dummy,dummy,multipliersTriangular]=quadprog(KTriangular(nodes_internal,nodes_internal),...
                                      -(bTriangular(nodes_internal)+bTriangular_correction(nodes_internal)),...
                                         [],[],...
                                        [],[],phi(nodes_internal),[],[],options);
                                    
    [vRecSol,dummy,dummy,dummy,multipliersRectangular]=quadprog(KRectangular(nodes_internal,nodes_internal),...
     -(bRectangular(nodes_internal)+bRectangular_correction(nodes_internal)),...
     [],[], ...
     [],[],phi(nodes_internal),[],[],options);
                                 
    vTriangular=zeros(size(coordinates,1),1);
    vRectangular=zeros(size(coordinates,1),1);
    vTriangular(nodes_dirichlet)=uDirichlet(nodes_dirichlet);
    vRectangular(nodes_dirichlet)=uDirichlet(nodes_dirichlet);
    vTriangular(nodes_internal)=vTriSol;
    vRectangular(nodes_internal)=vRecSol;
    vRectangular_finer=prolongate(vRectangular,elementsRectangular,edge2nodes);
    vRectangular_finest=prolongate(vRectangular_finer,elementsRectangular_finer,edge2nodes_finer);

    vRectangularGradient_midpoint=evaluate_GRADvQ1(vRectangular,elementsRectangular,hx,hy);

    muRectangular_nodal=zeros(size(coordinates,1),1); %muRectangular_nodal(nodes_internal)=multipliersRectangular.ineqlin;
    muRectangular_nodal(nodes_internal)=multipliersRectangular.lower;
    
    %energy distribution
    vRectangular_3Dvector=conv_ma2av(vRectangular(elementsRectangular));
    energyRectangular_distribution=(1/2)*avtamav(vRectangular_3Dvector,KRectangular_3Dmatrix,vRectangular_3Dvector) ...
                                    - avtav(bRectangular_3Dmatrix,vRectangular_3Dvector); %zde je chyba v prave strane !!!!!
    %finta :)                            
    %[bRectangular, bRectangular_3Dmatrix]=linearForm_fQ1_vQ1(f_nodal,elementsRectangular,hx,hy); %using nodal value of f
    %bRectangular2=MRectangular*f_nodal;                     

    energyRectangular=(1/2)*vRectangular'*KRectangular*vRectangular-bRectangular'*vRectangular;
    energyRectangular_finer=(1/2)*vRectangular_finer'*KRectangular_finer*vRectangular_finer-bRectangular_finer'*vRectangular_finer;
    energyRectangular_finest=(1/2)*vRectangular_finest'*KRectangular_finest*vRectangular_finest-bRectangular_finest'*vRectangular_finest;
    energyTriangular=(1/2)*vTriangular'*KTriangular*vTriangular-bTriangular'*vTriangular;

    errorNodalRectangular=u-vRectangular;
    errorNodalRectangular_finer=u_finer-vRectangular_finer;
    errorNodalRectangular_finest=u_finest-vRectangular_finest;
    errorNodalTriangular=u-vTriangular;

    u_bubbleRectangular=u_midpointRectangular-sum(u(elementsRectangular),2)/4;
    u_bubbleRectangular_finer=u_midpointRectangular_finer-sum(u_finer(elementsRectangular_finer),2)/4;
    u_bubbleRectangular_finest=u_midpointRectangular_finest-sum(u_finest(elementsRectangular_finest),2)/4;
    u_bubbleTriangular=u_midpointTriangular-sum(u(elementsTriangular),2)/3;

    if strfind(benchmark,'exact_ring')
       %errorRectangular_correction=u_H1seminorm_squared-(uDerivativeX'*MRectangular*uDerivativeX ...
       %                                                       +uDerivativeY'*MRectangular*uDerivativeY);
       errorRectangular_correction=u_H1seminorm_squared-u'*KRectangular*u;
       errorRectangular_correction_finer=u_H1seminorm_squared_finer-u_finer'*KRectangular_finer*u_finer;
       errorRectangular_correction_finest=u_H1seminorm_squared_finest-u_finest'*KRectangular_finest*u_finest;
       errorTriangular_correction= u_H1seminorm_squared-u'*KTriangular*u;

       %u_x'*MRectangular*u_x+u_y'*MRectangular*u_y

    else
       errorRectangular_correction=0;
       errorRectangular_correction_finer=0;
       errorRectangular_correction_finest=0;
       errorTriangular_correction= 0; 
    end
    
    fprintf('\n');
    fprintf('%d: u seminorm squared (using the same mesh) \n', u'*KRectangular*u);
    fprintf('%d: u seminorm squared (using the finer mesh) \n', u_finer'*KRectangular_finer*u_finer);
    fprintf('%d: u seminorm squared (using the finest mesh) \n', u_finest'*KRectangular_finest*u_finest);
    fprintf('%d: u seminorm squared (exact) \n', u_H1seminorm_squared);

    %u_bubbleRectangular=0*u_bubbleRectangular;
    errorRectangular_squared_distribution=quadratic_form_elementwise(errorNodalRectangular,elementsRectangular,KRectangular_3Dmatrix) ... + ;
                                          +(128/45)*((hx^2+hy^2)/(hx*hy))*u_bubbleRectangular.^2;                            
    errorRectangular_squared_density=errorRectangular_squared_distribution/(hx*hy);
    errorRectangular_squared=sum(errorRectangular_squared_distribution)+errorRectangular_correction;
    
    errorRectangular_squared_distribution_finer=quadratic_form_elementwise(errorNodalRectangular_finer,elementsRectangular_finer,KRectangular_finer_3Dmatrix) ... + ;
                                          +(128/45)*((hx_finer^2+hy_finer^2)/(hx_finer*hy_finer))*u_bubbleRectangular_finer.^2;                            
    errorRectangular_squared_density_finer=errorRectangular_squared_distribution_finer/(hx_finer*hy_finer);
    errorRectangular_squared_finer=sum(errorRectangular_squared_distribution_finer)+errorRectangular_correction_finer;
    
    errorRectangular_squared_distribution_finest=quadratic_form_elementwise(errorNodalRectangular_finest,elementsRectangular_finest,KRectangular_finest_3Dmatrix) ... + ;
                                          +(128/45)*((hx_finest^2+hy_finest^2)/(hx_finest*hy_finest))*u_bubbleRectangular_finest.^2;                            
    errorRectangular_squared_density_finest=errorRectangular_squared_distribution_finest/(hx_finest*hy_finest);
    errorRectangular_squared_finest=sum(errorRectangular_squared_distribution_finest)+errorRectangular_correction_finest;
    
    errorTriangular_squared=errorNodalTriangular'*KTriangular*errorNodalTriangular+errorTriangular_correction;

    energyDifferenceRectangular=energyRectangular-energy_exact;
    energyDifferenceRectangular_finer=energyRectangular_finer-energy_exact;
    energyDifferenceRectangular_finest=energyRectangular_finest-energy_exact;
    energyDifferenceTriangular=energyTriangular-energy_exact;

    Index1Rectangular=sqrt(energyDifferenceRectangular/((1/2)*errorRectangular_squared_finest));
    Index1Rectangular_finer=sqrt(energyDifferenceRectangular_finer/((1/2)*errorRectangular_squared_finest));
    Index1Rectangular_finest=sqrt(energyDifferenceRectangular_finest/((1/2)*errorRectangular_squared_finest));
    Index1Triangular=sqrt(energyDifferenceTriangular/((1/2)*errorTriangular_squared));

    fprintf('\n');
    fprintf('contact radius = %d \n', r_contact);
    fprintf('\n');
    fprintf('energy: rectangular solution = %d (using the same mesh) \n', energyRectangular);
    fprintf('energy: rectangular solution = %d (using the finer mesh) \n', energyRectangular_finer);
    fprintf('energy: rectangular solution = %d (using the finest mesh) \n', energyRectangular_finest);
    fprintf('energy: triangular  solution = %d \n', energyTriangular);
    fprintf('energy: exact       solution = %d \n', energy_exact);

    fprintf('\n');
    fprintf('(1/2) of squared error: rectangular solution = %d (bubble corrected)\n', (1/2)*errorRectangular_squared);
    fprintf('(1/2) of squared error: triangular  solution = %d \n', (1/2)*errorTriangular_squared);
    fprintf('\n');
    fprintf('energy difference: rectangular solution = %d \n', energyDifferenceRectangular);
    fprintf('energy difference:  triangular solution = %d \n', energyDifferenceTriangular);
    fprintf('\n');
    fprintf('I_eff 1 (square root of difference of energies over (1/2) of squared error): rectangular solution = %d (using the same mesh) \n',Index1Rectangular);
    fprintf('I_eff 1 (square root of difference of energies over (1/2) of squared error): rectangular solution = %d (using the finer mesh) \n',Index1Rectangular_finer);
    fprintf('I_eff 1 (square root of difference of energies over (1/2) of squared error): rectangular solution = %d (using the finest mesh) \n',Index1Rectangular_finest);
    fprintf('I_eff 1 (square root of difference of energies over (1/2) of squared error): triangular  solution = %d \n',Index1Triangular);

    %for the next paper! 
    gapRectangular=integral_constant_times_nodal_2D(lambda_midpointRectangular,vRectangular-phi,elementsRectangular,areasRectangular);
    gapTriangular=integral_constant_times_nodal_2D(lambda_midpointTriangular,vTriangular-phi,elementsTriangular,areasTriangular);

    fprintf('\n');

    visualize_solution;
              
    if strfind(benchmark,'exact_ring')
        vRectangular_circumscribed=zeros(size(coordinates_circumscribed,1),1);
        vRectangular_circumscribed(nodes2nodes_circumscribed)=vRectangular;

       f_midpointRectangular_circumscribed=f_amplitude*ones(size(elementsRectangular_circumscribed,1),1);  %works only for constant rhs!!!!!
       [dummy,dummy,elementsRectangular2edges_circumscribed]=getEdges_rectangles(elementsRectangular_circumscribed);

       [dummy,dummy,dummy,dummy,phi_nodal_circumscribed,dummy,dummy,dummy]=setup_f_and_u(benchmark,coordinates_circumscribed,f_amplitude,phi_amplitude,Lx,Ly);

       muRectangular_initial= evaluate_elements_average(elementsRectangular,muRectangular_nodal)./areasRectangular; 
       muRectangular_circumscribed_initial=zeros(size(elementsRectangular_circumscribed,1),1);
       muRectangular_circumscribed_initial(elements2elements_circumscribed)=muRectangular_initial;
       
       if bitget(bin2dec(int2str(draw)),5)          
            figure(5)   
            subplot(1,2,1)
            show_nodal_scalar(muRectangular_nodal,coordinates,elementsRectangular)
            title('nodal \mu from quadprog');
            view(-28,80);
            axis equal
            axis off
            colorbar('Southoutside');
            %screen2jpeg('squareNochetto_exactSolution')

            subplot(1,2,2);
            show_constant_scalar(muRectangular_circumscribed_initial,coordinates_circumscribed,elementsRectangular_circumscribed)
            title('nodal \mu rescalled and averaged');
            view(-28,80);
            axis equal
            axis off
            colorbar('Southoutside');
            %screen2jpeg('squareNochetto_exactSolution')   
        end                          
       
       
       [majorantRectangular,mu_midpointRectangular_circumscribed]=compute_majorant(vRectangular_circumscribed,...
                                         phi_nodal_circumscribed, ...
                                         coordinates_circumscribed,...
                                         elementsRectangular_circumscribed,...
                                         elementsRectangular2edges_circumscribed,...    
                                         muRectangular_circumscribed_initial, ...
                                         draw, ...
                                         order_f, ...
                                         f_midpointRectangular_circumscribed,...
                                         f_nodal,... 
                                         hx,hy,iterations_majorant,C_Omega,elementsRectangular_to_remove);

    else
        muRectangular_initial= evaluate_elements_average(elementsRectangular,muRectangular_nodal)./areasRectangular;   %scaling
        
        if bitget(bin2dec(int2str(draw)),5)          
            figure(5)   
            subplot(1,2,1)
            show_nodal_scalar(muRectangular_nodal,coordinates,elementsRectangular)
            title('nodal \mu from quadprog');
            view(-28,80);
            axis equal
            axis off
            colorbar('Southoutside');
            %screen2jpeg('squareNochetto_exactSolution')

            subplot(1,2,2);
            show_constant_scalar(muRectangular_initial,coordinates,elementsRectangular)
            title('nodal \mu rescalled and averaged');
            view(-28,80);
            axis equal
            axis off
            colorbar('Southoutside');
            %screen2jpeg('squareNochetto_exactSolution')

%             subplot(1,3,3)
%             show_constant_scalar(muRectangular,coordinates,elementsRectangular);
%             title('\mu majorant optimized');
%             view(-28,80);
%             axis equal
%             axis off
%             colorbar('Southoutside');       
        end                          
             
        [majorantRectangular,muRectangular]=compute_majorant(vRectangular,...
                                         phi, ...
                                         coordinates,...
                                         elementsRectangular,...
                                         elementsRectangular2edges,...
                                         muRectangular_initial, ...
                                         draw, ...
                                         order_f, ...
                                         f_averageRectangular,...
                                         f_nodal, ...
                                         hx,hy,iterations_majorant,C_Omega);      

             
    end

    Index2Rectangular=sqrt(majorantRectangular/energyDifferenceRectangular);
    Index2Rectangular_finer=sqrt(majorantRectangular/energyDifferenceRectangular_finer);
    Index2Rectangular_finest=sqrt(majorantRectangular/energyDifferenceRectangular_finest);

    fprintf('\n'); 
    fprintf('I_eff 2 (square root of majorant over difference of energies): rectangular solution = %d (using the same mesh) \n',Index2Rectangular);
    fprintf('I_eff 2 (square root of majorant over difference of energies): rectangular solution = %d (using the finer mesh) \n',Index2Rectangular_finer);
    fprintf('I_eff 2 (square root of majorant over difference of energies): rectangular solution = %d (using the finest mesh) \n',Index2Rectangular_finest);
    fprintf('------------ \n');

    errorRectangular_squared_all(level)=errorRectangular_squared;
    errorRectangular_squared_finer_all(level)=errorRectangular_squared_finer;
    errorRectangular_squared_finest_all(level)=errorRectangular_squared_finest;
    energyDifferenceRectangular_all(level)=energyDifferenceRectangular;
    energyDifferenceRectangular_finer_all(level)=energyDifferenceRectangular_finer;
    energyDifferenceRectangular_finest_all(level)=energyDifferenceRectangular_finest;
    majorantRectangular_all(level)=majorantRectangular;
    Index1Rectangular_all(level)=Index1Rectangular;
    Index1Rectangular_finer_all(level)=Index1Rectangular_finer;
    Index1Rectangular_finest_all(level)=Index1Rectangular_finest;
    Index2Rectangular_all(level)=Index2Rectangular;
    problemSize_all(level)=size(coordinates,1);    
end

fprintf('\n')
for l=min_level_refinement:level
        
        fprintf('level %d, ',l)
        fprintf('1/2 of squared error %d, ',errorRectangular_squared_finest_all(l)/2)
        fprintf('difference of energies %d, ',energyDifferenceRectangular_all(l))
        fprintf('majorant %d, ',majorantRectangular_all(l))
        fprintf('nodes %d ',problemSize_all(l))
        fprintf('\n')
end

if bitget(bin2dec(int2str(draw)),8)  
    figure(8)
    loglog(problemSize_all,majorantRectangular_all,'x-', ...
           problemSize_all,energyDifferenceRectangular_all,'x-.',...
           problemSize_all,errorRectangular_squared_finest_all/2,'x:','LineWidth',2,'MarkerSize',12);
    legend('majorant','diffence of energies','1/2 of squared error')
    axis tight
    xlabel('problem size = number of mesh nodes')
    %screen2jpeg('squareNochetto_convergence')
end



end

