function [majorantRectangular,muRectangular]=...
              compute_majorant(vRectangular, phi_nodal, ...
                                              coordinates, ...
                                              elementsRectangular, ...
                                              elementsRectangular2edges,...
                                              muRectangular_initial, ...
                                              draw, ...
                                              order_f, ...
                                              f_averageRectangular, ...
                                              f_nodal, ...     
                                              hx,hy,iterations_majorant,C_Omega,elements_to_remove)                                        
        
iteration_display=ceil(iterations_majorant/5);
%iteration_display=1;
iteration_betaupdate=iterations_majorant;
                                          
beta=1; 
muRectangular=muRectangular_initial;

[KRectangular, MRectangular, KRectangular_3Dmatrix, MRectangular_3Dmatrix]=assembly_FEM_matrices_rectangle_bilinear(elementsRectangular,hx,hy); %redundant, used for computation of majorant_left
vRectangular_H1seminorm_squared=vRectangular'*KRectangular*vRectangular;
vRectangular_H1seminorm_squared_distribution=quadratic_form_elementwise(vRectangular,elementsRectangular,KRectangular_3Dmatrix);
                                   
[K_RT0, M_RT0, dummy, M_RT0_3Dmatrix]=assembly_FEM_matrices_rectangle_RT0(elementsRectangular2edges,hx,hy);
     
vRectangular_gradient_times_tauRectangularRT0=vRectangular(elementsRectangular)*(1/4)* ...
            [-hx -hy -hx -hy; ...
             -hx  hy -hx  hy; ...
              hx  hy  hx  hy; ...
              hx -hy  hx -hy];

b_RT0=sparse(elementsRectangular2edges,ones(size(elementsRectangular2edges)),vRectangular_gradient_times_tauRectangularRT0);  

for iteration=1:iterations_majorant
    
    %tau computation
    c_RT0=sparse(elementsRectangular2edges,ones(size(elementsRectangular2edges)),kron(f_averageRectangular+muRectangular,[-hx hy hx -hy]));      %there might be a bug!
    tauRectangularRT0=((1+beta)*M_RT0+(1+1/beta)*C_Omega^2*K_RT0) \  ((1+beta)*b_RT0-(1+1/beta)*C_Omega^2*c_RT0); 
    tauRectangularRT0=full(tauRectangularRT0);
    tauRectangular_div=tauRectangularRT0(elementsRectangular2edges)*[-1/hy; 1/hx; 1/hy; -1/hx];  %divergence of elements
    
    if  ~mod(iteration,iteration_display) 
        [majorantRectangular, majorantRectangular_density,...
                  majorantRectangular_left,majorantRectangular_right,majorantRectangular_nonlinear,...
                  majorantRectangular_left_density,majorantRectangular_right_density,majorantRectangular_nonlinear_density]=evaluate_majorant_value(1);
    end
    
    %multiplicator computation 
    muRectangular=-(f_averageRectangular+tauRectangular_div)-(sum(vRectangular(elementsRectangular)-phi_nodal(elementsRectangular),2)/4)/((1+1/beta)*(C_Omega^2));
    muRectangular=max(muRectangular,zeros(size(muRectangular))); 
    
    if  ~mod(iteration,iteration_display) 
        [majorantRectangular, majorantRectangular_density, ...
                  majorantRectangular_left,majorantRectangular_right,majorantRectangular_nonlinear,...
                  majorantRectangular_left_density,majorantRectangular_right_density,majorantRectangular_nonlinear_density]=evaluate_majorant_value(1);
    end
      
    %switching beta
    if  ~mod(iteration,iteration_betaupdate) && (iteration<=iterations_majorant)
          %[majorantRectangular, majorantRectangular_density, ...
          %        majorantRectangular_left,majorantRectangular_right,majorantRectangular_nonlinear,...
          %        majorantRectangular_left_density,majorantRectangular_right_density,majorantRectangular_nonlinear_density]=evaluate_majorant_value(0);     
        
        beta=sqrt(C_Omega^2*majorantRectangular_right/majorantRectangular_left);
        
          [majorantRectangular, majorantRectangular_density, ...
                  majorantRectangular_left,majorantRectangular_right,majorantRectangular_nonlinear,...
                  majorantRectangular_left_density,majorantRectangular_right_density,majorantRectangular_nonlinear_density]=evaluate_majorant_value(1);
    end  
   
    % if  ~mod(iteration,iteration_display) 
         
     %   visualize_majorant(tauRectangularRT0,coordinates,elementsRectangular,elementsRectangular2edges);
    
     %end
     
        
end

 %[majorantRectangular, majorantRectangular_density, ...
 %                 majorantRectangular_left,majorantRectangular_right,majorantRectangular_nonlinear,...
 %                 majorantRectangular_left_density,majorantRectangular_right_density,majorantRectangular_nonlinear_density]=evaluate_majorant_value(1);

visualize_majorant(tauRectangularRT0,coordinates,elementsRectangular,elementsRectangular2edges);

 if bitget(bin2dec(int2str(draw)),9)           
            figure(105)
            fill3(X',Y',tau_RT0_element2nodevaluesX',tau_RT0_element2nodevaluesX','FaceColor','interp','LineStyle','none');
            set_figure
            view(-24,78)
            colorbar
            % axis on    
            %screen2jpeg('squareNochetto_discreteGradientX')       
            
            figure(106)
            show_constant_scalar(muRectangular,coordinates,elementsRectangular);
            set_figure
            view(-24,78)
            colorbar
            %axis on
            
          end






    function [majorantRectangular, majorantRectangular_density,...
              majorantRectangular_left,majorantRectangular_right,majorantRectangular_nonlinear,...
              majorantRectangular_left_density,majorantRectangular_right_density,majorantRectangular_nonlinear_density]=evaluate_majorant_value(visualize)
        
        tau_RT0_L2norm_squared_distribution=quadratic_form_elementwise(tauRectangularRT0,elementsRectangular2edges,M_RT0_3Dmatrix);


        majorantRectangular_left_distribution=vRectangular_H1seminorm_squared_distribution ...
                                              + tau_RT0_L2norm_squared_distribution  ...
                                              -2*sum(vRectangular_gradient_times_tauRectangularRT0.*tauRectangularRT0(elementsRectangular2edges),2);
        majorantRectangular_left_density=majorantRectangular_left_distribution/(hx*hy);
        majorantRectangular_left=(vRectangular_H1seminorm_squared+tauRectangularRT0'*M_RT0*tauRectangularRT0-2*b_RT0'*tauRectangularRT0);

        switch order_f
            case 0
                majorantRectangular_right_density=(f_averageRectangular+muRectangular+tauRectangular_div).^2;  
            case {1,2}
                [dummy, f_plus_mu_L2normSquared_density]=L2norm_fQ1_plus_muQ0_squared(f_nodal(elementsRectangular),muRectangular,elementsRectangular,hx,hy);
                %majorantRectangular_right2=(sum(hx*hy*(f_averageRectangular+mu).^2) ...
                %                  +2*c_RT0'*tauRectangularRT0 ...
                %                  +tauRectangularRT0'*K_RT0*tauRectangularRT0); %probably something wrong with c_RT0 term
                majorantRectangular_right_density=f_plus_mu_L2normSquared_density...
                                                  +tauRectangular_div.^2 ...
                                                  +2*(f_averageRectangular+muRectangular).*tauRectangular_div;
        
        end     
        majorantRectangular_right=sum(hx*hy*majorantRectangular_right_density);

        %pause

        majorantRectangular_nonlinear_density=(sum(vRectangular(elementsRectangular)-phi_nodal(elementsRectangular),2)/4).*muRectangular;
        majorantRectangular_nonlinear=hx*hy*sum(majorantRectangular_nonlinear_density);

        
        majorantRectangular_density=(1/2)*( (1+beta)* majorantRectangular_left_density ...
                                    +(1+1/beta)*C_Omega^2*majorantRectangular_right_density) ...
                                    +majorantRectangular_nonlinear_density;
        majorantRectangular=hx*hy*sum(majorantRectangular_density);
                                
                                
                                
        %display of iterations
        if visualize
            fprintf('majorant: rectangular solution =%10.3e, ',majorantRectangular);
            fprintf('left =%10.3e, ',majorantRectangular_left);
            fprintf('right =%10.3e, ',C_Omega^2*majorantRectangular_right);
            fprintf('nonlinear =%10.3e, ',majorantRectangular_nonlinear);
            fprintf('beta=%10.3e, ', beta);
            fprintf('iteration: %d \n', iteration);              
        end
        
    end




    function visualize_majorant(tau_RT0,coordinates,elementsRectangular,elementsRectangular2edges)
        
        coordinates_x=coordinates(:,1)-1;
        coordinates_y=coordinates(:,2)-1;
        X=coordinates_x(elementsRectangular);
        Y=coordinates_y(elementsRectangular);      
        tau_RT0_element2nodevaluesX=tau_RT0(elementsRectangular2edges)*[0 0 0 0; 0 1 1 0; 0 0 0 0; 1 0 0 1];
        tau_RT0_element2nodevaluesY=tau_RT0(elementsRectangular2edges)*[1 1 0 0; 0 0 0 0; 0 0 1 1; 0 0 0 0];
      
        if bitget(bin2dec(int2str(draw)),6)  
            figure(6)
            subplot(2,2,1)
            show_constant_scalar(majorantRectangular_density,coordinates,elementsRectangular);
            %fill3(X',Y',(1/max(max(majorantRectangular_left_density)))*kron(ones(1,4),majorantRectangular_density)',kron(ones(1,4),majorantRectangular_density)');
            %title(strcat('|| \nabla v - \tau ||^2 =',num2str(majorantRectangular_left)),'FontSize', 20);        
            set_figure
            view(2)
            colorbar('Southoutside');
            title(strcat('majorant =',num2str(majorantRectangular)),'FontSize', 20);
            %screen2jpeg(strcat('majorant_total_iteration',32,num2str(iteration)))      
            
            subplot(2,2,2)
            show_constant_scalar(majorantRectangular_left_density,coordinates,elementsRectangular);
            title(strcat('|| \nabla v - \tau ||^2 =',num2str(majorantRectangular_left)),'FontSize', 20);        
            set_figure
            view(2)
            colorbar('Southoutside');

            subplot(2,2,3)
            show_constant_scalar(majorantRectangular_right_density,coordinates,elementsRectangular);
            title(strcat('|| div \tau + f + \mu||^2 =',num2str(majorantRectangular_right)),'FontSize', 20);
            set_figure
             view(2)
            colorbar('Southoutside');

            subplot(2,2,4)
            show_constant_scalar(majorantRectangular_nonlinear_density,coordinates,elementsRectangular);
            title(strcat('\int_{\Omega} \mu (v-\phi) d\Omega =',num2str(majorantRectangular_nonlinear)),'FontSize', 20);
            set_figure
             view(2)
            colorbar('Southoutside');  
            
             %screen2jpeg(strcat('majorant_distributions_iteration',32,num2str(iteration)))     
        end
        
        
        if bitget(bin2dec(int2str(draw)),7)  
            figure(7) 
            subplot(1,3,1);    
            fill3(X',Y',tau_RT0_element2nodevaluesX',tau_RT0_element2nodevaluesX','FaceColor','interp','LineStyle','none');
            title('\tau^*_x');
            set_figure
            colorbar('Southoutside');

            subplot(1,3,2);
            fill3(X',Y',tau_RT0_element2nodevaluesY',tau_RT0_element2nodevaluesY','FaceColor','interp','LineStyle','none');
            title('\tau^*_y');
            set_figure
            colorbar('Southoutside');

            subplot(1,3,3);
            show_constant_scalar(muRectangular,coordinates,elementsRectangular);
            title('\mu');
            set_figure
            colorbar('Southoutside');      
        end
            
        
            
            if exist('elements_to_remove','var')            
                X(elements_to_remove,:)=[];
                Y(elements_to_remove,:)=[];
                tau_RT0_element2nodevaluesX(elements_to_remove,:)=[];
                tau_RT0_element2nodevaluesY(elements_to_remove,:)=[];
                %muRectangular(elements_to_remove)=[];  %attention!!

%                 figure(888) 
%                 subplot(1,2,1);    
%                 fill3(X',Y',tau_RT0_element2nodevaluesX',tau_RT0_element2nodevaluesX');
%                 title('\tau x-component');
%                 set_figure
%                 colorbar('Southoutside');
% 
%                 subplot(1,2,2);
%                 fill3(X',Y',tau_RT0_element2nodevaluesY',tau_RT0_element2nodevaluesY');
%                 title('\tau y-component');
%                 set_figure
%                 colorbar('Southoutside');

%                 subplot(1,3,3);
%                 fill3(X',Y',kron(ones(1,4),muRectangular)',kron(ones(1,4),muRectangular)');
%                 title('\mu');
%                 set_figure
%                 colorbar('Southoutside');      
            end        
        end    
end

%v_test=v_test_function(coordinates);
%v_test_compX=v_test(:,1);
%v_test_compY=v_test(:,2);

% B_nodal_RT0CompX=(1/12)*[0 0 0 0; ...
%                          hx*hy 2*hx*hy 2*hx*hy hx*hy; ...
%                          0 0 0 0; ...
%                          2*hx*hy hx*hy hx*hy 2*hx*hy];
% 
% B_nodal_RT0CompY=(1/12)*[2*hx*hy 2*hx*hy hx*hy hx*hy; ...
%                          0 0 0 0; ...
%                          hx*hy hx*hy 2*hx*hy 2*hx*hy; ...
%                          0 0 0 0];

%Z=v_test_compX(elementsRectangular)*B_nodal_RT0CompX+v_test_compY(elementsRectangular)*B_nodal_RT0CompY;
%sparse(elementsRectangular2edges,ones(size(elementsRectangular2edges)),Z,size(edge2nodes,1),1);     

%v_test_RT=full(M_RT0\sparse(elementsRectangular2edges,ones(size(elementsRectangular2edges)),Z,size(edge2nodes,1),1)); 

% gradX_in_IP=(1/hx)*[-1 -1/2  0 -1/2 -1/2; ...
%                      1  1/2  0  1/2  1/2; ...
%                      0  1/2  1  1/2  1/2; ...
%                      0 -1/2 -1 -1/2 -1/2];
%           
% gradY_in_IP=(1/hy)*[-1/2  0  -1/2  -1  -1/2; ...
%                     -1/2 -1  -1/2   0  -1/2; ...
%                      1/2  1   1/2   0   1/2; ...
%                      1/2  0   1/2   1   1/2];
%          
% valueX_in_IP=[  0   0  0   0  0;   ...
%               1/2   1  1/2 0  1/2; ...
%                 0   0  0   0  0;   ...
%               1/2  0  1/2 1  1/2];
%           
% valueY_in_IP=[ 1   1/2   0   1/2  1/2; ...
%                0     0   0    0     0; ...
%                0    1/2  1   1/2  1/2;   ...
%                0     0   0    0    0];
%            
% bilinear_form=(hx*hx/6)*(gradX_in_IP*diag([1 1 1 1 2])*valueX_in_IP'+gradY_in_IP*diag([1 1 1 1 2])*valueY_in_IP');
% 
% 
% b_RT0=sparse(elementsRectangular2edges,ones(size(elementsRectangular2edges)),vRectangular(elementsRectangular)*bilinear_form);
