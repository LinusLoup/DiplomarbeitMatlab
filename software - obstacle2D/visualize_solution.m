if bitget(bin2dec(int2str(draw)),1)  
    %exact solution and multiplier
    figure(1)  
    subplot(1,4,1)
    show_nodal_scalar(u,coordinates,elementsRectangular);
    hold on
    show_nodal_scalar(phi,coordinates,elementsRectangular)
    title('u');
    set_figure; view(-24,54);
    colorbar('Southoutside');
    hold off

    subplot(1,4,2)
    %trisurf(elementsRectangular,coordinates(:,1),coordinates(:,2),lambda);
    show_constant_scalar(lambda_midpointRectangular,coordinates,elementsRectangular);
    title('\lambda');
    set_figure
    colorbar('Southoutside');

    %exact flux
    subplot(1,4,3)
    show_nodal_scalar(uGradientX,coordinates,elementsRectangular)
    title('\nabla_x u');
    set_figure
    colorbar('Southoutside');

    subplot(1,4,4)
    show_nodal_scalar(uGradientY,coordinates,elementsRectangular)
    title('\nabla_y u');
    set_figure
    colorbar('Southoutside');
end

%discrete rectangular solution
if bitget(bin2dec(int2str(draw)),2)      
    figure(2) 
    subplot(1,3,1)
    show_nodal_scalar(vRectangular,coordinates,elementsRectangular)
    hold on
    show_nodal_scalar(phi,coordinates,elementsRectangular)
    title('v');
    set_figure; view(-24,54);
    colorbar('Southoutside');
    hold off

    %colorbar('Southoutside');

    subplot(1,3,2)
    show_constant_scalar(vRectangularGradient_midpoint(:,1),coordinates,elementsRectangular);
    title('\nabla_x v');
    set_figure
    colorbar('Southoutside');

    subplot(1,3,3)
    show_constant_scalar(vRectangularGradient_midpoint(:,2),coordinates,elementsRectangular);
    title('\nabla_y v');
    %view(2);
    set_figure
    colorbar('Southoutside');

    %discrete triangular solution
    %figure(3)     
    %show_nodal_scalar(vTriangular,coordinates,elementsTriangular)
    %title('triangular solution v');
    %set_figure
    %colorbar('Southoutside');
end

%loading 
if bitget(bin2dec(int2str(draw)),3)  
    figure(3)
    subplot(1,2,1)
    show_nodal_scalar(f_nodal,coordinates,elementsRectangular)
    title('rhs f');
    set_figure; view(-24,82);
    colorbar('Southoutside');

    subplot(1,2,2)
    show_constant_scalar(f_midpointRectangular,coordinates,elementsRectangular);
    title('rhs f midpoint value');
    set_figure; view(-24,82);
    colorbar('Southoutside');

%     subplot(1,3,3)
%     show_constant_scalar(f_averageRectangular,coordinates,elementsRectangular);
%     title('rhs f average value');
%     set_figure; view(-24,82);
%     colorbar('Southoutside');
end

%exact error distribution
if bitget(bin2dec(int2str(draw)),4)  
      figure(4)
    subplot(1,3,1)
    show_constant_scalar(errorRectangular_squared_density/2,coordinates,elementsRectangular);
    title(strcat('|| \nabla (u_h - v) ||^2 /2 =',num2str(errorRectangular_squared/2)),'FontSize', 20);
    set_figure
    colorbar('Southoutside');
    
    subplot(1,3,2)
    show_constant_scalar(errorRectangular_squared_density_finer/2,coordinates_finer,elementsRectangular_finer);
    title(strcat('|| \nabla (u_{h/2} - v) ||^2 /2 =',num2str(errorRectangular_squared_finer/2)),'FontSize', 20);
    set_figure
    colorbar('Southoutside');
    
    subplot(1,3,3)
    show_constant_scalar(errorRectangular_squared_density_finest/2,coordinates_finest,elementsRectangular_finest);
    title(strcat('|| \nabla (u_{h/4} - v) ||^2 /2 =',num2str(errorRectangular_squared_finest/2)),'FontSize', 20);
    set_figure
    colorbar('Southoutside'); 

%       show_constant_scalar(errorRectangular_squared_density/2,coordinates,elementsRectangular);
%       title(strcat('|| \nabla (u - v) ||^2 /2 =',num2str(errorRectangular_squared_finest/2)),'FontSize', 20);
%       set_figure
%       view(2)
%       colorbar('Southoutside')    
end

if bitget(bin2dec(int2str(draw)),9)  
    figure(101)
    show_nodal_scalar(u,coordinates,elementsRectangular);
    hold on
    show_nodal_scalar(phi,coordinates,elementsRectangular)
    %title('exact solution u');
    set_figure
    view(19,36)
    colorbar
    hold off
    
    figure(102)
    show_nodal_scalar(uGradientX,coordinates,elementsRectangular)
    %title('exact flux \sigma y-component');
    set_figure
    view(-24,78)
    colorbar
    %axis on
    %colorbar('Southoutside');
    %screen2jpeg('squareNochetto_exactGradientX')
    
    figure(103)
    show_constant_scalar(lambda_midpointRectangular,coordinates,elementsRectangular);
    set_figure
    view(-24,78)
    colorbar
    %axis on
    %colorbar('Southoutside');
    %screen2jpeg('squareNochetto_exactMultiplier')
    
    figure(104)
    show_nodal_scalar(vRectangular,coordinates,elementsRectangular)
    hold on
    show_nodal_scalar(phi,coordinates,elementsRectangular)
    set_figure
    view(19,36)
    colorbar
    %axis on
    hold off
    %screen2jpeg('squareNochetto_exactSolution')
    
    
end






