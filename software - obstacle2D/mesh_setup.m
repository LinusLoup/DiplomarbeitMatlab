nex=(2^level); 
ney=(2^level); 
nx=nex+1; % number of nodes in x direction (minimum 2!)
ny=ney+1; % number of nodes in y direction (minimum 2!)

%mesh generation
hx=Lx/(nx-1); %mesh size in x direction
hy=Ly/(ny-1); %mesh size in x direction

hx_finer=hx/2;
hy_finer=hy/2;

hx_finest=hx_finer/2;
hy_finest=hy_finer/2;


%coordinates 
X=(kron(0:hx:Lx,ones(ny,1)))';      %x coordinates of nodes
Y=(kron(ones(1,nx),(0:hy:Ly)'))';   %y coordinates of nodex
%coordinates=[X(:) Y(:)];

factor=1;
coordinates=[(X(:)-Lx/2)*factor+Lx/2 (Y(:)-Ly/2)*factor+Ly/2];
hx=hx*factor;
hy=hy*factor;


%centres coordinates 
%XC=(kron((hx/2):hx:(Lx-(hx/2)),ones(ney,1)))';      
%YC=(kron(ones(1,nex),((hy/2):hy:(Ly-(hy/2)))'))';   
%elementsRectangular2midpoint=[XC(:) YC(:)];

%elements rectangular
element_ref=[1 2 nx+2 nx+1];
element_x=[];
for i=1:(nx-1)
    element_x=[element_x; element_ref+i-1];
end
elementsRectangular=[];
for j=1:(ny-1)
    elementsRectangular=[elementsRectangular; element_x+(nx)*(j-1)];
end
clear element_x element_ref

% modify rectangular mesh to a ring mesh

if strfind(benchmark,'exact_ring')
    %large mesh than ring
    distance_to_center=sqrt(sum([coordinates(:,1)-Lx/2 coordinates(:,2)-Ly/2].^2,2)); 
    elementsRectangular_to_remove=find(sum(distance_to_center(elementsRectangular)>1,2)>=4); 
    [elementsRectangular_circumscribed,coordinates_circumscribed]=remove_elements_and_coordinates(elementsRectangular,coordinates,elementsRectangular_to_remove);
       
    %smaller mesh than ring
    distance_to_center=sqrt(sum([coordinates_circumscribed(:,1)-Lx/2 coordinates_circumscribed(:,2)-Ly/2].^2,2)); 
    elementsRectangular_to_remove=find(sum(distance_to_center(elementsRectangular_circumscribed)>1,2)>=1); 
    [elementsRectangular,coordinates,nodes2nodes_circumscribed,elements2elements_circumscribed,dummy]=remove_elements_and_coordinates(elementsRectangular_circumscribed,coordinates_circumscribed,elementsRectangular_to_remove);
end

nn=size(coordinates,1);
neRectangular=size(elementsRectangular,1);

%underlining triangular mesh
elementsTriangular1=elementsRectangular(:,1:3);
elementsTriangular2=elementsRectangular(:,[3 4 1]);
elementsTriangular=zeros(2*neRectangular,3);
elementsTriangular(1:2:2*neRectangular,:)=elementsTriangular1;
elementsTriangular(2:2:2*neRectangular,:)=elementsTriangular2;
clear elementsTriangular1 elementsTriangular2

neTriangular=2*neRectangular;

[coordinates_finer,elementsRectangular_finer]=refinement_uniform_2D(coordinates,elementsRectangular,[],[],[]);  
[coordinates_finest,elementsRectangular_finest]=refinement_uniform_2D(coordinates_finer,elementsRectangular_finer,[],[],[]);
    
elementsRectangular2midpoint=evaluate_elements_average(elementsRectangular,coordinates);
elementsRectangular2midpoint_finer=evaluate_elements_average(elementsRectangular_finer,coordinates_finer);
elementsRectangular2midpoint_finest=evaluate_elements_average(elementsRectangular_finest,coordinates_finest);
elementsTriangular2midpoint=evaluate_elements_average(elementsTriangular,coordinates);

%edges in elements ordered as: bottom, right, top, left, this is important!!!!!
[edge2nodes, edge2elementsRectangular,elementsRectangular2edges]=getEdges_rectangles(elementsRectangular);
[edge2nodes_finer, edge2elementsRectangular_finer,elementsRectangular2edges_finer]=getEdges_rectangles(elementsRectangular_finer);
%edges_internal=find(edge2elementsRectangular(:,2));   %internal edges  
edges_boundary=(find(~edge2elementsRectangular(:,2))); %boundary edges
nodes_dirichlet=unique(edge2nodes(edges_boundary,:));
nodes_internal=setdiff((1:nn)',nodes_dirichlet); 





