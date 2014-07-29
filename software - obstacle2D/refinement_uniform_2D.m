function [coordinates_finer,elements_finer,dirichlet_finer,neumann_finer,element_parameters_finer]=refinement_uniform_2D(coordinates,elements,dirichlet,neumann,element_parameters) 
%function: [coordinates,elements3,dirichlet]=refinement_uniform(coordinates,elements3,dirichlet) 
%requires: getEdges, symetrizeMatrix  
%uniform refinement of a 2D triangulation 

%uniform refinement   
[elements2edges, edge2nodes, dummy, nodes2edge]=getEdges(elements);    

%dirichlet edges of uniformly refined mesh
if ~isempty(dirichlet)
    dirichlet_edges=diag(nodes2edge(dirichlet(:,1),dirichlet(:,2)));
    if (size(dirichlet,2)==2)
        dirichlet_finer=[dirichlet(:,1) dirichlet_edges+size(coordinates,1); dirichlet_edges+size(coordinates,1) dirichlet(:,2)];     
    else  %3d not working yet, prepared
        dirichlet_finer=[dirichlet(:,1) dirichlet_edges+size(coordinates,1) dirichlet(:,3); ...
                   dirichlet_edges+size(coordinates,1) dirichlet(:,2) dirichlet(:,3)];
    end
end

%neumann edges of uniformly refined mesh
if ~isempty(neumann)
    neumann_edges=diag(nodes2edge(neumann(:,1),neumann(:,2)));
    neumann_finer=[neumann(:,1) neumann_edges+size(coordinates,1); neumann_edges+size(coordinates,1) neumann(:,2)];
end

%elements on uniformly refined mesh
if (size(coordinates,2)==2)  
        elements2edgesMidnodes=elements2edges+size(coordinates,1);
        
        %triangles
        if (size(elements,2)==3)       
            elements_refin1= [elements(:,1) elements2edgesMidnodes(:,3) elements2edgesMidnodes(:,2)];
            elements_refin2= [elements(:,2) elements2edgesMidnodes(:,1) elements2edgesMidnodes(:,3)];
            elements_refin3= [elements(:,3) elements2edgesMidnodes(:,2) elements2edgesMidnodes(:,1)];    
            elements_finer=[elements2midnodes; elements_refin1; elements_refin2; elements_refin3];  
        end
        
        %rectangles
        if (size(elements,2)==4)
            elements2Midnode=(1+size(coordinates,1)+size(edge2nodes,1):size(elements,1)+size(coordinates,1)+size(edge2nodes,1))';
            elements_refin1= [elements(:,1) elements2edgesMidnodes(:,1) elements2Midnode elements2edgesMidnodes(:,4)];
            elements_refin2= [elements2edgesMidnodes(:,1) elements(:,2) elements2edgesMidnodes(:,2) elements2Midnode];
            elements_refin3= [elements2Midnode elements2edgesMidnodes(:,2) elements(:,3) elements2edgesMidnodes(:,3)];  
            elements_refin4= [elements2edgesMidnodes(:,4) elements2Midnode elements2edgesMidnodes(:,3) elements(:,4)];  
            elements_finer=[elements_refin1; elements_refin2; elements_refin3; elements_refin4];
        end     
end

if (nargout==5) && (nargin==5)
    element_parameters_finer=kron(ones(4,1),element_parameters);
end

%coordinates of uniformly refined mesh
edges2midpoint=(coordinates(edge2nodes(:,1),:)+coordinates(edge2nodes(:,2),:))/2;
if (size(coordinates,2)==2)  
        %triangles
        if (size(elements,2)==3)   
            coordinates_finer=[coordinates; edges2midpoint];  
        end
        
        %rectangles
        if (size(elements,2)==4)   
            elements2midpoint=evaluate_elements_average(elements,coordinates);
            coordinates_finer=[coordinates; edges2midpoint; elements2midpoint];  
        end
end



        

 

    function [element2edges, edge2nodes, edge2elements,nodes2edge]=getEdges(elements)
    %function: [element2edges, edge2nodes]=edge_numbering(elements)
    %requires: deleterepeatedrows
    %generates edges of (triangular) triangulation defined in elements
    %elements is matrix, whose rows contain numbers of its element nodes 
    %element2edges returns edges numbers of each triangular element
    %edge2nodes returns two node numbers of each edge
    %example in 2D: [element2edges, edge2nodes]=getEdges([1 2 3; 2 4 3])
    %example in 3D: [element2edges, edge2nodes]=getEdges([1 2 3 4; 1 2 3 5; 1 2 4 6])

    %2D case
    if (size(coordinates,2)==2)  
        %triangles
        if (size(elements,2)==3)
            %extracts sets of edges 
            edges1=elements(:,[2 3]);
            edges2=elements(:,[3 1]);
            edges3=elements(:,[1 2]);

            %as sets of their nodes (vertices)
            vertices=zeros(size(elements,1)*3,2);
            vertices(1:3:end,:)=edges1;
            vertices(2:3:end,:)=edges2;
            vertices(3:3:end,:)=edges3;

            %repeated sets of nodes (joint edges) are eliminated 
            [edge2nodes,element2edges]=deleterepeatedrows(vertices);
            element2edges=reshape(element2edges,3,size(elements,1))';
        end
        
        %rectangles
        if (size(elements,2)==4)
            %extracts sets of edges 
            edges1=elements(:,[1 2]);
            edges2=elements(:,[2 3]);
            edges3=elements(:,[3 4]);
            edges4=elements(:,[4 1]);

            %as sets of their nodes (vertices)
            vertices=zeros(size(elements,1)*4,2);
            vertices(1:4:end,:)=edges1;
            vertices(2:4:end,:)=edges2;
            vertices(3:4:end,:)=edges3;
            vertices(4:4:end,:)=edges4;
            
            %repeated sets of nodes (joint edges) are eliminated 
            [edge2nodes,element2edges]=deleterepeatedrows(vertices);
            element2edges=reshape(element2edges,4,size(elements,1))';
        end
        
    end
    
    %3D case
    if (size(coordinates,2)==3)  
        if (size(elements,2)==4)
            %extracts sets of edges 
            edges1=elements(:,[1 2]);
            edges2=elements(:,[1 3]);
            edges3=elements(:,[1 4]);
            edges4=elements(:,[2 3]);
            edges5=elements(:,[2 4]);
            edges6=elements(:,[3 4]);

            %as sets of their nodes (vertices)
            vertices=zeros(size(elements,1)*6,2);
            vertices(1:6:end,:)=edges1;
            vertices(2:6:end,:)=edges2;
            vertices(3:6:end,:)=edges3;
            vertices(4:6:end,:)=edges4;
            vertices(5:6:end,:)=edges5;
            vertices(6:6:end,:)=edges6;

            %repeated sets of nodes (joint edges) are eliminated 
            [edge2nodes,element2edges]=deleterepeatedrows(vertices);
            element2edges=reshape(element2edges,6,size(elements,1))';
        end
    end
    
    edge2elements=entryInWhichRows(element2edges); 
    if (size(edge2elements,2)==1)
        edge2elements=[edge2elements 0*edge2elements];    %all edges are boundary edges!!!
    end
    
    nodes2edge=sparse(edge2nodes(:,1),edge2nodes(:,2),1:size(edge2nodes,1),size(coordinates,1),size(coordinates,1));
    nodes2edge=symetrizeMatrix(nodes2edge); 
    
    end

    function A_sym = symetrizeMatrix(A)
        [i,j,k]=find(A);
        W=sparse([i; j], [j; i], ones(size(k,1)*2,1));
        A_help=sparse([i; j], [j; i], [k; k]);
        [i,j,k]=find(A_help);
        [i,j,kk]=find(W);
        A_sym=sparse(i,j,(kk.^(-1)).*k); %Now Kantennr_sym is a symetric form of Kantennr
    end
   
    function [matrix,I]=deleterepeatedrows(matrix)
        %function: [element2edges, edge2nodes]=edge_numbering(elements)
        %requires: deleterepeatedrows
        %generates edges of (triangular) triangulation defined in elements
        %elements is matrix, whose rows contain numbers of its element nodes 
        %element2edges returns edges numbers of each triangular element
        %edge2nodes returns two node numbers of each edge
        %example: [element2edges, edge2nodes]=edge_numbering([1 2 3; 2 4 3])


        %fast and short way suggested by John D'Ericco working in both 2D and 3D
        matrixs=sort(matrix,2);
        [dummy,J,I] = unique(matrixs,'rows');
        %I=reshape(I,size(matrixs,2),size(I,1)/size(matrixs,2));
        matrix=matrix(J,:);
    end

    function entryrows=entryInWhichRows(A)
        %function: entryrows=entryInWhichRows(A)
        %requires: none
        %for every entry of integer matrix A, 
        %its rows indices are stored in output matrix,
        %zeros entries indicate no more occurence
        %example: entryrows=entryInWhichRows([1 2; 1 3; 2 2]) returns
        %         entryrows=[1   2   0; 
        %                    1   3   3; 
        %                    2   0   0]   
        %meaning: entry 1 appears in rows 1 and 2
        %         entry 2 appears in rows 1 and 3 (twice)
        %         entry 3 appears in row  2 only

        %size computation; 
        r=max(max(A));
        repetition=accumarray(A(:),ones(numel(A),1));
        c=max(repetition);

        %filling rows occurences
        %this part should be somehow vectorized!
        entryrows=zeros(r,c);
        repetition=zeros(r,1);
        for i=1:size(A,1)
            for j=1:size(A,2)
               index=A(i,j); 
               repetition(index)=repetition(index)+1; 
               entryrows(index,repetition(index))=i;
            end
        end
    end

end


