function [elements,coordinates,nodes2nodes_old,elements2elements_old, nodes_to_remove]=remove_elements_and_coordinates(elements,coordinates,elements_to_remove)

        if numel(elements_to_remove)>0
            elements2elements_old=setdiff((1:size(elements,1))',elements_to_remove);
                 

            elements(elements_to_remove,:)=[];

            nodes_to_remove=setdiff( (1:size(coordinates,1))',unique(elements));
            
            nodes2nodes_old=setdiff((1:size(coordinates,1))',nodes_to_remove);
            

            shift_current=0;
            shift=zeros(size(coordinates,1),1);
            for i=1:size(nodes_to_remove,1)-1
                shift_current=shift_current+1;
                for j=nodes_to_remove(i):nodes_to_remove(i+1)-1
                    shift(j)=shift_current;
                end
            end
            shift_current=shift_current+1;
            for j=nodes_to_remove(end):size(coordinates,1)
                shift(j)=shift_current;
            end
            index_nodes_new=(1:size(coordinates,1))'-shift;

            elements=index_nodes_new(elements);

                 
            
            coordinates(nodes_to_remove,:)=[];
            
            
            
        else
            nodes2nodes_old=(1:size(coordinates,1))';
            elements2elements_old=(1:size(elements,1))';
        end
        
        
        
end
