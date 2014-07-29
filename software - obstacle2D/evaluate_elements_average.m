function  elements2midpoint=evaluate_elements_average(elements,coordinates)
        

dummy=coordinates(elements(:,1),:);
for i=2:size(elements,2)
    dummy=dummy+coordinates(elements(:,i),:);
end
   elements2midpoint=dummy/size(elements,2);
   
end
