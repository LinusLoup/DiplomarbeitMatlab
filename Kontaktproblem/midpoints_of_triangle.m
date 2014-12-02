function [midpoints,mid_triangle] = midpoints_of_triangle(triangle,points)
%MIDPOINTS_OF_TRIANGLE evaluates the midpoints of the edges of the triangle T. The first and second row stores the x- and y-values of the midpoints and the third and fourth row stores the indices of the nodes, that are the start and end of the edge. The matrix mid_triangle restores the indices of the midpoints in the rows of every triangle by the columns (by positive mathematical orientation).

[nt] = size(triangle,2);
mid_triangle = zeros(3,nt);
midpoints = zeros(4,1);
ind_counter = 1;

for i = 1:nt
    % determination of the nodes of a triangle:
    tri = triangle(1:3,i);
    poi = points(:,tri);
    
    % evaluation of the midpoints to this triangle:
    mid_poi = [1/2*(poi(:,1)+poi(:,2)), 1/2*(poi(:,2)+poi(:,3)),...
        1/2*(poi(:,1)+poi(:,3));tri([1,2]),tri([2,3]),tri([1,3])];
    
    % verification, if the midpoints have been already calculated:
    [glob,loc] = ismember(midpoints(1:2,:)',mid_poi(1:2,:)','rows');
    global_ind = find(glob);
    local_ind = loc(global_ind);
    
    % case distinction and determination of the vector of midpoints:
    switch length(global_ind)
        case 1
            ind = setdiff([1,2,3],local_ind);
            mid_triangle(local_ind,i) = global_ind;
            mid_triangle(ind,i) = [ind_counter,ind_counter+1];
            midpoints(:,[ind_counter,ind_counter+1]) = mid_poi(:,ind);
            ind_counter = ind_counter + 2;
        case 2
            ind = setdiff([1,2,3],local_ind);
            mid_triangle(local_ind,i) = global_ind;
            mid_triangle(ind,i) = ind_counter;
            midpoints(:,ind_counter) = mid_poi(:,ind);
            ind_counter = ind_counter + 1;
        case 3
            mid_triangle(local_ind,i) = global_ind;
        otherwise
        mid_triangle(:,i) = [ind_counter,ind_counter+1,ind_counter+2];
        midpoints(:,ind_counter:ind_counter+2) = mid_poi;
        ind_counter = ind_counter + 3;
    end
end

end