function [midpoints,mid_triangle] = midpoints_of_triangle(triangle,points)
%MIDPOINTS_OF_TRIANGLE berechnet die Mittelpunkte der Kanten eines Dreiecks
%T und gibt diese in midpoints aus. Dabei werden in den ersten beiden Zeilen
%die x- und y-Koordinate der Punkte gespeichert und in den letzten beiden
%Zeilen die Indizes der Eckpunkte der Seite, auf der sich der Mittelpunkt
%befindet. mid_triangle speichert für ein Dreieck T (Spaltenindex) in den
%Zeilen die Indizes der Mittelpunkte im Dreieck in mathematisch positiver
%Drehrichtung.

[nt] = size(triangle,2);
mid_triangle = zeros(3,nt);
midpoints = zeros(4,1);
ind_counter = 1;

for i = 1:nt
    % Berechnung der Eckpunkte zu einem Dreieck:
    tri = triangle(1:3,i);
    poi = points(:,tri);
    
    % Berechnung der Mittelpunkte zu diesem Dreieck:
    mid_poi = [1/2*(poi(:,1)+poi(:,2)), 1/2*(poi(:,2)+poi(:,3)),...
        1/2*(poi(:,1)+poi(:,3));tri([1,2]),tri([2,3]),tri([1,3])];
    
    % Überprüfung, ob Mittelpunkte schonmal berechnet wurden:
    [glob,loc] = ismember(midpoints(1:2,:)',mid_poi(1:2,:)','rows');
    global_ind = find(glob);
    local_ind = loc(global_ind);
    
    % Die Fallunterscheidung und Einarbeitung der Mittelpunkte in den
    % Vektor:
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