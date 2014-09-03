function [wi,gauss,ansatz_values] = quad_tri(points,ansatz_fun,num_of_nodes)
%QUAD_TRI berechnet die für die Quadratur benötigten Gewichte,
%Gaußpunkte, sowie die Funktionswerte der Ansatzfunktionen und übergibt
%diese.
%Hierbei sind x und y die x- bzw. y-Werte der Eckpunkte des Dreiecks,
%ansatz_fun ein Vektor aus Ansatzfunktionen und num_of_nodes die Anzahl der
%Knoten für die Quadratur.


%% x- und y-Werte der Punkte:
x = points(1,:);
y = points(2,:);

switch num_of_nodes
        
        case 3
            %% Bestimmung der Gewichte und der lokalen Punkte xi,eta.
            wi = [1/6,1/6,1/6];
            local_poi = [1/2,1/2,0;0,1/2,1/2];
            xi = local_poi(1,:);
            eta = local_poi(2,:);
            
        case 7
            %% Hilfsgewichte:
            wh = [(155-sqrt(15))/2400,(155+sqrt(15))/2400];
            
            %% Bestimmung der Gewichte und der lokalen Punkte xi,eta.
            wi = [9/80,wh(1),wh(1),wh(1),wh(2),wh(2),wh(2)];
            local_poi = [1/3,(6-sqrt(15))/21,(9+2*sqrt(15))/21,(6-sqrt(15))/21,...
                        (6+sqrt(15))/21,(9-2*sqrt(15))/21, (6+sqrt(15))/21
                        1/3,(6-sqrt(15))/21,(6-sqrt(15))/21,(9+2*sqrt(15))/21,...
                        (6+sqrt(15))/21,(6+sqrt(15))/21,(9-2*sqrt(15))/21];
            xi = local_poi(1,:);
            eta = local_poi(2,:);
            
        otherwise
            error('keine Quadraturformel bekannt');
            
end  

%% Berechnung der Gaußpunkte aus xi,eta mit den Eckpunkten x,y.
gauss(1,:) = x(1)+(x(2)-x(1))*xi+(x(3)-x(1))*eta;
gauss(2,:) = y(1)+(y(2)-y(1))*xi+(y(3)-y(1))*eta;

%% Berechnung der Ansatz-fkt.-Werte in den lokalen Koordinaten. 
ansatz_values = ansatz_fun(xi,eta);

end

