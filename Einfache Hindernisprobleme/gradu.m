function grad_u = gradu(points,zvalues)
%GRADU berechnet für gegebene (x,y,z)-Werte den Gradienten der Funktion
%u=a*x+b*y+c, d.h. grad u=(a,b).


%% x- und y-Werte der Punkte:
x = points(1,:);
y = points(2,:);


%% Funktionaldeterminante und weitere Transformationsfaktoren, um vom 
%% Referenzelement auf ein lokales Dreieck T mit Ecken (x1,y1),(x2,y2),
%% (x3,y3) zu transformieren:
J = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));


%% Berechnung von grad_u:
grad_u = 1/J * ([y(3)-y(1), y(1)-y(2); x(1)-x(3), x(2)-x(1)]*...
    [zvalues(2)-zvalues(1);zvalues(3)-zvalues(1)])';

end