function [J,S] = local_mat(points,options)
%LOCAL_MAT berechnet die Funktionaldeterminante J und eine lokale
%Steifigkeitsmatrix S.
%
%Als Eingabeparameter werden die Punkte des lokalen Dreiecks points und
%ein String options, der angibt, welche lokale Steifigkeitsmatrix
%zurückgegeben werden soll, erwartet. options kann die Werte 'linear' und
%'bubble' bzw. 'quadratic' annehmen.


%% x- und y-Werte der Punkte:
x = points(1,:);
y = points(2,:);


%% Funktionaldeterminante und weitere Transformationsfaktoren, um vom 
%% Referenzelement auf ein lokales Dreieck T mit Ecken (x1,y1),(x2,y2),
%% (x3,y3) zu transformieren:
J = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));
a = ((x(3)-x(1))^2 + (y(3)-y(1))^2);
b = -((y(2)-y(1))*(y(3)-y(1)) + (x(2)-x(1))*(x(3)-x(1)));
c = ((x(2)-x(1))^2 + (y(2)-y(1))^2);


%% Berechnung von grad_u:
grad_u = 1/J * ([y(3)-y(1), y(1)-y(2); x(1)-x(3), x(2)-x(1)]*[z(2)-z(1);...
    z(3)-z(1)])';


switch lower(options)
    case {'linear'}
    %% lokale Steifigkeitsmatrix für lin. Ansatz auf einem Dreieck T der 
    %% Triangulierung:
    S = 1/J * [ a/2+b+c/2, -a/2-b/2, -b/2-c/2;
                    -a/2-b/2,    a/2,      b/2;
                    -b/2-c/2,    b/2,      c/2     ];
    
    case {'bubble','quadratic'}
    %% Die allg. Steifigkeitsmatrix für quadr. Ansatzfunktionen:
    S_quad = 4/3*1/J * [ a+b+c, -b-c, c;
                        -b-c,  a+b+c, -a-c;
                          c,    -a-c, a+b+c ];
                  
    S = diag(diag(S_quad));
    %S = S_quad;
    
    
    otherwise
        error('Unknown option');
end


end

