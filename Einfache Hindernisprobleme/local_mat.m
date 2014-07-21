function [S,f_local,J] = local_mat(points,uS_local,options)
%LOCAL_MAT berechnet die Funktionaldeterminante J und eine lokale
%Steifigkeitsmatrix S.
%
%Als Eingabeparameter werden die Punkte des lokalen Dreiecks points und
%ein String options, der angibt, welche lokale Steifigkeitsmatrix
%zur�ckgegeben werden soll, erwartet. options kann die Werte 'linear' und
%'bubble' bzw. 'quadratic' annehmen.

%% Initialisierung von f_local:
f_local = zeros(3,1);

%% x- und y-Werte der Punkte:
x = points(1,:);
y = points(2,:);

%% Funktionaldeterminante und weitere Transformationsfaktoren, um vom 
%% Referenzelement auf ein lokales Dreieck T mit Ecken (x1,y1),(x2,y2),
%% (x3,y3) zu transformieren:
J = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));
a = 1/J * ((x(3)-x(1))^2 + (y(3)-y(1))^2);
b = -1/J * ((y(2)-y(1))*(y(3)-y(1)) + (x(2)-x(1))*(x(3)-x(1)));
c = 1/J * ((x(2)-x(1))^2 + (y(2)-y(1))^2);


switch lower(options)
    case {'linear'}
    %% lokale Steifigkeitsmatrix f�r lin. Ansatz auf einem Dreieck T der 
    %% Triangulierung:
    S = [ a/2+b+c/2, -a/2-b/2, -b/2-c/2;
        -a/2-b/2,    a/2,      b/2;
        -b/2-c/2,    b/2,      c/2     ];
    
    case {'bubble','quadratic'}
    %% Die allg. Steifigkeitsmatrix f�r quadr. Ansatzfunktionen:
    S_quad = 4/3* [ a+b+c, -b-c, b;
                    -b-c,  a+b+c, -a-b;
                      b,    -a-b, a+b+c ];
                  
    S = diag(diag(S_quad));
    %S = S_quad;
    
    %% Gradient von u_S:
    gradu_S = gradu(points,uS_local);
    % Transformation auf das Referenzelement:
    gradu_S = [x(2)-x(1), y(2)-y(1); x(3)-x(1), y(3)-y(1)]*gradu_S';
    
    f_local = -(a*gradu_S(1)+b*gradu_S(2))*[0; 2/3; -2/3]...
        -(b*gradu_S(1)+c*gradu_S(2))*[-2/3; 2/3; 0];
    
    otherwise
        error('Unknown option');
end


end

