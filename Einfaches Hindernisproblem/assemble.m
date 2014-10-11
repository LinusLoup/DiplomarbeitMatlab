function [A,f] = assemble(points,triangle,fun,num_of_nodes,option,u_S)
%ASSEMBLE Summary of this function goes here
%   Detailed explanation goes here


%% neue Zuordnung: Dreicke <-> Mittelpunkte
[midpoints,midtriangle] = midpoints_of_triangle(triangle,points);

%% Initialisierungen der Dimensionen und Endergebnisse:
np = size(points,2);
nt = size(triangle,2);
nmp = size(midpoints,2);

if nargin == 4
    option = 'linear';
    u_S = zeros(np,1);
end

if nargin == 5
    u_S = zeros(np,1);
end


%% Beginn der Assemblierung:
switch lower(option)
    case {'linear'}
    %% lineare Hutfunktionen auf dem Referenzelement:
    hat = @(xi,eta) [1-xi-eta; xi; eta];
    
    
    %% Initialisierung der globalen Gr��en:
    A = sparse(np,np);
    f = sparse(np,1);
    my_tri = triangle(1:3,:);
    
    case {'bubble','quadratic'}

    %% bubble-Funktionen auf dem Referenzelement:
    hat = @(xi,eta) [4*xi.*(1-xi-eta); 4*xi.*eta; 4*eta.*(1-xi-eta)];
    
    
    %% Initialisierung der globalen Gr��en:
    A = sparse(nmp,nmp);
    f = sparse(nmp,1);
    my_tri = midtriangle;
    
end

%% Schleife �ber die Dreiecke zur Assemblierung:
for i = 1:nt
    poi = points(:,triangle(1:3,i));
    u_S_loc = u_S(triangle(1:3,i));
    tri = my_tri(:,i);
    
    %% Berechnung der lok.,lin. Steifigkeitsmatrix und der Fkt.-Determinante:
    [S,fl,J] = local_mat(poi,u_S_loc,option);
        
    %% Berechnung des lokalen Vektors fl:
    % Quadraturformeln von https://www-user.tu-chemnitz.de/~rens/lehre/...
    % archiv/numerik1_11SS/folien/gaussxd.pdf
    [wi,gauss,ansatz_values] = quad_tri(poi,hat,num_of_nodes);
            
    %Quadratur auf Dreieck:
    for k = 1:3
        fl(k) = fl(k)+J*sum(wi.*fun(gauss(1,:),gauss(2,:)).*...
            ansatz_values(k,:));
    end
    
    %% Assemblierung in die globale Steifigkeitsmatrix A und den Vektor f:
    for k = 1:3
        for l = 1:3
            A(tri(k),tri(l)) = A(tri(k),tri(l))+S(k,l);
        end
    
        f(tri(k)) = f(tri(k))+fl(k);
    end
end

end
