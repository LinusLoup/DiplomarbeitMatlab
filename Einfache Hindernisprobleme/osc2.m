function osc2_val = osc2(Nplusplus_set,N0minus_set,nodes,triangles,...
    edges,f_load)
%OSC2 berechnet den zweiten Oszillationsterm osc2(u_S,psi,f).


%% Initialisierung:
osc2_val = 0;
osc2_sum1 = 0;
osc2_sum2 = 0;
np = length(nodes);
N_set = 1:np;
N_sum2 = setdiff(N_set,union(Nplusplus_set,N0minus_set));
h_P1 = zeros(size(Nplusplus_set));
h_P2 = zeros(size(N_sum2));

%% Berechnung der h_P für die beiden Summen:
for i = 1:length(h_P1)
    % Berechnnung der Kanten, auf denen P liegt:
    [~,E_p] = find(edges(3:4,:)==Nplusplus_set(i));
    
    % Berechnung der Längen von E \in E_p:
    points_1 = nodes(:,edges(3,E_p));
    points_2 = nodes(:,edges(4,E_p));
    points_diff = points_1 - points_2;
    h_P1(i) = max(sqrt(points_diff(1,:).^2 + points_diff(2,:).^2));
end

for i = 1:length(h_P2)
    % Berechnnung der Kanten, auf denen P liegt:
    [~,E_p] = find(edges(3:4,:)==N_sum2(i));
    
    % Berechnung der Längen von E \in E_p:
    points_1 = nodes(:,edges(3,E_p));
    points_2 = nodes(:,edges(4,E_p));
    points_diff = points_1 - points_2;
    h_P2(i) = max(sqrt(points_diff(1,:).^2 + points_diff(2,:).^2));
end


%% Berechnung der Summen:
for i = 1:length(Nplusplus_set)
    % Träger von phi_P, d.h. Dreiecksindizes
    [~,w_p] = find(triangles(1:3,:)==N_sum2(i));
end

for i = 1:length(N_sum2)
    int_h = 0;
    % Träger von phi_P, d.h. Dreiecksindizes
    [~,w_p] = find(triangles(1:3,:)==N_sum2(i));
    
    % Berechnung der Norm:
    for j = 1:length(w_p)
        % Punkte des jeweiligen Dreiecks, auf dem wir uns befinden, mit
        % Jacobi-Determinante:
        mypoi = nodes(:,triangles(1:3,w_p(j)));
        x = mypoi(1,:);
        y = mypoi(2,:);
        J = (x(2)-x(1))*(y(3)-y(1))-(x(3)-x(1))*(y(2)-y(1));
        % Berechnung der Gaußpunkte für die Integration:
        [wi,gauss,~] = quad_tri(mypoi,@(x,y) 0,7);
        % Quadratur:
        int_h = int_h + J * sum(wi.*f_load(gauss(1,:),gauss(2,:)).^2);
    end
    
    osc2_sum2 = osc2_sum2 + h_P2(i)^2*int_h;
end

osc2_val = osc2_val + sqrt(osc2_sum1 + osc2_sum2);

end