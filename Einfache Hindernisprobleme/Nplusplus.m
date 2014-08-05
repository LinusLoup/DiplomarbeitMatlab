function Nplusplus_set = Nplusplus(Nplus_set,edges,rho_E,d_E)
%NPLUSPLUS Summary of this function goes here
%   Detailed explanation goes here

%% Initialisierung:
Nplusplus_set = zeros(size(Nplus_set));

%% Berechnung der Einträge von N++:
for i = 1:length(Nplus_set)
    [~,E_p] = find(edges == Nplus_set(i));
    
    % Überprüfung, ob alle für alle E \in E_p gilt: rho_E >= -d_E
    if all(rho_E(E_p) + d_E(E_p) >= 0)
        Nplusplus_set(i) = Nplus_set(i);
    end
end

Nplusplus_set = setdiff(Nplusplus_set,0);

end