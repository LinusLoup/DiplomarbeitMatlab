function Nplusplus_set = Nplusplus(Nplus_set,edges,rho_E,d_E)
%NPLUSPLUS calculates the set N^{++} of non-contact nodes, where the
%approximate error eps_V is not in contact. 

% Initialization:
Nplusplus_set = zeros(size(Nplus_set));

% evaluation of the entries of N^{++}:
for i = 1:length(Nplus_set)
    [~,E_p] = find(edges == Nplus_set(i));
    
    % verificate, if for all E \in E_p is valid: rho_E >= -d_E
    if all(rho_E(E_p) + d_E(E_p) >= 0)
        Nplusplus_set(i) = Nplus_set(i);
    end
end

% elimination of the zeros:
Nplusplus_set = setdiff(Nplusplus_set,0);

end