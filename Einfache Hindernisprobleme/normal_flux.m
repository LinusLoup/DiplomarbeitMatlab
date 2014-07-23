function j_E = normal_flux(gradu_S,n)
%NORMAL_FLUX berechnet den Normalenfluß über eine Kante E mit
%Normalenvektor n. Hierbei ist gradu_S eine (2x2)-Matrix, die in jeder
%Zeile i den Gradient von u_S bzgl. des Dreiecks T_i enthält. Der
%Normalenvektor n ist hier ein Spalteneinheitsvektor.

normal = gradu_S*n;
j_E = normal(2)-normal(1);

end

