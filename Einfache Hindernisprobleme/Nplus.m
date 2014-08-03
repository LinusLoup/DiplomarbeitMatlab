function Nplus_set = Nplus(contact_set,nodes)
%NPLUS berechnet N^+, die Menge der Knotenindizes, die nicht in Kontakt
%stehen.
%
%Dabei ist contact_set die Indexmenge N^0 der Kontaktknoten und nodes die
%Menge der Knoten.

np = length(nodes);
Nplus_set = setdiff(1:np,contact_set);

end