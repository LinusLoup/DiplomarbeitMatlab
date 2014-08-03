function N0_set = N0(uS_values,obstacle_values)
%N0 berechnet die Menge N^0 der Punkte, die in Kontakt mit dem Hindernis
%stehen.
%
%Hierbei sind uS_values die Funktionswerte der Approximation u_S an den
%Punkten P und obstacle_values die zu P gehörenden Funktionswerte bzgl. des
%Hindernisses.
%
%N0_set ist die Menge der Indizes der Punkte P, in denen u_S(P) = phi(P)
%gilt.

N0_set = find(uS_values == obstacle_values);

end