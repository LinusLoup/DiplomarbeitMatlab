function [osc1_val,osc1_vec] = osc1(N0plus_set,obstacle_values,nodes,...
    triangles,uS_values)
%OSC1 berechnet die Hindernisoszillation osc_1(u_S,phi), d.h. von der
%approximativen Lösung u_S und dem Hindernis psi abhängig.
%
%Hierbei setzen wir voraus, dass das Hindernis \psi affin ist, ansonsten
%stellt die Berechnung nur eine Annäherung an den eigentlichen Term dar.


%% Initialierung:
osc1_vec = zeros(length(nodes),1);

%% Berechnung der Summe über die Punkte aus N0+:
for i = 1:length(N0plus_set)
    % Hilfsvariable für das Integral (die Norm) der Gradienten:
    int_h = 0;
    % Träger der Ansatzfunktion & Index der lokalen Ansatzfunktion:
    [~,w_p] = find(triangles(1:3,:)==N0plus_set(i));
    
    % Berechnung der einzelnen Normen (Integrale):
    for j = 1:length(w_p)
        % Eckpunkte des Dreiecks, Funktionaldet. und Faktoren für Trafo:
        p_index = triangles(1:3,w_p(j));
        mypoi = nodes(:,p_index);
        x = mypoi(1,:);
        y = mypoi(2,:);
        J = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));
        a = 1/J * ((x(3)-x(1))^2 + (y(3)-y(1))^2);
        b = -1/J * ((y(2)-y(1))*(y(3)-y(1)) + (x(2)-x(1))*(x(3)-x(1)));
        c = 1/J * ((x(2)-x(1))^2 + (y(2)-y(1))^2);
        
        % Berechnung der Gradienten von psi und u_S:
        grad_psi = gradu(mypoi,obstacle_values(p_index));
        grad_u = gradu(mypoi,uS_values(p_index));
        
        % Transformation der Gradienten auf das Referenzdreieck:
        grad_psi = [x(2)-x(1), y(2)-y(1); x(3)-x(1), y(3)-y(1)]*grad_psi';
        grad_u = [x(2)-x(1), y(2)-y(1); x(3)-x(1), y(3)-y(1)]*grad_u';
        
            % Da das Integral über eine konstante Funktion gezogen wird,
            % integrieren wir einfach durch Multiplikation des Funktions-
            % wertes (Höhe) mal die Grundfläche (Fläche des Ref-Dreiecks):
            int_h = int_h + 1/2*(a*grad_psi(1)^2+2*b*grad_psi(1)...
                *grad_psi(2)+c*grad_psi(2)^2-2*(a*grad_psi(1)*grad_u(1)...
                +b*(grad_psi(1)*grad_u(2)+grad_psi(2)*grad_u(1))...
                +c*grad_psi(2)*grad_u(1))+a*grad_u(1)^2+2*b*(grad_u(1)...
                *grad_u(2))+c*grad_u(2)^2);
    end
    
    osc1_vec(N0plus_set(i)) = int_h;
end

% Wurzel ziehen aus der Summe der Integrale:
osc1_val = sqrt(sum(osc1_vec));

end