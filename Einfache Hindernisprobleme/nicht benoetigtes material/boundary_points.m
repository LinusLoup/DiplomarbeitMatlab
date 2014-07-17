function boundary = boundary_points(geom)
%BOUNDARY_POINTS Summary of this function goes here
%   Detailed explanation goes here

[m,n] = size(geom);
pos_count = 0;
first_param = 0:0.05:1;

%% Zur Vermeidung doppelter Punkte im Übergang von den Randstücken:
param = first_param(1:length(first_param)-1);
param_length = length(param);

%% Berechnung der Randpunkte:
for i = 1:n
    if geom(1,i) == 2
        boundary(1,pos_count*param_length+1:(pos_count + 1)*param_length)...
            = (1.-param).*geom(2,i)+param.*geom(3,i);
        boundary(2,pos_count*param_length+1:(pos_count + 1)*param_length)...
            = (1.-param).*geom(4,i)+param.*geom(5,i);
        pos_count = pos_count + 1;
    elseif geom(1,i) == 4
        a = geom(10,i);
        b = geom(11,i);
        x0 = geom(8,i); 
        y0 = geom(9,i);
        boundary(1,pos_count*param_length+1:(pos_count + 1)*param_length)...
            = a.*cos(param)+ x0;
        boundary(2,pos_count*param_length+1:(pos_count + 1)*param_length)...
            = b.*sin(param)+ y0;
        pos_count = pos_count + 1;
    else
        error('the given geometry is not known');
    end
end

end

