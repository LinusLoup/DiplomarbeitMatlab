function z = my_fun( x,y )
%MY_FUN Summary of this function goes here
%   Detailed explanation goes here

z = zeros(size(x));

for i = 1:length(x)
r = sqrt(x(i)^2+y(i)^2);

if r >= 1
    z(i) = r.^2/2-log(r)-1/2;
else
    z(i) = 0;
end

end

end

