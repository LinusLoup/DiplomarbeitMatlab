function N0_set = N0(uS_values,inner_points,obstacle_values)
%N0 determines the set N^0 of the inner points, which are in contact with 
%the obsacle. It will be done by simply compare the function values of u_S 
%and the obstacle.

index = find(abs(uS_values - obstacle_values) < 0.00001);
N0_set = intersect(index,inner_points);

end