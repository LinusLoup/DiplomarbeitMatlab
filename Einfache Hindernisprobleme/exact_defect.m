function [eps_V_exact,rho_E,d_E,a_phi] = exact_defect(points,triangle,...
    mid_triangle,rho_s,u_S_mid,z_obs_midpoints)
    
    ntri = size(triangle,2);
    nmp = size(u_S_mid,1);

    dxi_bubble = @(xi,eta) [4-8*xi-4*eta; 4*eta;-4*eta];
    deta_bubble = @(xi,eta) [-4*xi; 4*xi; 4-4*xi-8*eta];
    mymidtri = [mid_triangle;zeros(1,ntri)];
    a_phi = zeros(nmp,1);

for i = 1 : nmp
    index = find(mymidtri == i);
    bubble_ind = mod(index,4);
    triangle_ind = ceil(index/4);
    
    poi = points(:,triangle(1:3,triangle_ind));
    x = poi(1,:);
    y = poi(2,:);
    J = (x(2)-x(1))*(y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1));
    a = ((x(3)-x(1))^2 + (y(3)-y(1))^2);
    b = -((y(2)-y(1))*(y(3)-y(1)) + (x(2)-x(1))*(x(3)-x(1)));
    c = ((x(2)-x(1))^2 + (y(2)-y(1))^2);
    
    [wi,~,val_dxi_bubble] = quad_tri(poi,dxi_bubble,7);
    [~,~,val_deta_bubble] = quad_tri(poi,deta_bubble,7);
    
    for j = 1 : length(bubble_ind)
        a_phi(i) = a_phi(i) + 1/J* sum(wi.*(a*val_dxi_bubble...
            (bubble_ind(j),:).^2 + b*val_dxi_bubble(bubble_ind(j),:).*...
            val_deta_bubble(bubble_ind(j),:) * 2 ...
            +c*val_deta_bubble(bubble_ind(j),:).^2));
    end
    
end

a_phi = a_phi.^(1/2);
rho_E = rho_s./a_phi;
d_E = (u_S_mid-z_obs_midpoints).*a_phi;

eps_V_exact = max(-d_E,rho_E)./a_phi;

end