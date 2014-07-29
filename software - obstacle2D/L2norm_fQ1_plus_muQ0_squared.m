function [L2norm, L2norm_density]=L2norm_fQ1_plus_muQ0_squared(fQ1,muQ0,elementsRectangular,hx,hy)

[dummy,MRectangular, dummy, MRectangular_3Dmatrix]=assembly_FEM_matrices_rectangle_bilinear(elementsRectangular,hx,hy);
f_plus_mu_3Dvector=conv_ma2av(fQ1+kron(ones(1,size(fQ1,2)),muQ0));
L2norm_distribution=avtamav(f_plus_mu_3Dvector,MRectangular_3Dmatrix,f_plus_mu_3Dvector);
L2norm_density=L2norm_distribution/(hx*hy);
L2norm=sum(L2norm_distribution);
end
