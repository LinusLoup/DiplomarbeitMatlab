function [b, b_3Dmatrix]=linearForm_fQ0_vQ1(fQ0,elementsRectangular,hx,hy)

area=(1/4)*hx*hy;
b=sparse(elementsRectangular,ones(size(elementsRectangular)),area*kron(fQ0,ones(1,4)));
b_3Dmatrix=conv_ma2av(area*kron(fQ0,ones(1,4)));

end
