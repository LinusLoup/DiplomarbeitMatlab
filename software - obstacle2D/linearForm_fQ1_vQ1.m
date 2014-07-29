function [b, b_3Dmatrix]=linearForm_fQ1_vQ1(fQ1,elementsRectangular,hx,hy)

[dummy,MRectangular, dummy, MRectangular_3Dmatrix]=assembly_FEM_matrices_rectangle_bilinear(elementsRectangular,hx,hy);
b=MRectangular*fQ1;
b_3Dmatrix=amt(avtam(conv_ma2av(fQ1(elementsRectangular)),MRectangular_3Dmatrix));


end
