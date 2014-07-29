function [b, b_3Dmatrix]=linearForm_fQ1plusbubble_vQ1(f_nodalRectangleExtended,elementsRectangularExtended,hx,hy)
elementsRectangular=elementsRectangularExtended(:,1:4);
NN=max(max(elementsRectangular));

[dummy,MRectangular, dummy, MRectangular_3Dmatrix]=assembly_FEM_matrices_rectangle_bilinearplusbubble(elementsRectangularExtended,hx,hy);
bRectangleExtended=MRectangular*f_nodalRectangleExtended;
b=bRectangleExtended(1:NN);

b_3Dmatrix=amt(avtam(conv_ma2av(f_nodalRectangleExtended(elementsRectangularExtended)),MRectangular_3Dmatrix));
b_3Dmatrix(5,:,:)=[];
end
