function [A,areas,dphi]=stifness_matrixP1_2D(elements,coordinates)
NE=size(elements,1); %number of elements
DIM=size(coordinates,2); %problem dimension

%particular part for a given element in a given dimension
NLB=3; %number of local basic functions, it must be known!
coord=zeros(DIM,NLB,NE);
for d=1:DIM
    for i=1:NLB
        coord(d,i,:)=coordinates(elements(:,i),d);
    end
end   
IP=[1/3;1/3];
[dphi,jac] = phider(coord,IP,'P1'); %integration rule, it must be known! 
dphi = squeeze(dphi); 
areas=abs(squeeze(jac))/factorial(DIM);
Alocal_3Dmatrix=astam(areas',amtam(dphi,dphi));
Y_3Dmatrix=reshape(repmat(elements,1,NLB)',NLB,NLB,NE);

%copy this part for a creation of a new element
X_3Dmatrix=permute(Y_3Dmatrix,[2 1 3]);
A=sparse(X_3Dmatrix(:),Y_3Dmatrix(:),Alocal_3Dmatrix(:));  