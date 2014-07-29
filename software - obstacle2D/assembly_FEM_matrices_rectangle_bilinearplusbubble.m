function[K, M, K_3Dmatrix, M_3Dmatrix]=assembly_FEM_matrices_rectangle_bilinearplusbubble(elementsExtended,hx,hy)

    NE=size(elementsExtended,1); 

    DK=(1/3)*(hx/hy+hy/hx);
    Klocal=zeros(5);
    for i = 1 : 4
        Klocal(i,i)=DK;
    end
    %bubble addon
    Klocal(5,5)=(128/45)*(hx/hy+hy/hx);
 
    MDK1=-(1/3)*(hy/hx)+(1/6)*(hx/hy);
    Klocal(2,1)=MDK1;
    Klocal(1,2)=MDK1;
    Klocal(4,3)=MDK1;
    Klocal(3,4)=MDK1;
    MDK2=-(1/6)*(hy/hx)-(1/6)*(hx/hy);
    Klocal(3,1)=MDK2;
    Klocal(1,3)=MDK2;
    Klocal(4,2)=MDK2;
    Klocal(2,4)=MDK2;
    MDK3=(1/6)*(hy/hx)-(1/3)*(hx/hy);
    Klocal(4,1)=MDK3;
    Klocal(1,4)=MDK3;
    Klocal(3,2)=MDK3;
    Klocal(2,3)=MDK3;

    DM=(1/9)*hx*hy;
    Mlocal=zeros(4);
    for i = 1 : 4
        Mlocal(i,i)=DM;
    end
    
    MDM1=(1/18)*hx*hy;
    Mlocal(2,1)=MDM1;
    Mlocal(1,2)=MDM1;
    Mlocal(4,3)=MDM1;
    Mlocal(3,4)=MDM1;
    Mlocal(3,2)=MDM1;
    Mlocal(2,3)=MDM1;
    Mlocal(4,1)=MDM1;
    Mlocal(1,4)=MDM1;

    MDM2=(1/36)*hx*hy;
    Mlocal(3,1)=MDM2;
    Mlocal(1,3)=MDM2;
    Mlocal(4,2)=MDM2;
    Mlocal(2,4)=MDM2;  
    
    %bubble addon
    Mlocal(1,5)=(1/9)*hx*hy;
    Mlocal(2,5)=(1/9)*hx*hy;
    Mlocal(3,5)=(1/9)*hx*hy;
    Mlocal(4,5)=(1/9)*hx*hy;
    Mlocal(5,1)=(1/9)*hx*hy;
    Mlocal(5,2)=(1/9)*hx*hy;
    Mlocal(5,3)=(1/9)*hx*hy;
    Mlocal(5,4)=(1/9)*hx*hy;
    Mlocal(5,5)=(64/225)*hx*hy;
    
    r=kron(elementsExtended(:,:),ones(1,5));
    c=kron(ones(1,5),elementsExtended(:,:));
    v_K=kron(ones(NE,1),reshape(Klocal,1,5^2));
    v_M=kron(ones(NE,1),reshape(Mlocal,1,5^2));
    
     K_3Dmatrix=zeros(5,5,NE);
     M_3Dmatrix=zeros(5,5,NE);
     for i=1:5
         for j=1:5
            K_3Dmatrix(i,j,:)=Klocal(i,j);
            M_3Dmatrix(i,j,:)=Mlocal(i,j);
         end
     end
     
    K=sparse(r,c,v_K);
    M=sparse(r,c,v_M);
end
