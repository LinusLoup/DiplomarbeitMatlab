function[K, M, K_3Dmatrix, M_3Dmatrix]=assembly_FEM_matrices_rectangle_bilinear(elements,hx,hy)

    NE=size(elements,1); 

    DK=(1/3)*(hx/hy+hy/hx);
    Klocal=zeros(4);
    for i = 1 : 4
        Klocal(i,i)=DK;
    end
 
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

    r=kron(elements(:,:),ones(1,4));
    c=kron(ones(1,4),elements(:,:));
    v_K=kron(ones(NE,1),reshape(Klocal,1,16));
    v_M=kron(ones(NE,1),reshape(Mlocal,1,16));
    
     K_3Dmatrix=zeros(4,4,NE);
     M_3Dmatrix=zeros(4,4,NE);
     for i=1:4
         for j=1:4
            K_3Dmatrix(i,j,:)=Klocal(i,j);
            M_3Dmatrix(i,j,:)=Mlocal(i,j);
         end
     end
     
    K=sparse(r,c,v_K);
    M=sparse(r,c,v_M);

    
    
    %stifness matrix for biquadratic elements
    DKQ1=(1/3)*(hx/hy+hy/hx);
    DKQ2=(8/45)*(3*hx/hy+10*hy/hx);
    DKQ3=(8/45)*(10*hx/hy+3*hy/hx);
    
    KQlocal=zeros(9);
    for i = 1 : 4
        KQlocal(i,i)=DKQ1;
    end
    
    KQlocal(5,5)=DKQ2;
    KQlocal(6,6)=DKQ3;
    KQlocal(7,7)=DKQ2;
    KQlocal(8,8)=DKQ3;
    KQlocal(9,9)=(128/45)*(hx/hy+hy/hx);
 
    MDKQ1=-(1/3)*(hy/hx)+(1/6)*(hx/hy);
    KQlocal(2,1)=MDKQ1;
    KQlocal(1,2)=MDKQ1;
    KQlocal(4,3)=MDKQ1;
    KQlocal(3,4)=MDKQ1;
    MDKQ2=-(1/6)*(hy/hx)-(1/6)*(hx/hy);
    KQlocal(3,1)=MDKQ2;
    KQlocal(1,3)=MDKQ2;
    KQlocal(4,2)=MDKQ2;
    KQlocal(2,4)=MDKQ2;
    MDKQ3=(1/6)*(hy/hx)-(1/3)*(hx/hy);
    KQlocal(4,1)=MDKQ3;
    KQlocal(1,4)=MDKQ3;
    KQlocal(3,2)=MDKQ3;
    KQlocal(2,3)=MDKQ3;
    
    KQlocal(1,5)=(1/3)*hx/hy;
    KQlocal(5,1)=(1/3)*hx/hy;
    KQlocal(1,6)=-(1/3)*hy/hx;
    KQlocal(6,1)=-(1/3)*hy/hx;
    KQlocal(1,7)=-(1/3)*hx/hy;
    KQlocal(7,1)=-(1/3)*hx/hy;
    KQlocal(1,8)=(1/3)*hy/hx;
    KQlocal(8,1)=(1/3)*hy/hx;
    KQlocal(1,9)=0;
    KQlocal(9,1)=0;
    
    KQlocal(2,5)=(1/3)*hx/hy;
    KQlocal(5,2)=(1/3)*hx/hy;
    KQlocal(2,6)=(1/3)*hy/hx;
    KQlocal(6,2)=(1/3)*hy/hx;
    KQlocal(2,7)=-(1/3)*hx/hy;
    KQlocal(7,2)=-(1/3)*hx/hy;
    KQlocal(2,8)=-(1/3)*hy/hx;
    KQlocal(8,2)=-(1/3)*hy/hx;
    KQlocal(2,9)=0;
    KQlocal(9,2)=0;
    
    KQlocal(3,5)=-(1/3)*hx/hy;
    KQlocal(5,3)=-(1/3)*hx/hy;
    KQlocal(3,6)=(1/3)*hy/hx;
    KQlocal(6,3)=(1/3)*hy/hx;
    KQlocal(3,7)=(1/3)*hx/hy;
    KQlocal(7,3)=(1/3)*hx/hy;
    KQlocal(3,8)=-(1/3)*hy/hx;
    KQlocal(8,3)=-(1/3)*hy/hx;
    KQlocal(3,9)=0;
    KQlocal(9,3)=0;
    
    KQlocal(4,5)=-(1/3)*hx/hy;
    KQlocal(5,4)=-(1/3)*hx/hy;
    KQlocal(4,6)=-(1/3)*hy/hx;
    KQlocal(6,4)=-(1/3)*hy/hx;
    KQlocal(4,7)=(1/3)*hx/hy;
    KQlocal(7,4)=(1/3)*hx/hy;
    KQlocal(4,8)=(1/3)*hy/hx;
    KQlocal(8,4)=(1/3)*hy/hx;
    KQlocal(4,9)=0;
    KQlocal(9,4)=0;
    
    KQlocal(5,6)=0;
    KQlocal(6,5)=0;
    KQlocal(5,7)=-(8/45)*(-5*hy/hx+3*hx/hy);
    KQlocal(7,5)=-(8/45)*(-5*hy/hx+3*hx/hy);
    KQlocal(5,8)=0;
    KQlocal(8,5)=0;
    KQlocal(5,9)=(16/9)*hy/hx;
    KQlocal(9,5)=(16/9)*hy/hx;
    
    KQlocal(6,7)=0;
    KQlocal(7,6)=0;
    KQlocal(6,8)=(8/45)*(-3*hy/hx+5*hx/hy);
    KQlocal(8,6)=(8/45)*(-3*hy/hx+5*hx/hy);
    KQlocal(6,9)=(16/9)*hx/hy;
    KQlocal(9,6)=(16/9)*hx/hy;
    
    KQlocal(7,8)=0;
    KQlocal(8,7)=0;
    KQlocal(7,9)=(16/9)*hy/hx;
    KQlocal(9,7)=(16/9)*hy/hx;
    
    KQlocal(8,9)=(16/9)*hx/hy;
    KQlocal(9,8)=(16/9)*hx/hy;
    
    
    %mass matrix for biquadratic elements
    DMQ1=(1/9)*hx*hy;
    DMQ2=(8/45)*hx*hy;
    
    MQlocal=zeros(9);
    for i = 1 : 4
        MQlocal(i,i)=DMQ1;
    end
    
    for i = 5 : 8
        MQlocal(i,i)=DMQ2;
    end
   
    MQlocal(9,9)=(64/225)*hx*hy;
 
    MDMQ1=(1/18)*hx*hy;
    MQlocal(2,1)=MDMQ1;
    MQlocal(1,2)=MDMQ1;
    MQlocal(4,3)=MDMQ1;
    MQlocal(3,4)=MDMQ1;
    MQlocal(3,2)=MDMQ1;
    MQlocal(2,3)=MDMQ1;
    MQlocal(4,1)=MDMQ1;
    MQlocal(1,4)=MDMQ1;

    MDMQ2=(1/36)*hx*hy;
    MQlocal(3,1)=MDMQ2;
    MQlocal(1,3)=MDMQ2;
    MQlocal(4,2)=MDMQ2;
    MQlocal(2,4)=MDMQ2;      
    
    MQlocal(1,5)=(1/9)*hx*hy;
    MQlocal(5,1)=(1/9)*hx*hy;
    MQlocal(1,6)=(1/18)*hx*hy;
    MQlocal(6,1)=(1/18)*hx*hy;
    MQlocal(1,7)=(1/18)*hx*hy;
    MQlocal(7,1)=(1/18)*hx*hy;
    MQlocal(1,8)=(1/9)*hx*hy;
    MQlocal(8,1)=(1/9)*hx*hy;
    MQlocal(1,9)=(1/9)*hx*hy;
    MQlocal(9,1)=(1/9)*hx*hy;
    
    MQlocal(2,5)=(1/9)*hx*hy;
    MQlocal(5,2)=(1/9)*hx*hy;
    MQlocal(2,6)=(1/9)*hx*hy;
    MQlocal(6,2)=(1/9)*hx*hy;
    MQlocal(2,7)=(1/18)*hx*hy;
    MQlocal(7,2)=(1/18)*hx*hy;
    MQlocal(2,8)=(1/18)*hx*hy;
    MQlocal(8,2)=(1/18)*hx*hy;
    MQlocal(2,9)=(1/9)*hx*hy;
    MQlocal(9,2)=(1/9)*hx*hy;
    
    MQlocal(3,5)=(1/18)*hx*hy;
    MQlocal(5,3)=(1/18)*hx*hy;
    MQlocal(3,6)=(1/9)*hx*hy;
    MQlocal(6,3)=(1/9)*hx*hy;
    MQlocal(3,7)=(1/9)*hx*hy;
    MQlocal(7,3)=(1/9)*hx*hy;
    MQlocal(3,8)=(1/18)*hx*hy;
    MQlocal(8,3)=(1/18)*hx*hy;
    MQlocal(3,9)=(1/9)*hx*hy;
    MQlocal(9,3)=(1/9)*hx*hy;
    
    MQlocal(4,5)=(1/18)*hx*hy;
    MQlocal(5,4)=(1/18)*hx*hy;
    MQlocal(4,6)=(1/18)*hx*hy;
    MQlocal(6,4)=(1/18)*hx*hy;
    MQlocal(4,7)=(1/9)*hx*hy;
    MQlocal(7,4)=(1/9)*hx*hy;
    MQlocal(4,8)=(1/9)*hx*hy;
    MQlocal(8,4)=(1/9)*hx*hy;
    MQlocal(4,9)=(1/9)*hx*hy;
    MQlocal(9,4)=(1/9)*hx*hy;
    
    MQlocal(5,6)=(1/9)*hx*hy;
    MQlocal(6,5)=(1/9)*hx*hy;
    MQlocal(5,7)=(4/45)*hx*hy;
    MQlocal(7,5)=(4/45)*hx*hy;
    MQlocal(5,8)=(1/9)*hx*hy;
    MQlocal(8,5)=(1/9)*hx*hy;
    MQlocal(5,9)=(8/45)*hx*hy;
    MQlocal(9,5)=(8/45)*hx*hy;
    
    MQlocal(6,7)=(1/9)*hx*hy;
    MQlocal(7,6)=(1/9)*hx*hy;
    MQlocal(6,8)=(4/45)*hx*hy;
    MQlocal(8,6)=(4/45)*hx*hy;
    MQlocal(6,9)=(8/45)*hx*hy;
    MQlocal(9,6)=(8/45)*hx*hy;
    
    MQlocal(7,8)=(1/9)*hx*hy;
    MQlocal(8,7)=(1/9)*hx*hy;
    MQlocal(7,9)=(8/45)*hx*hy;
    MQlocal(9,7)=(8/45)*hx*hy;
    
    MQlocal(8,9)=(8/45)*hx*hy;
    MQlocal(9,8)=(8/45)*hx*hy;
    
    
    
    
end
