function [K_RT0, M_RT0, K_RT0_3Dmatrix, M_RT0_3Dmatrix]=assembly_FEM_matrices_rectangle_RT0(elements2edges,hx,hy)

    NE=size(elements2edges,1); 

    %Raviart Thomas elements of order 1 - RT0
    Klocal_RT0=[hx/hy -1 -hx/hy 1; ...
                -1 hy/hx 1 -hy/hx; ...
                -hx/hy 1 hx/hy -1; ...
                1 -hy/hx -1 hy/hx];

    Mlocal_RT0=(1/6)*[2*hx*hy 0 hx*hy 0; ...
                      0 2*hx*hy 0 hx*hy; ...
                      hx*hy 0 2*hx*hy 0; ...
                      0 hx*hy 0 2*hx*hy];

                  
    K_RT0_3Dmatrix=zeros(4,4,NE);
    M_RT0_3Dmatrix=zeros(4,4,NE);
    for i=1:4
        for j=1:4
            K_RT0_3Dmatrix(i,j,:)=Klocal_RT0(i,j);
            M_RT0_3Dmatrix(i,j,:)=Mlocal_RT0(i,j);
         end
     end
                             
                  
    r_RT0=kron(elements2edges(:,:),ones(1,4));
    c_RT0=kron(ones(1,4),elements2edges(:,:));
    v_K_RT0=kron(ones(NE,1),reshape(Klocal_RT0,1,16));
    v_M_RT0=kron(ones(NE,1),reshape(Mlocal_RT0,1,16));

    K_RT0=sparse(r_RT0,c_RT0,v_K_RT0);
    M_RT0=sparse(r_RT0,c_RT0,v_M_RT0);

end
