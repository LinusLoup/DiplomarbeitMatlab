function GRADvRectangular_in_midpoint=evaluate_GRADvQ1(vRectangular,elementsRectangular,hx,hy)

GRADmatrix=[-1/2/hx -1/2/hy; ...
            1/2/hx  -1/2/hy; ...
            1/2/hx   1/2/hy; ...
            -1/2/hx   1/2/hy];
        
GRADvRectangular_in_midpoint=vRectangular(elementsRectangular)*GRADmatrix;

end
