function vNodal_finer=prolongate(vNodal,elements,edge2nodes)
%prolongates bilinear function defined on rectangles
vNodal_finer=[vNodal; sum(vNodal(edge2nodes),2)/2; sum(vNodal(elements),2)/4];
