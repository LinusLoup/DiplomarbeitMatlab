%% initiales Gitter:
[ P0,EB0,T0 ] = initmesh('cirsg');

%% globale Verfeinerung:
[ P1,EB1,T1 ] = refinemesh('cirsg',P0,EB0,T0);
subplot(3,2,1); pdemesh(P0,EB0,T0);
subplot(3,2,2); pdemesh(P1,EB1,T1);

%% lokale Verfeinerung:
% Punkte, die den Nullpunkt enthalten:
ind = find(find(P0(1,:) == 0) == find(P0(2,:) == 0));
find_ind = find(P0(1,:) == 0);
ind_point = find_ind(ind);
TR0 = [find(T0(1,:) == ind_point),find(T0(2,:) == ind_point),find(T0(3,:) == ind_point)]';
[ P2,EB2,T2 ] = refinemesh('cirsg',P0,EB0,T0,TR0);

subplot(3,2,3); pdemesh(P0,EB0,T0);
subplot(3,2,4); pdemesh(P2,EB2,T2);

%% Triangulierung mit Delaunay:

X = [-1 0 1 1 1 0 -1 -1
    -1 -1 -1 0 1 1 1 0]';
tri = delaunay(X);

points = [ -1,  1, 1, -1
         -1, -1, 1,  1 ];
geom(1:2) = [2, 4];
geom(3:6) = points(1,:);
geom(7:10) = points(2,:);
    
dl = decsg(geom');
    
[p,e,t] = initmesh(dl,'Hmax',0.3);
subplot(3,2,5); pdemesh(p,e,t);

pointsl = [ -1,  1, 1, 0, 0, -1
         -1, -1, 1,  1, 0,  0];
geoml(1:2) = [3, 6];
geoml(3:8) = pointsl(1,:);
geoml(9:14) = pointsl(2,:);
    
dll = decsg(geoml');

[pl,el,tl] = initmesh(dll,'Hmax',0.1);
subplot(3,2,6); pdemesh(pl,el,tl);