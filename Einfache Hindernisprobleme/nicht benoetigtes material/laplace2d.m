data = load('mycircle.mat');
%data = load('mysquare.mat');

axis on;

%% Lastfunktion f:
%fun = @(x,y) zeros(1,length(x));
fun = @(x,y) -1*ones(1,length(x));
%fun = @(x,y) 5*(-x.^2-y.^2);
%fun = @(x,y) -3*x.^2-10*y.^2;


%% Initialisierung und Plotten des Gitters:
h = 0.2;
[p,e,t] = initmesh(data.mycircleg,'Hmax',h);
%[p,e,t] = refinemesh(data.mycircleg,p,e,t);
subplot(2,2,[1,2]);pdemesh(p,e,t);
[ntri] = size(t,2);
[np] = size(p,2);

for j = 1 : ntri
    tri = t(1:3,j);
    poi = p(:,tri);
    adv = (poi(:,1)+poi(:,2)+poi(:,3))/3;
    text(adv(1),adv(2), num2str(j), 'FontSize',9);
end

for i = 1 : np
    text(p(1,i),p(2,i), num2str(i), 'FontSize',9);
end

%% Assemblierung der Daten mit eigenem Assemblierer:
[Q,G,H,R]=assemb(data.mycircleb,p,e);
[A,f] = assemble(p,t,fun,7);

%% Berechnnung der Indizes für Rand- bzw. Nicht-Randpunkte:
non_bound_ind = find(prod(H == zeros(size(H))));
l = length(non_bound_ind);
bound_ind = find(sum(H));

%% Elimination der Zeilen und Spalten, die durch Dirichlet-Rand vorgegeben:
A_new = A(non_bound_ind,non_bound_ind);
f_new = f(non_bound_ind);

%% Lösung des LGS und Einbau der Dirichletdaten:
uu = A_new\f_new;
u = dirichlet_boundary(uu,H,R);

subplot(2,2,3); pdesurf(p,t,u);

%% Vergleich der Daten anhand des Assemblierers von Matlab:
U=assempde(data.mycircleb,p,e,t,1,0,-1);
[K,M,F,a,b,c,d]=assempde(data.mycircleb,p,e,t,1,0,-1);
subplot(2,2,4); pdesurf(p,t,U);