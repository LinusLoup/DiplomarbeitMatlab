function [u,fval] = exact_solution(g,b,fun,obstacle,h)
%EXACT_SOLUTION Summary of this function goes here
%   Detailed explanation goes here

tic;

[p,e,t] = initmesh(g,'Hmax',h);
% [p,e,t] = refinemesh(g,p,e,t);
% [p,e,t] = refinemesh(g,p,e,t);
% [p,e,t] = refinemesh(g,p,e,t);
% [p,e,t] = refinemesh(g,p,e,t);
% [p,e,t] = refinemesh(g,p,e,t);
% [p,e,t] = refinemesh(g,p,e,t);

z_obs_prob = obstacle(p(1,:),p(2,:));

[A,f] = assemble(p,t,fun,7);
[~,~,H,R]=assemb(b,p,e);

[u,fval] = quadprog(A,-f,[],[],H,R,z_obs_prob,[]);



toc;

subplot(2,1,[1,2]);
pdemesh(p,e,t,u);

length(u)
1/2*u'*A*u-u'*f

end