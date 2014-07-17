%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program solves the obstacle problem
% 
%               -div(grad(u)) >= f
%                           u >= g
%    (div(grad(u)) + f)(u - g) = 0
%
% on the 2D unit square for homogeneous Dirichlet 
% boundary conditions with finite element discretisation. 
% The solution of the resulting variational inequality 
% of the first kind is performed by
% 
% Command: "obstacle(item)"
% - a monotone Gauss-Seidel method      if item == 1
% - an active set method                if item == 2
% - an Augmented Lagrangean method      if item == 3
%
% Referenes: Atkinson/Han: 
% "Theoretical numerical analysis : a functional analysis framework"
% For active set method we refer to: Hüeber/Wohlmuth:
% "A primal-dual active set strategy for non-linear multibody
%  contact probelms" (Computer Methods in Applied Mechanics and
%  Engineering 194:3147-3166, 2005)
% 
% Author: Stefan Hüeber
% Email:  hueeber@mathematik.uni-stuttgart.de
% Date:   23.05.2007
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obstacle_new(item)

%%%% User specified Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 24; % number of discretisation points in each direction
ffunc = @(x,y) zeros(N,1);%inline('-16*x.*(1-x).*y.*(1-y)','x','y'); % right hand side f
gfunc = @(x,y) -(x-0.5).^2-(y-0.5).^2+0.1; %inline('-0.008*(1+2*x+2*y)','x','y');     % obstacle function g
maxit = 1000; % maximum number of iterations
% Augmented Lagrangean method
eps = 0.001;
alpha = 1.0;   % 0.0: Pure Penalty
               % 1.0: Augmented Lagrangean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = 1/(N-1);
[X,Y] = meshgrid(0:h:1); % discretisation of unit square

% Assemble finite element matrix
A = gallery('poisson',N); % system matrix

% Assemble right hand side
ff = h*h*(ffunc(X,Y)); 
ft = ff';
f = ft(:)

% Discrete obstacle values
gh = gfunc(X,Y);
gt = gh';
g = gt(:); 

% Insert Dirichlet datas on the whole boundary
for (i=1:N)
	A(i,:)   = zeros(1,N*N); A(i,i)     = 1; f(i)   = 0;
	ind = i+(N-1)*N;
	A(ind,:) = zeros(1,N*N); A(ind,ind) = 1; f(ind) = 0;
	ind = (i-1)*N+1;
	A(ind,:) = zeros(1,N*N); A(ind,ind) = 1; f(ind) = 0;
	ind = i*N;
	A(ind,:) = zeros(1,N*N); A(ind,ind) = 1; f(ind) = 0;
end

% initial values
uold = zeros(N*N,1);
lambda = zeros(N*N,1);
u = uold;
it = 0;

switch item
	
case 1 % monotone Gauss-Seidel method
	while (it < maxit)
		for (i=1:N*N)
			u(i) =  (f(i) - A(i,1:i-1)*u(1:i-1) ...
							 - A(i,i+1:N*N)*uold(i+1:N*N)) ...
							/ A(i,i);
			if u(i) < g(i);
				u(i) = g(i);
			end
		end
		it = it+1;
		relerr = norm(u-uold)/norm(u)	;
		
		if (relerr < 1e-06) 
			break;
		else
			uold = u;
		end	
	end
	
case 2 % active set method
	lambda = zeros(N*N,1);
	I = ones(N*N,1); % active indeces
	while (it < maxit)
		it = it+1;
		lambda0 = lambda;
		u0 = u;
		I0 = I;
		I = lambda + g - u > 0; % update active indeces
		if (I == I0) && (it > 1)
			break;
		end
		Mat = A;
		rhs = f;
		id = find(I);
		for i = 1:length(id)
			Mat(id(i),:) = 0;
			Mat(id(i),id(i)) = 1;
			rhs(id(i)) = g(id(i)); 
		end
		u = Mat\rhs;
		lambda = A*u-f;
    end
    
case 3 % Augmented Lagrangean method (Grossman/Roos, p. 453/454)
    itimax = 0;
    rho = 1.0/eps;
    while (it < maxit)
        iti = 0;
        uoldi = u;
        while (iti < maxit)
            tmp = lambda + rho*(g-uoldi);
            Ait = A + rho*diag(tmp > 0);
            fit = f + (lambda + rho*g).*(tmp>0);
            u = Ait\fit;
            iti = iti + 1;
            relerri = norm(u-uoldi)/norm(u);
            if (relerri < 1e-11) 
                if (iti > itimax)
                    itimax = iti;
                end
                break;
            else
                uoldi = u;
            end	
        end
	    it = it+1;
        relerr(it) = norm(u-uold)/norm(u);
        if (relerr(it) < 1e-10) 
            break;
        else
            lambda = alpha*(lambda + rho*(g-u));
            lambda = lambda.*(lambda>0);
            uold = u;
        end	
    end    
end


% display solution
for (i=1:N)
	U(i,1:N) = u((i-1)*N+1:i*N);
end

figure(1)
surf(X,Y,U);