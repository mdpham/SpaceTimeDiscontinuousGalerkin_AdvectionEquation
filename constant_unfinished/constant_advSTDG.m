function [U,err] = advSTDG(N,Nt,T,movingmesh,testcase,plot)
% TEST CASE
% f = u_t + u_x, from the PDE
if strcmp(testcase,'constant') % CONSTANT SOLUTION
    uexact = @(t,x) 1;
    f = @(t,x) 0;
elseif strcmp(testcase,'linear') % LINEAR SOLUTION
    uexact = @(t,x) t + x;
    f = @(t,x) 2;
elseif strcmp(testcase,'nonlinear') % NONLINEAR SOLUTION
%     uexact = @(t,x) sin(2*pi*x) + t;
%     f = @(t,x) 2*pi*cos(2*pi*x) + 1;
    uexact = @(t,x) cos(2*pi*(x-t));
    f = @(t,x) 0;
elseif strcmp(testcase,'full')
    uexact = @(t,x) sin(2*pi*x)*sin(3*pi*t); % PROBLEM SOLUTION
    f = @(t,x) 3*pi*sin(2*pi*x)*cos(3*pi*t) + 2*pi*cos(2*pi*x)*sin(3*pi*t);    
else
    error('no such testcase');
end

% TIME VARIABLE
t = 0;
dt = 1/Nt;

% COMPUTATIONAL DOMAIN
% Flow domain
if strcmp(movingmesh,'fixed')
    a = @(t) 0;
    b = @(t) 1;
elseif strcmp(movingmesh,'linear')
    a = @(t) -0.2*t;
    b = @(t) 1+0.2*t;
elseif strcmp(movingmesh,'full')
    a = @(t) sin(2*pi*t)/10;
    b = @(t) exp(-t);
else
    error('no such mesh')
end
flowdomain = @(t) linspace(a(t),b(t),N+1);
% Construct the top, bottom flow domains for a timeslab at time t
    function [flowdom_n, flowdom_n1] = constructtimeslab(t)
        flowdom_n = flowdomain(t);
        flowdom_n1 = flowdomain(t+dt);
    end
% Grid velocity for the side facet at point x_j
    function dx = griddx(t,j)
        [flowdom_n, flowdom_n1] = constructtimeslab(t);
        dx = flowdom_n1(j) - flowdom_n(j);
    end
    function vg = gridvel(t,j)
        dx = griddx(t,j);
        vg = dx/dt;
        if 1 < vg; error(' check a - vg > 0'); end
    end
% Get the four points of a space-time element as well as top/bottom lengths
    function [x_jn, x_j1n, x_jn1, x_j1n1, dx_jn, dx_jn1] = stvertices(j,t)
        [flowdom_n, flowdom_n1] = constructtimeslab(t);
        x_jn = flowdom_n(j);
        x_j1n = flowdom_n(j+1);
        x_jn1 = flowdom_n1(j);
        x_j1n1 = flowdom_n1(j+1);
        dx_jn = x_j1n - x_jn;
        dx_jn1 = x_j1n1 - x_jn1;
    end

% REFERENCE MAPPINGS: returns anonymous functions to be used in quadrature
% (x0,x1) physical time and space respectively
% (xi0,xi1) reference time and space respectively
%   Spatial reference for element j
%     function x = F_K_n(j,t)
%         [x_jn, x_j1n, x_jn1, x_j1n1, dx_jn, dx_jn1] = stvertices(j,t);
%         x = @(xi1) x_jn*0.5*(1-xi1) + x_j1n*0.5*(1+xi1);
%     end
%   Space-time reference for element j
    function [x0,x1] = G_K_n(j,tn)
        [x_jn, x_j1n, x_jn1, x_j1n1, dx_jn, dx_jn1] = stvertices(j,t);
        x0 = @(xi0,xi1) tn + 0.5*dt*(1+xi0); 
        x1 = @(xi0,xi1) ...
            0.5*(1-xi0)*(x_jn*0.5*(1-xi1) + x_j1n*0.5*(1+xi1))...
            + 0.5*(1+xi0)*(x_jn1*0.5*(1-xi1) + x_j1n1*0.5*(1+xi1));
    end
    % Determinant of Jacobian of transformation
    function J = jacobianfun(j,tn)
        [x_jn,x_j1n,x_jn1,x_j1n1,dx_jn,dx_jn1] = stvertices(j,tn);
        J = @(xi0,xi1) (0.5*dt)*0.25*((1-xi0)*dx_jn + (1+xi0)*dx_jn1);
    end

% BOTTOM OF SPACE-TIME ELEMENT FACE INTEGRAL
    function kj_n = facet_kjn(u_prev,j,t)
        [x_jn, x_j1n, x_jn1, x_j1n1, dx_n, dx_n1] = stvertices(j,t);
        kj_n = -quadrature1D(@(xi) 0.5*dx_n*u_prev(j));
    end
% TOP OF SPACE-TIME ELEMENT FACE INTEGRAL
    function kj_n1 = facet_kjn1(j,t)
        [x_jn, x_j1n, x_jn1, x_j1n1, dx_n, dx_n1] = stvertices(j,t);
        kj_n1 = quadrature1D(@(xi) 0.5*dx_n1);
    end
% SIDES OF SPACE-TIME ELEMENT FACET INTEGRALS
% Project manufactured solution on boundary to left side of first element
% The flux that accounts for grid motion and advection is computed here
% Note that on facet S_j1 that normal is (0,1)
    function sj1_v = facet_sj1v(j,t)
        vg = gridvel(t,j+1); %j+1 since on right of jth element
        dx = griddx(t,j+1);
        sj1_v = quadrature1D(@(xi0) (1-vg)*(0.5*dt));
    end
% FIRST ELEMENT LEFT FACET BOUNDARY CONDITION
    function s1_v = facet_s1v(t)
        % Left hand side of element: BC imposed by solution
        [x0,x1] = G_K_n(1,t); % 1 since first element
        vg = gridvel(t,1);
        dx = griddx(t,1);
        u_a = @(xi0) -(1-vg)*uexact(x0(xi0,-1),x1(xi0,-1))*(0.5*dt);
        s1_v = quadrature1D(u_a);
    end
% CONSTRUCT GLOBAL MATRIX FOR TIMESLAB t AND SOLVE
    function u_dof = globalassembly(u_prev,t)
        [flowdom_n, flowdom_n1] = constructtimeslab(t);
        M = zeros(N,N);
        b = zeros(N,1);
        % GLOBAL MATRIX FOR BILINEAR FORM
        for j = 1:N
            % BOTTOM OF SPACE-TIME ELEMENT
            % information advected up, move known information to right hand side
            b(j) = b(j) - facet_kjn(u_prev, j, t);;
            
            % TOP OF SPACE-TIME ELEMENT
            M(j,j) = M(j,j) + facet_kjn1(j, t);
            
            % SIDES OF SPACE-TIME ELEMENT
            if j == 1
                % Left hand side facet: move known info to RHS of matrix eq             
                b(j) = b(j) - facet_s1v(t);
                % Right hand side facet: (vL-vR)f(uL,uR) = (vL-vR)uL
                M(j,j) = M(j,j) + facet_sj1v(j,t);
                M(j+1,j) = M(j+1,j) - facet_sj1v(j,t);; % Subtract since vL-vR
            elseif j == N
                % Right hand BC (outflow): vLf(uL,uR) = vLuL
                M(j,j) = M(j,j) + facet_sj1v(j,t);
            else
                % S_j1 facet for interior spacetime elements
                % For each element we need only do the right facets
                M(j,j) = M(j,j) + facet_sj1v(j,t);
                M(j+1,j) = M(j+1,j) - facet_sj1v(j,t);
            end
        end 
        % SOURCE TERM FOR LINEAR FORM
        % See f term above for problem source term
        S = zeros(N,1);
        for j = 1:N
            J = jacobianfun(j,t);
            [x0,x1] = G_K_n(j,t);
            f_source = @(xi0, xi1) f(x0(xi0,xi1),x1(xi0,xi1))*J(xi0,xi1);
            S(j) = quadrature2D(f_source); 
        end
%         display(M)
%         display(b)
        u_dof = M\(b+S);
    end

% INITIAL CONDITIONS
    function U0 = initialconditions(N)
        [flowdom_n, flowdom_n1] = constructtimeslab(0);
        U0 = zeros(N,1);
        for j = 1:N
            [x0,x1] = G_K_n(j,0);
            % Integral term from LHS cancels out element length |K_jn| term on RHS
            u0 = @(xi1) uexact(x0(-1,xi1),x1(-1,xi1))*0.5;
            U0(j) = quadrature1D(u0);
        end
    end

% COMPUTE L2 ERRORS BY QUADRATURE
    function err = finalerror(u)
        err = 0;
        for j = 1:N
            u_h = u(j);
            approx_err = @(x) abs(u_h - uexact(T,x))^2;
%             J = jacobianfun(j,T);
            % Timeslab at time=T and use dx_jn since bottom 
            [x_jn,x_j1n,x_jn1,x_j1n1,dx_jn,dx_jn1] = stvertices(j,T-dt);
            [x0,x1] = G_K_n(j,T-dt);
            err = err + quadrature1D(@(xi1) approx_err(x1(1,xi1))*(0.5*dx_jn1));
        end
        sqrt(err);
    end

% WRITE OUT EXPLICITLY WHAT SOURCE TERM IS AND FOR FULL CASE

% MAIN PROGRAM
% Project initial conditions onto initial degrees of freedom and solve each
% time slab sequentially using the previous degrees of freedom as initial
% conditions for the current slab
U = initialconditions(N);
reconstruct(U,t,dt,flowdomain,uexact);
while t < T
    Unext = globalassembly(U,t);
    t = t+dt;
    if plot; reconstruct(Unext,t,dt,flowdomain,uexact); end;
    U = Unext;
end
reconstruct(U,t,dt,flowdomain,uexact);
err = finalerror(U);
end