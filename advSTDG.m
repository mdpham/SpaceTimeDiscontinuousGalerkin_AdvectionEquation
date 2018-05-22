% AMATH 990 A3: Space-time Discontinuous Galerkin for Advection Equation
% Usage:
%   N : number of space-time elements per timeslab
%   Nt : number of timeslabs
%   T : final time
%   flow : different flowdomains
%           'fixed'
%           'linear'
%           'full'
%   testcase : different manufactured solutions
%           'constant'
%           'linear'
%           'nonlinear1'
%           'nonlinear2'
%           'full'
%   draw: boolean for plotting
% The nice thing about DG is Geometric Conservation Law
% So try constant solution on moving mesh, should stay constant
function [U,approx_err] = advSTDG(N,Nt,T,flow,testcase,draw)
% FLOW DOMAIN
if strcmp(flow,'fixed')
    a = @(t) 0;
    b = @(t) 1;
    flowdomainaxis = [-0.5 1.5 -1 1.5];
elseif strcmp(flow,'linear')
    a = @(t) -0.2*t;
    b = @(t) 1-0.2*t;
    flowdomainaxis = [-0.4 1.4 -1 1.5];
elseif strcmp(flow,'full')
    a = @(t) sin(2*pi*t)/10;
    b = @(t) exp(-t);
    flowdomainaxis = [-0.2 1.2 -1.2 1.2];
else
    error('no such mesh')
end
% TEST CASES
% The source term is given by
% f = u_t + u_x, from the PDE
if strcmp(testcase,'constant') % CONSTANT SOLUTION
    uexact = @(t,x) 1;
    f = @(t,x) 0;
elseif strcmp(testcase,'linear') % LINEAR SOLUTION
    uexact = @(t,x) t + x;
    f = @(t,x) 2;
elseif strcmp(testcase,'nonlinear1') % NONLINEAR SOLUTION
    uexact = @(t,x) sin(2*pi*(x-t));
    f = @(t,x) 0;
elseif strcmp(testcase,'nonlinear2')
    uexact = @(t,x) sin(2*pi*x) + t;
    f = @(t,x) 2*pi*cos(2*pi*x) + 1;
elseif strcmp(testcase,'full')
    uexact = @(t,x) sin(2*pi*x)*sin(3*pi*t); % PROBLEM SOLUTION
    f = @(t,x) 3*pi*sin(2*pi*x)*cos(3*pi*t) + 2*pi*cos(2*pi*x)*sin(3*pi*t);    
else
    error('no such testcase');
end



% COMPUTATIONAL DOMAIN
% TIME VARIABLE
dt = T/Nt;
% MESH
X = zeros(Nt+1,N+1);
for i = 1:(Nt+1);
    t = (i-1)*dt;
    X(i,:) = linspace(a(t),b(t),N+1);
end
% GRID VELOCITY
% jnode the facet index (i.e. 1...N+1)
% t the time slab (i.e. 1...Nt)
    function vg = gridvel(t,jnode)
       dx = X(t+1,jnode)-X(t,jnode);
       vg = dx/dt;
    end
% REFERENCE MAPPING (Reference to Space-time for element j at slab t)
    function [x0,x1] = G_K_n(t,j)
        F1 = @(xi1) X(t,j)*0.5*(1-xi1) + X(t,j+1)*0.5*(1+xi1);
        F2 = @(xi1) X(t+1,j)*0.5*(1-xi1) + X(t+1,j+1)*0.5*(1+xi1);
        % Time
        x0 = @(xi0) (t-1)*dt + 0.5*dt*(1+xi0); 
        % Space        
        x1 = @(xi0,xi1) 0.5*(1-xi0)*F1(xi1) + 0.5*(1+xi0)*F2(xi1);
    end
% JACOBIAN
% Determinant of Jacobian matrix on t-th timeslab for space-time element j
    function J = jacobianfun(t,j)
        K_jn = X(t,j+1) - X(t,j);
        K_jn1 = X(t+1,j+1) - X(t+1,j+1);
        % For some reason there shouldn't be a 0.5
%         J = @(xi0) abs((0.5*dt)*0.25*((1-xi0)*K_jn + (1+xi0)*K_jn1));
        J = @(xi0) abs(dt*0.25*((1-xi0)*K_jn + (1+xi0)*K_jn1));
    end

% Function for getting global indices for jth element
j_idx = @(j) (j-1)*3 + (1:3);

% SIDE FACETS
% Right facet of last element
function sN1 = facet_sN1(t);
    vg = gridvel(t,N+1);
    [x0,x1] = G_K_n(t,N);
    sN1 = [
        quad1D(@(xi0) (1-vg)*(0.5*dt)) quad1D(@(xi0) (1-vg)*(0.5*dt)) quad1D(@(xi0) (1-vg)*xi0*(0.5*dt));
        quad1D(@(xi0) (1-vg)*(0.5*dt)) quad1D(@(xi0) (1-vg)*(0.5*dt)) quad1D(@(xi0) (1-vg)*xi0*(0.5*dt));
        quad1D(@(xi0) (1-vg)*xi0*(0.5*dt)) quad1D(@(xi0) (1-vg)*xi0*(0.5*dt)) quad1D(@(xi0) (1-vg)*xi0^2*(0.5*dt));
        ];
end
% Left facet of first element (BC imposed by solution)
function s1 = facet_s1(t);
    vg = gridvel(t,1);
    [x0,x1] = G_K_n(t,1);
    s1 = [
        quad1D(@(xi0) -(1-vg)*uexact(x0(xi0),x1(xi0,-1))*(0.5*dt));
        -quad1D(@(xi0) -(1-vg)*uexact(x0(xi0),x1(xi0,-1))*(0.5*dt));
        quad1D(@(xi0) -(1-vg)*uexact(x0(xi0),x1(xi0,-1))*xi0*(0.5*dt));
        ];
end
% Interior facets
% j the element
% Since we may be on vL or vR, need a switch for reference space
% lr = 0 if left, 1 if right
function sj1 = facet_sj1(t,j,lr);
    jnode = j+1;
    vg = gridvel(t,jnode);
    sj1 = [
        quad1D(@(xi0) (1-vg)*(0.5*dt)) quad1D(@(xi0) (1-vg)*(0.5*dt)) quad1D(@(xi0) (1-vg)*(0.5*dt)*xi0);
        quad1D(@(xi0) (1-vg)*(0.5*dt)) quad1D(@(xi0) (1-vg)*(0.5*dt)) quad1D(@(xi0) (1-vg)*(0.5*dt)*xi0);
        quad1D(@(xi0) (1-vg)*xi0*(0.5*dt)) quad1D(@(xi0) (1-vg)*xi0*(0.5*dt)) quad1D(@(xi0) (1-vg)*xi0^2*(0.5*dt));
        ];
    if lr
        sj1(2,:) = -sj1(2,:);
    end
%     disp(full(sj1));
end
% BOTTOM FACET
% Will be moved to the RHS since known information from previous timeslab
% is advected up and used as initial conditions for the current timeslab
function kjn = facet_kjn(t,j,u_prev)
    u_prev = u_prev(j_idx(j));
    u_jprev = @(xi1) u_prev(1) + u_prev(2)*xi1 + u_prev(3);
    kjn = -[
        quad1D(@(xi1) u_jprev(xi1)*(0.5*(X(t,j+1)-X(t,j))));
        quad1D(@(xi1) u_jprev(xi1)*xi1*(0.5*(X(t,j+1)-X(t,j))));
        -quad1D(@(xi1) u_jprev(xi1)*(0.5*(X(t,j+1)-X(t,j))));
        ];
%     disp(full(kjn));
end
% TOP FACET
function kjn1 = facet_kjn1(t,j)
    kjn1 = [
        quad1D(@(xi1) (0.5*(X(t+1,j+1)-X(t+1,j)))) ... % first row
        quad1D(@(xi1) xi1*(0.5*(X(t+1,j+1)-X(t+1,j)))) ...
        quad1D(@(xi1) (0.5*(X(t+1,j+1)-X(t+1,j))));
        quad1D(@(xi1) xi1*(0.5*(X(t+1,j+1)-X(t+1,j)))) ... % second row
        quad1D(@(xi1) xi1^2*(0.5*(X(t+1,j+1)-X(t+1,j)))) ...
        quad1D(@(xi1) xi1*(0.5*(X(t+1,j+1)-X(t+1,j))));
        quad1D(@(xi1) (0.5*(X(t+1,j+1)-X(t+1,j)))) ... % third row
        quad1D(@(xi1) xi1*(0.5*(X(t+1,j+1)-X(t+1,j)))) ...
        quad1D(@(xi1) (0.5*(X(t+1,j+1)-X(t+1,j))));
        ];
%     disp(full(kjn1));
end
% X DERIVATIVE VOLUME INTEGRAL
function Kvx = volume_K_x()
    Kvx = [ 0 0 0;
        quad2D(@(xi0,xi1) (0.5*dt)) ...
        quad2D(@(xi0,xi1) xi1*(0.5*dt)) ...
        quad2D(@(xi0,xi1) xi0*(0.5*dt));
        0 0 0;];
%     disp(full(Kvx));
end
% T DERIVATIVE VOLUME INTEGRAL
function Kvt = volume_K_t(t,j)
    Fn = @(xi1) 0.5*((1-xi1)*X(t,j) + (1+xi1)*X(t,j+1));
    Fn1 = @(xi1) 0.5*((1-xi1)*X(t+1,j) + (1+xi1)*X(t+1,j+1));
    Kjn = X(t,j+1) - X(t,j);
    Kjn1 = X(t+1,j+1) - X(t+1,j);
    Kvt = [ 0 0 0;
        -quad2D(@(xi0,xi1) 0.5*(Fn1(xi1)-Fn(xi1))) ...
        -quad2D(@(xi0,xi1) 0.5*(Fn1(xi1)-Fn(xi1))*xi1) ...
        -quad2D(@(xi0,xi1) 0.5*(Fn1(xi1)-Fn(xi1))*xi0);
        quad2D(@(xi0,xi1) 0.25*((1-xi0)*Kjn + (1+xi0)*Kjn1)) ...
        quad2D(@(xi0,xi1) 0.25*((1-xi0)*Kjn + (1+xi0)*Kjn1)*xi1) ...
        quad2D(@(xi0,xi1) 0.25*((1-xi0)*Kjn + (1+xi0)*Kjn1)*xi0); 
        ];
%     disp(full(Kvt));
end

% GLOBAL LINEAR SYSTEM
% U_prev the previous timeslab degrees of freedom used as IC for current
% t the timeslab index 
function U_next = globalsystem(U_prev,t)
    M = zeros(3*N,3*N); B = zeros(3*N,1);
    % BILINEAR FORM    
    for j = 1:N
        % VOLUME INTEGRALS
        M(j_idx(j),j_idx(j)) = M(j_idx(j),j_idx(j)) - (volume_K_x()+volume_K_t(t,j));
        % BOTTOM FACET: move to RHS
        B(j_idx(j)) = B(j_idx(j)) - facet_kjn(t,j,U_prev);
        % TOP FACET
        M(j_idx(j),j_idx(j)) = M(j_idx(j),j_idx(j)) + facet_kjn1(t,j);
        % SIDE FACETS
        if j == 1
            % Known information coming from left facet, move to RHS
            B(j_idx(1)) = B(j_idx(1)) - facet_s1(t);
            % Right facet treated as interior
            %  subtract for vR since jump is vL-vR            
            M(j_idx(1),j_idx(1)) = M(j_idx(1),j_idx(1)) + facet_sj1(t,1,0);
            M(j_idx(2),j_idx(1)) = M(j_idx(2),j_idx(1)) - facet_sj1(t,1,1);
        elseif j == N
            M(j_idx(N),j_idx(N)) = M(j_idx(j),j_idx(j)) + facet_sN1(t);
        else
            M(j_idx(j),j_idx(j)) = M(j_idx(j),j_idx(j)) + facet_sj1(t,j,0);
            M(j_idx(j+1),j_idx(j)) = M(j_idx(j+1),j_idx(j)) - facet_sj1(t,j,1);
        end
    end
    % LINEAR FORM
    S = zeros(3*N,1);   
    for j = 1:N
        [x0,x1] = G_K_n(t,j);
        J = jacobianfun(t,j);
        S(j_idx(j)) = [
            quad2D(@(xi0,xi1) f(x0(xi0),x1(xi0,xi1))*J(xi0));
            quad2D(@(xi0,xi1) f(x0(xi0),x1(xi0,xi1))*J(xi0)*xi1);
            quad2D(@(xi0,xi1) f(x0(xi0),x1(xi0,xi1))*J(xi0)*xi0);
            ];
    end
%     disp(full(M(1:6,1:6)))
%     disp(full(B))
    U_next = M\(B+S);
end

% INITIAL CONDITIONS
% Project initial conditions onto approximation space
function U0 = initialconditions(N)
    U0 = zeros(3*N,1);
    for j = 1:N
        [x0,x1] = G_K_n(1,j);
        J = jacobianfun(1,j);
        V = [
            quad2D(@(xi0,xi1) J(xi0)) quad2D(@(xi0,xi1) xi1*J(xi0)) quad2D(@(xi0,xi1) xi0*J(xi0))
            quad2D(@(xi0,xi1) xi1*J(xi0)) quad2D(@(xi0,xi1) xi1^2*J(xi0)) quad2D(@(xi0,xi1) xi0*xi1*J(xi0))
            quad2D(@(xi0,xi1) xi0*J(xi0)) quad2D(@(xi0,xi1) xi0*xi1*J(xi0)) quad2D(@(xi0,xi1) xi0^2*J(xi0))
            ];
        F = [
            quad2D(@(xi0,xi1) uexact(0,x1(-1,xi1))*J(xi0));
            quad2D(@(xi0,xi1) uexact(0,x1(-1,xi1))*xi1*J(xi0));
            quad2D(@(xi0,xi1) uexact(0,x1(-1,xi1))*xi0*J(xi0));
            ];
        U0(j_idx(j)) = V\F;
    end
end

% PLOTTING
% Plots degrees of freedom for t-th timeslab at top of timeslab
function [] = reconstruct(U,t)
    clf; xplot = []; yplot = [];
    % For each element plot the end points in order to
    % construct the global solution for plotting
    hold on
    for j = 1:N
        x1 = X(t+1,j); x2 = X(t+1,j+1); xplot = [x1 x2];
        % Approximate solution
        idx = j_idx(j);
        Uc = U(idx(1)); Us = U(idx(2)); Ut = U(idx(3));
        % by property of Legendre basis
        u1 = Uc - Us + Ut;
        u2 = Uc + Us + Ut;
        yplot = [u1 u2];
        plot(xplot,yplot,'r-','LineWidth',2);
        % Exact solution         
        x = linspace(x1,x2,10);
        y = uexact(t*dt,x);
        plot(x,y,'k-','LineWidth',1);
    end
    hold off
    % Plot options
    title(['elements = ',num2str(N), ', time = ', num2str(t*dt)]);
    xlabel('x'); ylabel('u');
    axis(flowdomainaxis);
    drawnow
end
% ERROR
% Compute error with exact solution by quadrature
function err = finalerror(u)
    err = 0;
    for j = 1:N
        [x0,x1] = G_K_n(Nt,j);
        u_el = u(j_idx(j));
        uc = u_el(1); us = u_el(2); ut = u_el(3);
        u_h = @(x) uc + us*x + ut; %xi0 = 1 at top of slab
        approx_err = @(x) abs(u_h(x) - uexact(T,x))^2;
        err = err + quad1D(@(xi1) approx_err(x1(1,xi1))*(0.5*(X(Nt+1,j+1)-X(Nt+1,j))));
    end
    sqrt(err);
end

% SOLVE TIMESLABS SEQUENTIALLY 
% % % % % % % % % % % % % % % 
U = initialconditions(N);
if draw; reconstruct(U,0); end;
for t = 1:(Nt)
    % Compute solution for t-th timeslab (i.e. 1...Nt)
    U_next = globalsystem(U,t);
    U = U_next;
    if draw; reconstruct(U,t); end;
end
% % % % % % % % % % % % % % % 
approx_err = finalerror(U);
end