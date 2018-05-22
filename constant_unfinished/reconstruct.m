% u is vector 2N
function [] = reconstruct(U,t,dt,flowdomain,uexact)
    N = length(U);
    
    yplot = [];
    xplot = [];
    
    x_n = flowdomain(t);
%     x_n1 = flowdomain(t+dt);
    % For each element plot the end points in order to
    % construct the global solution for plotting
    clf
    hold on
    for i = 1:N
        x1 = x_n(i);
        x2 = x_n(i+1);
        xplot = [x1 x2];
        % by property of Legendre polynomials
        y1 = U(i);
        y2 = U(i);
        yplot = [y1 y2];
        plot(xplot,yplot,'r-','LineWidth',2);
        title(['elements = ',num2str(N), ', time = ', num2str(t)]);
        xlabel('x');
        ylabel('u');
        %         
        x = linspace(x1,x2,10);
        y = uexact(t,x);
        plot(x,y,'r-','LineWidth',2);
    end
    hold off
    axis([-0.2 1 -1 3]);
    drawnow
end