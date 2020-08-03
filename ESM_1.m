% Filename: PDE_solver.m
% Author: Maud El-Hachem
% Queensland University of Technology, Brisbane, Australia, July 2020
% Reference:  M. El-Hachem, S.W. McCue, M.J. Simpson (2020) 
% Invading and receding sharp-fronted travelling waves.
% The script contains:
%   - two calls to the function fisherStefanInvadingReceding to generate 
%     Figures 2(a) and (e).
%	- the function fisherStefanInvadingReceding
%   - the function tridia
% This function generates the solutions to Equations (1)-(3) with the 
% following parameters, and plot the density profiles at  t=0,10,20,30
% L(0) = 200, alpha = 0.5, kappa = {0.5859, -0.5387}
% dt = 0.1, dxi = 0.0001, discretisation steps used for |c| <<< 1
% if |c| >>> 1, as in figure 2 (d) and (h), use dt = 0.01, dxi = 0.000001

% Generating Figure 2(b)
% Calling the function fisherStefanInvadingReceding(kappa) with the
% following input: kappa = 0.5859
fisherStefanInvadingReceding(0.5859);
% setting the view of the axis of the figure from 150 to 250
xlim([150 250])
% Generating Figure 2(e)
% Calling the function fisherStefanInvadingReceding(kappa) with the
% following input: kappa = -0.5387
fisherStefanInvadingReceding(-0.5387);
% setting the view of the axis of the figure from 0 to 200
xlim([0 200])

function fisherStefanInvadingReceding(kappa)

    % for displaying purposes
    % array of colors used to print density profiles at required times
    % purple, green, orange
    colors = [126 47 142; 55 120 0;217 118 0]/255;
    % if kappa is positive, print in purple,
    % otherwise, print in green
    if (kappa >= 0)
        colorLine = colors(1,1:3);
    else
        colorLine = colors(2,1:3);
    end
    
    % \Delta t, Equation (4)
    dt = 0.1;
    % \Delta \xi, Equation (4)
    dxi = 0.0001;
    % spatial domains dicretised
    xi = 0:dxi:1;
    % number of nodes in the spatial domain
    nodes_xi = size(xi,2);
    % tolerance \epsilon used in Newton-Raphson algorithm
    tol = 1e-08;
    % total time
    total_time = 30;
    % total number of steps
    ts = round(total_time/dt+1);
    % times when is it required to print the density profiles
    time_toPrint = [0 10 20 30];
    % initial total length of the domain L(0) in Eq.(7)
    lt = 200;
    % parameter \alpha in the initial conditions in Eq.(7)
    alpha = 0.5;
    
    % previous position of the moving boundary
    lt_p = lt;
    % array of positions of the moving boundary for all times steps
    lt_array = zeros(1,ts); 
    % initialisation of variables used in Newton-Raphson algorithm
    % function F
    Fu = zeros(1,nodes_xi);
    % current densities
    u = zeros(1,nodes_xi);
    % coefficients a b c of the tridiagonal matrix
    % Jacobian
    % J(u)
    coeffA_u = zeros(1,nodes_xi);
    coeffB_u = zeros(1,nodes_xi);
    coeffC_u = zeros(1,nodes_xi);

    % Initialisation of density u(xi,0)
    for i = 1:nodes_xi
        u(1,i) = alpha;
    end
    % initialisation at the nodes where x=L(0)
    u(1,nodes_xi) = 0;
    % previous densities
    u_p = u;

    figure
    hold on
    % print the initial conditions
    if (time_toPrint(1)==0)
        plot(xi(1:nodes_xi)*lt, u, 'Color', colors(3,1:3), 'LineWidth',2 ,'DisplayName', '$u(x,0)$');
    end

    % for displaying purposes
    % formatting the text in latex font, size 18
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 18);
    
    %% Newton-Raphson algorithm - main loop
    % for all time steps
    for j = 1:ts
        
        %current time
        t = j * dt;
        disp(t);
        
        % storing the current position of the moving boundary
        lt_array(j) = lt;

        condition = 1;

        % while the tolerance is not reached
        while (condition)
            % boundary condition at x = 0, Eq.(5)
            coeffA_u(1,1) = 0.0;
            coeffB_u(1,1) = -1.0;
            coeffC_u(1,1) = 1.0;
            Fu(1) = -1.0*(u(1,2)-u(1,1));

            % boundary condition at x = L(t), Eq.(5)
            coeffA_u(1,nodes_xi) = 0;
            coeffB_u(1,nodes_xi) = 1.0;
            coeffC_u(1,nodes_xi) = 0;
            Fu(1,nodes_xi) = -1.0 * (u(1,nodes_xi));

            % J(u) delta u = -F(u)
            % Equation (4)
            for i = 2:nodes_xi-1
               coeffA_u(1,i) = 1.0/(dxi^2*lt^2) - (i-1)*dxi/lt * (lt-lt_p)/(2*dt*dxi);
               coeffB_u(1,i) = - 2.0/(dxi^2*lt^2)  - 1.0/dt + (1 - 2.0*u(1,i));
               coeffC_u(1,i) = 1.0/(dxi^2*lt^2) + (i-1)*dxi/lt * (lt-lt_p)/(2*dt*dxi);
               Fu(1,i) = -(u(1,i+1) - 2*u(1,i) + u(1,i-1))/(dxi^2*lt^2) ...
                   - (i-1)*dxi/lt * (u(1,i+1) - u(1,i-1)) * (lt-lt_p)/(2*dt*dxi) ...
                   + (u(1,i)-u_p(1,i)) / dt - u(1,i)*(1-u(1,i));
            end 
            delta_u = tridia(coeffA_u, coeffB_u, coeffC_u, Fu, nodes_xi);

            % correction of u
            for i = 1:nodes_xi
                u(1,i) = u(1,i) + delta_u(1,i);
            end

            if (norm(delta_u,Inf) <= tol)
                condition = 0;
            end
            %Stefan condition Eq.(6)
            lt = lt_p + dt*(-kappa*(-u(1,nodes_xi-1))/(dxi*lt));
        end

        % print density profile u(x,t) at each required time t
        if (isempty(find(time_toPrint == t,1)) == 0)
            plot(xi(1:nodes_xi)*lt, u, 'k-','Color',colorLine,'LineWidth',2);
        end
        % updating current L(t)
        lt_p = lt;
        % updating current u(x,t)
        u_p = u;
    end
    %%
    % computing final wave speed
    c = (lt_array(end) - lt_array(end-1))/dt;
    textc = strcat(strcat('$c = ',num2str(round(c,2))), '$');
    % printing final wave speed
    text(lt-50,0.4,textc,'interpreter','latex','fontsize',18)
    % printing axis labels
    ylabel('$u(x,t)$','interpreter','latex','fontsize',18);
    xlabel('$x$','interpreter','latex','fontsize',18);
    box on
    hold off
end
%% Function tridia
% This function implements Thomas algorithm that solves a system
% of equations Ax = d, where A is a tridiagonal matrix. The parameters 
% a,b and c are the three diagonals of the matrix A. N is the size of 
% the vector solution x.
function x = tridia(a,b,c,d,N)
    x=zeros(1,N);
    bb=b;
    dd=d;
    for i=2:N
        ff=a(i)/bb(i-1);
        bb(i)=bb(i)-c(i-1)*ff;
        dd(i)=dd(i)-dd(i-1)*ff;
    end

    for i=1:N-1
    x(N)=dd(N)/bb(N);    
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
    end
end 
