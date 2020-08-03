% Filename: phasePlane.m
% Author: Maud El-Hachem
% Queensland University of Technology, Brisbane, Australia, July 2020
% Reference:  M. El-Hachem, S.W. McCue, M.J. Simpson (2020) 
% Invading and receding sharp-fronted travelling waves.
% The script contains:
%   - two calls to the function plotPhasePlane to generate 
%     Figures 3(a) and (e).
%	- the function plotPhasePlane
%   - the function HeunsSolver

% Generating Figure 3(a)
% Calling the function plotPhasePlane(c,dz,z_begin,z_end) with the
% following inputs: c = 0.25, dz = 0.0001, z_begin = 0, z_end = 1000

dz = 0.0001;
z2 = 100;
z1 = 0;
c = 0.25;
[intersect]= plotPhasePlane(c,dz,z1,z2);

kappa= -c/intersect;
disp('c=')
disp(c)
disp('kappa=')
disp(kappa)

% Generating Figure 3(e)
% Calling the function plotPhasePlane(c,dz,z_begin,z_end) with the
% following inputs: c = 0.25, dz = 0.0001, z_begin = 0, z_end = 500

dz = 0.0001;
z2 = 100;
z1 = 0;
c = -0.5;
[intersect]= plotPhasePlane(c,dz,z1,z2);

kappa= -c/intersect;
disp('c = ')
disp(c)
disp('kappa = ')
disp(kappa)

% Function plotPhasePlane
% This function solves Equations (11) and (12) in the phase plane
% by Heun's method and plots the solution on the plane V(z) versus U(z).
% The same plot shows also the equilibrium points (0,0) and (1,0) 
% and the intersection point of the solution with the 
% corresponding axis U(z)=0.
% INPUT ARGUMENTS:
% ** c, the wave speed, if positive, 0 < c < 2, if negative, c < 0 
% ** dz, the step size used to discretise the domain of z.
% ** z_begin and z_end, the lower and upper limit of the domain of z, used
% to integrate numerically Equations (11) and (12) by Heun's method,
% such as z_begin <= z <= z_end. The initial conditions are applied at
% z = z_begin.
function [intersect]= plotPhasePlane(c,dz,z_begin,z_end)
    
    % for displaying purposes
    % color orange
    colors = [217 118 0]/255;
    figure
    % for displaying purposes
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(gca,'fontsize', 18);
    
    % depending on c, 
    % - calculating the eigenvalues and the eigenvectors of the solution 
    %   around the equilibrium point (1,0)
    % - determining the initial conditions close to the equilibrium point
    %   (1,0) along the eigenvector of the solution
    
    Us = 1;
    A = [0 1;(2*Us-1) -c];

    [v,d]=eig(A);
    IC = [0;0];
    if (c>0)
        if (d(1,1) > 0)
            IC = v(:,1);
        end
        if (d(2,2) > 0)
            IC = v(:,2);
        end
    end
    if (c<0)
        if (d(1,1) < 0)
            IC = v(:,1);
        end
        if (d(2,2) < 0)
            IC = v(:,2);
        end
    end
    
    % Setting the initial conditions
    IC2 = IC(1,1)*0.001;
    IC1 = Us+IC(2,1)*0.001;

    % Solving Equations Equations (11) and (12) in the phase plane with
    % Heun's method
    [U, V] = heunSolver(c, dz, z_begin, z_end, IC1, IC2);
     
    % Finding the intersection point with the axis U(z)=0
    intersects = [];
    for ii = 1:length(U)-1
        if (U(ii,1)>= 0 && U(ii+1,1)<= 0 || U(ii,1)<= 0 && U(ii+1,1)>= 0 )
            var1 = [U(ii,1) U(ii+1,1)];
            var2 = [V(ii,1) V(ii+1,1)];
            intersects = [ intersects interp1(var1,var2,0,'linear')];
        end
    end

    intersect = min(intersects);
    
    % for displaying purposes
    % computing the inferior limit of the V axis to be displayed
    if (intersect < -0.5)
        Vlimit = -0.05 + intersect;
    else
        Vlimit = -0.5;
    end
    
    
    % setting and computing the field vectors of the solution
    y1 = linspace(-0.4,0.9,8);
    y2 = linspace(Vlimit,0.4,8);
    [x,y] = meshgrid(y1,y2);
    du = zeros(size(x));
    dv = zeros(size(x));
    if c < 0
        f1 = @(Y) [-Y(2);-(-c*Y(2)+(-Y(1)+Y(1).*Y(1)));];
    else
        f1 = @(Y) [Y(2);(-c*Y(2)+(-Y(1)+Y(1).*Y(1)));];
    end
    for i = 1:numel(x)
        Yprime = f1([x(i); y(i)]);
        du(i) = Yprime(1);
        dv(i) = Yprime(2);
    end

    % plotting the axis at U(z) and V(z)
    line([-0.5 1],[0 0],'Color','k','LineStyle','-','LineWidth',2);
    hold on
    line([0 0],[Vlimit 0.5],'Color','k','LineStyle','-','LineWidth',2);
    % plotting the solution 
    plot(U,V,'r:','LineWidth',2,'Color',colors(1,1:3));

    % plotting the field vectors
    scaleFactor = 0.15;
    if (c<-2)
        scaleFactor = 0.1/(abs(c));
    end
    quiver(x,y,du*scaleFactor,dv*scaleFactor,'b','LineWidth',1,'AutoScale','off'); 
    
    % plotting the intersection point with U(z) = 0
    plot(0,intersect,'mo','MarkerEdgeColor','m','MarkerFaceColor','m','LineWidth',4);
    % plotting the equilibrium points (0,0) and (1,0)
    plot(0,0,'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',4);
    plot(1,0,'o','MarkerEdgeColor','k','MarkerFaceColor','k','LineWidth',4);
    % plotting the value of the wave speed c
    textc = strcat(strcat('$c = ',num2str(round(c,2))), '$');
    text(0.4,0.45,textc,'interpreter','latex','fontsize',18)
    % plotting the labels of the corresponding axis
    ylabel('$V(z)$','interpreter','latex');
    xlabel('$U(z)$','interpreter','latex'); 
    % limiting the axis of U(z) and V(z)
    xlim([-0.5,1]);
    ylim([Vlimit,0.5]);
    box on;
    hold off
end


% Function heunSolver
% This function solves Equations (11) and (12) by Heun's method 
% INPUT ARGUMENTS:
% ** c, the wave speed, if positive, 0 < c < 2, if negative, c < 0 
% ** dz, the step size used to discretise the domain of z
% ** z_begin and z_end, the lower and upper limit of the numerical domain 
% of z such as z_begin <= z <= z_end. The initial conditions are applied at
% z = z_begin.
% ** V1, X1, the values of the initial conditions
% OUTPUT ARGUMENTS:
% ** Uout : The solution U(z)
% ** Vout : The solution V(z)
% The output size of Vout and Xout may not be equal to the original number
% of nodes correponding to (z_end-z_begin)/dz+1 in the discretised domain.
% The array Uout and Vout may be truncated. 

function [Uout, Vout] = heunSolver(c,dz,z_begin,z_end,U1,V1)

    % z domain
    z = z_begin:dz:z_end;
    % number of nodes in the domain
    sz = length(z);

    % initialisation 
    V = zeros(sz,1);
    U = zeros(sz,1);
    U(1) = U1;
    V(1) = V1;
   
    szout = sz;
    
    % for all steps in the domain
    for i = 1:sz-1
        Ubar = dz * V(i) + U(i); 
        Vbar = dz * (-c*V(i)- U(i)*(1-U(i))) + V(i);
        U(i+1) = dz/2 * (V(i)+Vbar) + U(i); 
        V(i+1) = dz/2 * ((-c*V(i) - U(i)*(1-U(i))) + (-c*Vbar - Ubar *(1-Ubar))) + V(i); 
        % if the solution value is too large, stop the solver
        if (V(i+1) < c-1 || V(i+1) > 2 || U(i+1) < -2 || U(i+1) > 2)
            szout = i;
            break;
        end
    end

    Uout = U(1:szout,1);
    Vout = V(1:szout,1);
end