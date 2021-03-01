% Q2
% Reyad ElMahdy
% 101064879

%% Question 2
% In this part, the code was modified so that there are two resistive boxes
% inside the region forming a 'bottle's neck' for the current to pass
% through. The code uses the same boundary confitions as the previous part.
% After remapping the G matrix and generating the V plot, the current is
% then examined with varying box with, conductivity, and mesh spacing to
% see how it changes. 

close all
clear
clc

W = 20;
L = 1.5*W;

Lb = 9;
Wb = 6;
s = 0.5;
nx = floor(L/s + 1);
ny = floor(W/s + 1);

cond = 1; % Conductive
res = 1e-2; % Resistive

% Generating the conductivity map
condMap = zeros(nx,ny);
for i = 1:nx
   for j = 1:ny
       if (i-1>0.5*(L-Lb)/s) && ((i-1)<0.5*(L+Lb)/s) && (((j-1)<Wb/s||(j-1)>(W-Wb)/s))
           condMap(i,j) = res;
       else
           condMap(i,j) = cond;
       end
   end
end

% Adding the map to the plot
figure(1)
imagesc([0 L],[0 W],condMap');
title('Map of Electrical Conductivity')
% Yellow region is conductive, blue is resistive

G = sparse(nx*ny);
Bc = zeros(1,nx*ny);
for i = 1:nx
    for j = 1:ny
        n = j +(i-1)*ny; %Mapping Eq
        %Setting Boundary conditions
        if i == 1 || i == nx
            G(n,n) = 1;
            if i == 1
                Bc(n) = 1;
            end
        elseif j == 1 || j == ny
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nyp = j+1+(i-1)*ny;
            nym = j-1+(i-1)*ny;
            
            %Resistances from the conductivity map
            rxm = (condMap(i,j)+condMap(i-1,j))/2;
            rxp = (condMap(i,j)+condMap(i+1,j))/2;
            
            if j == 1
                ryp = (condMap(i,j)+condMap(i,j+1))/2;
                G(n,n) = -(rxm+rxp+ryp);
                G(n,nyp) = ryp; 
            else
                rym = (condMap(i,j)+condMap(i,j-1))/2;
                G(n,n) = -(rxm+rxp+rym);
                G(n,nym) = rym;
            end
            G(n,nxm) = rxm;
            G(n,nxp) = rxp; 
            
        % Rest of the nodes
        else
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nym = j-1+(i-1)*ny;
            nyp = j+1+(i-1)*ny;
            
            % Getting the resistances from the conduction map
            rxm = (condMap(i,j)+condMap(i-1,j))/2;
            rxp = (condMap(i,j)+condMap(i+1,j))/2;
            ryp = (condMap(i,j)+condMap(i,j+1))/2;
            rym = (condMap(i,j)+condMap(i,j-1))/2;
            
            % Assigning the node equations
            G(n,n) = -(rxm+rxp+rym+ryp);
            G(n,nxm) = rxm;
            G(n,nxp) = rxp;
            G(n,nym) = rym;
            G(n,nyp) = ryp; 
        end
    end
end

%solving and plotting
V = G\Bc';
mappedV = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j+(i-1)*ny;
        mappedV(i,j) = V(n);
    end
end

[xaxis, yaxis] = meshgrid(0:s:L,0:s:W);
figure(2)
surf(xaxis',yaxis',mappedV)
hold on
imagesc([0 L],[0 W],mappedV')
title('Potential Over the Region')
hold off

[Ey, Ex] = gradient(-mappedV);

figure(3)
quiver(xaxis',yaxis',Ex,Ey)
xlim([0 L])
ylim([0 W])
title('E Field Quiver')
Jx = condMap.*Ex;
Jy = condMap.*Ey;

figure(4)
quiver(xaxis',yaxis',Jx,Jy)
xlim([0 L])
ylim([0 W])
title('Current Density Quiver')

% Currents at the two contacts
I1 = sum(Jx(1,:))
I2 = sum(Jx(nx,:))


% Generating the other plots
s = linspace(0.5, 10, 20);
for i = 1:20
    I(i) = Curr(s(i), res, Wb);
end
figure(5)
plot(s,I)
title('Plot of current vs mesh density')

cond = linspace(0, 1, 20);
for i = 1:20
    I(i) = Curr(0.5, cond(i), Wb);
end
figure(6)
plot(cond,I)
title('Plot of current vs conductivty of the boxes')

Wb = linspace(0,12,20);
for i = 1:20
    I(i) = Curr(0.5, 0.01, Wb(i));
end
figure(7)
plot(Wb,I)
title('Plot of current vs box width')

% Creating a reuseable function that outputs the current
% This will be used for the remaining plots

function [I] = Curr(s, cond, Wb)
W = 20;
L = 1.5*W;
Lb = 9;
nx = floor(L/s+1);
ny = floor(W/s+1);

    condMap = zeros(nx,ny);
    for i = 1:nx
        for j = 1:ny
            if (i-1>0.5*(L-Lb)/s) && ((i-1)<0.5*(L+Lb)/s) && (((j-1)<Wb/s||(j-1)>(W-Wb)/s))
                condMap(i,j) = cond;
            else
                condMap(i,j) = 1;
            end
        end
    end
   
    G = sparse(nx*ny);
    Bc = zeros(1,nx*ny);
    
    for i = 1:nx
        for j = 1:ny
            n = j +(i-1)*ny; %Mapping Eq
            %Setting Boundary conditions
            if i == 1 || i == nx
                G(n,n) = 1;
                if i == 1
                    Bc(n) = 1;
                end
            elseif j == 1 || j == ny
                nxm = j+(i-2)*ny;
                nxp = j+i*ny;
                nyp = j+1+(i-1)*ny;
                nym = j-1+(i-1)*ny;
            
                % Resistances from the conductivity map
                rxm = (condMap(i,j)+condMap(i-1,j))/2;
                rxp = (condMap(i,j)+condMap(i+1,j))/2;
            
                if j == 1
                    ryp = (condMap(i,j)+condMap(i,j+1))/2;
                    G(n,n) = -(rxm+rxp+ryp);
                    G(n,nyp) = ryp; 
                else
                    rym = (condMap(i,j)+condMap(i,j-1))/2;
                    G(n,n) = -(rxm+rxp+rym);
                    G(n,nym) = rym;
                end
                
            G(n,nxm) = rxm;
            G(n,nxp) = rxp; 
            
            % Rest of the nodes
            else
                nxm = j+(i-2)*ny;
                nxp = j+i*ny;
                nym = j-1+(i-1)*ny;
                nyp = j+1+(i-1)*ny;
            
                % Getting the resistances from the conduction map
                rxm = (condMap(i,j)+condMap(i-1,j))/2;
                rxp = (condMap(i,j)+condMap(i+1,j))/2;
                ryp = (condMap(i,j)+condMap(i,j+1))/2;
                rym = (condMap(i,j)+condMap(i,j-1))/2;
            
                % Assigning the node equations
                G(n,n) = -(rxm+rxp+rym+ryp);
                G(n,nxm) = rxm;
                G(n,nxp) = rxp;
                G(n,nym) = rym;
                G(n,nyp) = ryp; 
            end
        end
    end
    
    V = G\Bc';
    mappedV = zeros(nx,ny);
    for i = 1:nx
        for j = 1:ny
            n = j+(i-1)*ny;
            mappedV(i,j) = V(n);
        end
    end
    
    [Ey, Ex] = gradient(-mappedV);
    Jx = condMap.*Ex;
    Jy = condMap.*Ey;
    
    I = sum(Jx(1,:));
end