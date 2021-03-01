% Q1 a)

% Reyad ElMahdy
% 101064879

close all
clear
clc

W = 20;
L = 1.5*W;

% Getting matrix dimensions
s = 0.5; % Space between each point in the meshgrid
nx = floor(L/s + 1);
ny = floor(W/s + 1);

%Creating matricies and mapping the 2D geometry
G = sparse(nx*ny);
Bc = zeros(1,nx*ny);

for i = 1:nx
    for j = 1:ny
        n = j+(i-1)*ny; % Mapping Equation
        
        % Setting boundary conditions
        if i == 1
            G(n,n) = 1;
            Bc(n) = 1;
        elseif i == nx
            G(n,n) = 1;
        elseif j == 1 || j == ny
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            if j == 1
                nyp = j+1+(i-1)*ny;
                G(n,nyp) = 1;
            else
                nym = j-1+(i-1)*ny;
                G(n,nym) = 1;
            end
            G(n,n) = -3;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
        else
            % Rest of the nodes
            nxm = j+(i-2)*ny;
            nxp = j+i*ny;
            nym = j-1+(i-1)*ny;
            nyp = j+1+(i-1)*ny;
            G(n,n) = -4;
            G(n,nxm) = 1;
            G(n,nxp) = 1;
            G(n,nym) = 1;
            G(n,nyp) = 1;
        end
    end
end

% Find V and map it to the space
V = G\Bc';
mappedV = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        n = j+(i-1)*ny;
        mappedV(i,j) = V(n);
    end
end

% Plot
xaxis = 0:s:L;
figure(1)
plot(xaxis,mappedV(:,floor(ny/2)))
figure(2)
[xaxis, yaxis] = meshgrid(0:s:W, 0:s:L);
surf(xaxis,yaxis,mappedV)


