%% Q1 b)

% Reyad ElMahdy
% 101064879

% In this part the code is modified so that the boundary conditions are V=0
% for the top and bottom boundaries and V=1 for the right and left sides,
% this is then compared with an analytical solution of the problem to see
% if they produce the same results and to find out the pros and cons of
% using each method

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
        if i == 1 || i == nx
            G(n,n) = 1;
            Bc(n) = 1;
        elseif j == 1 || j == ny
            G(n,n) = 1;
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

%%Analytic Solution

sol = zeros(nx,ny);
iMAX = 100;

[x,y] = meshgrid(linspace(-L/2,L/2,nx), linspace(0,W,ny));
a = W;
b = L/2;

for i = 1:iMAX
     n = 2*i-1;
     sol = sol+(4/pi).*(1/n).*(cosh((n*pi).*x'./a)./cosh((n*pi).*b./a)).*sin((n*pi).*y'./a);
     figure(3)
     surf(xaxis,yaxis,sol)
     pause(0.001)
end

% The analytical solution and the numerical solution both produce very
% similar results. The numerical solution is much faster, but it requires
% large matrix operations that make the iterations and mapping complex and harder to
% follow.
% The numerical solution would probably consume a lot of memory due
% to the large amount of data points (This is somwewhat mitigated when
% using sparse matrices but the effect still exists). 
% The analytical solution is simpler to follow and could be more accurate,
% but it is much slower than the numerical solution due to the time and
% number of iterations it takes to converge. The movie helped to watch the
% solution to see how fast it converges.