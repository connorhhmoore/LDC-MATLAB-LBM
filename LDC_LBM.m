%% LBM Lid Driven Cavity Benchmark for Re=1000
% Connor Moore, <connor.moore@ontariotechu.net>
clear; clc; close all;

% Set number of cells
Nx = 30;
Ny = Nx;

% Max number of iterations
max_iter = 5000;

% Initialize geometry flags
free_cell   = 0;
wall_cell   = 1;
lid_cell    = 2;

% Initialize cell matrix
geom = zeros(Nx,Ny);

% Mark wall and lid cells
geom(:,1)   = wall_cell;
geom(:,end) = wall_cell;
geom(end,:) = wall_cell;
geom(1,:)   = lid_cell;

% Set initial density and the lid velocity
rho_ini      = 1;
lid_velocity = 0.1;

% Initialize more parameters
Re  = 1e3;
nu  = (Ny-1)*lid_velocity/Re;
tau = (6*nu+1)/2;

% The directions are as indicated:
%    9  5  6
%     \ | /
%   4 - 1 - 2
%     / | \
%    8  3  7

% Initialize required matrices
F   = (zeros(Nx,Ny,9));  % Distribution function for each cell
Feq = (zeros(size(F)));  % Equilibirum distr. function
rho = (zeros(Nx,Ny));    % Macroscopic density
u   = (zeros(Nx,Ny));    % Macroscopic x-velocity
v   = (zeros(Nx,Ny));    % Macroscopic y-velocity
sqr = (zeros(Nx,Ny));    % u² + y²

% Set initial distributions:
F(:,:,1)         = rho_ini*(4/9);   % Center node
F(:,:,[2,3,4,5]) = rho_ini*(1/9);   % Side/Up-Down nodes
F(:,:,[6,7,8,9]) = rho_ini*(1/36);  % Corner nodes

% Initialize some value keys
u_key = [0,1,0,-1,0,1,1,-1,-1];
v_key = [0,0,-1,0,1,1,-1,-1,1];
w     = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36];

% Used for bounce back scheme
wall_key = [1,6,7,8,9,2,3,4,5];

% Initialize GIF file
figure(Name="LBM Lid-driven cavity")
outputFile = 'LBM_LDC.gif';

for i = 1:max_iter
    
    rho = sum(F,3);
    u = sum(F.*shiftdim(repmat(u_key,[Nx,1,Ny]),2),3)./rho;
    v = sum(F.*shiftdim(repmat(v_key,[Nx,1,Ny]),2),3)./rho;

    % Set known velocities for the lid
    u(geom==lid_cell) = lid_velocity;
    v(geom==lid_cell) = 0;

    % Define sum of square velocities
    sqr = u.^2 + v.^2;

    % Update equilibrium distribution for each ordinate
    Feq(:,:,1) = w(1)*rho.*(1 - 1.5.*sqr);
    Feq(:,:,2) = w(2)*rho.*(1 + 3*u + 4.5*u.^2 - 1.5*sqr);
    Feq(:,:,3) = w(3)*rho.*(1 - 3*v + 4.5*v.^2 - 1.5*sqr);
    Feq(:,:,4) = w(4)*rho.*(1 - 3*u + 4.5*u.^2 - 1.5*sqr);
    Feq(:,:,5) = w(5)*rho.*(1 + 3*v + 4.5*v.^2 - 1.5*sqr);
    Feq(:,:,6) = w(6)*rho.*(1 + 3*(u+v) + 4.5*(u+v).^2 - 1.5*sqr);
    Feq(:,:,7) = w(7)*rho.*(1 + 3*(u-v) + 4.5*(u-v).^2 - 1.5*sqr);
    Feq(:,:,8) = w(8)*rho.*(1 + 3*(-u-v) + 4.5*(-u-v).^2 - 1.5*sqr);
    Feq(:,:,9) = w(9)*rho.*(1 + 3*(-u+v) + 4.5*(-u+v).^2 - 1.5*sqr);

    % Apply bounceback at wall positions
    for j = find(geom==wall_cell)'
        [row,col] = ind2sub([Nx,Ny],j);
        F(row,col,[1,2,3,4,5,6,7,8,9])=F(row,col,[1,4,5,2,3,8,9,6,7]);
    end

    % Apply condition to lid cell
    for j = find(geom==lid_cell)'
        [row,col] = ind2sub([Nx,Ny],j);
        F(row,col,:) = Feq(row,col,:);
    end

    % Apply fluid cell calculations
    for j = find(geom==free_cell)'
        [row,col] = ind2sub([Nx,Ny],j);
        F(row,col,:)=F(row,col,:).*(1-1/tau)+Feq(row,col,:)./tau;
    end

    % Particle streaming for each ordinate
    F(:,:,2) = circshift(F(:,:,2),1,2);
    F(:,:,3) = circshift(F(:,:,3),-1,1);
    F(:,:,4) = circshift(F(:,:,4),-1,2);
    F(:,:,5) = circshift(F(:,:,5),1,1);
    F(:,:,6) = circshift(F(:,:,6),[1,1]);
    F(:,:,7) = circshift(F(:,:,7),[-1,1]);
    F(:,:,8) = circshift(F(:,:,8),[-1,-1]);
    F(:,:,9) = circshift(F(:,:,9),[1,-1]);

    fprintf("%i/%i\n",i,max_iter);

    % Capture every 100 iterations for plotting
    if(mod(i,100)==0)
        imagesc(sqrt(sqr)/lid_velocity);
        title("Lid-driven cavity LBM solution Relative veclocity u/u_{lid}, iter: "+i)
        colorbar;
        xlabel("{\itx}-position");
        ylabel("{\ity}-position");
        obj=gca;
        exportgraphics(obj, outputFile, Append=true);
    end
end