%% PART II - SHELL MODELLING

% Useful commands to initialize script
clear
close all

%% DATA
% 1.1 Material properties and thickness for each different shell
% Define the problem's data (e.g. dimensions, material properties, etc.)
E   = 68.8e9;      % Young's modulus [Pa]
nu  = 0.33;        % Poisson ratio [-]
rho = 2700;        % Density [kg/m^3]
h1  = 30e-3;       % Thickness (front) [m]
h2  = 25e-3;       % Thickness (aft) [m]
h3  = 11e-3;        % Thickness (upper) [m]
h   = [h1,h2,h3];  % Thickness (vector) [m]
y1  = 0.345;       % Length from leading edge to front of wingbox [m]
y2  = 0.960;       % Length from leading edge to aft of wingbox [m]
yc  = 0.5791;      % Shear center position [m]

%Others
b = 5;              % Wingspan [m]


%% PREPROCESS
% 1.2 Mesh discretization data 
% Load mesh data
load('shell.mat','xn','Tn','Tm','indRoot','indPointA','indPointB','indSpar1','indSpar2');
N_elem = size(Tn,1);
N_nod = size(xn,1);
N_dof = 6*N_nod; % total number degrees of freedom
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [Nelem x 4]
% Tm : Material connectivities [Nelem x 1]
% indRoot   : Array of indices for root section nodes.
% indPointA : Index for node at point A.
% indPointB : Index for node at point B.
% indSpar1  : Array of indices for front spar centerline nodes.
% indSpar2  : Array of indices for rear spar centerline nodes.

% Define boundary conditions and forces data matrices: Up, Fe, Pe, Be
n = size(indRoot,1);
% 1.3 Boundary conditions
Up = [zeros(n*6,1), repelem(indRoot,6,1), repmat((1:6)',n,1)];
% 1.4 External forces
Pe = zeros(0,3);
Be = zeros(0,3);

case_load = 'uForce';  % Options: 'uForce', 'uTorque'
switch case_load
    case 'uForce' %F_z=1
        F_A = (y2-yc)/(y2-y1); 
        F_B = (yc-y1)/(y2-y1);
    case 'uTorque'
        F_A = -1/(y2-y1); %M_x=1
        F_B = 1/(y2-y1);
    otherwise
        error('Unsupported load type: %s', case_load);
end
Fe =[F_A, indPointA, 3; F_B, indPointB, 3];

%% SOLVER

% Obtain system matrices

% TIP: To avoid recomputing the system matrices use a save/load structure:
if 1 % Set to 1 to (re)compute the system matrices and '0' to load them
    
    % Compute system matrices (as long as parameters don't change there is 
    % no need to repeat the matrix assembly on every run)
  
    % 2) Assembly of global matrices:
    % 3) Compute artificial rotation stiffness matrix:
    [K,M,N,R,Me,S_4,Bb,Bmn,Bmt,Bs] = AssyMatrixShells(Tm,Tn,xn,h,E,nu,rho,0);

    % Once (re)computed, save them to a separate data file
    save('shell_matrices.mat','K','M'); 
    % TIP: Add other potential results that can be reused in other parts
    % (e.g. element's length 'l', elements rotations matrices 'R', etc.)
    
else
    
    % Load previously computed results
    load('shell_matrices.mat','K','M');
    
end


% 4) Compute global force vector
% 4.1 Point loads:
f_vec = zeros(N_dof,1);
for q = 1:size(Fe,1)
    f_vec(6*(Fe(q,2)-1)+Fe(q,3),1) = f_vec(6*(Fe(q,2)-1)+Fe(q,3),1) + Fe(q,1);
end %loop over rows in Fe

% 4.2 Nodal distributed forces
P = zeros(N_nod,6);
for r = 1:size(Pe,1)
    P(Pe(r,2),Pe(r,3)) = P(Pe(r,2),Pe(r,3)) + Pe(r,1);
end %loop over rows in Pe

% 4.3 Nodal body forces
B = zeros(N_nod,6);
for s = 1:size(Be,1)
    B(Be(s,2),Be(s,3)) = B(Be(s,2),Be(s,3)) + Be(s,1);
end %loop over rows in Be

% 4.4 Assembly process
f_vec = AssyForceShells(Tn,B,P,R,N,Me,S_4,f_vec);

% 5) Boundary conditions
% 5.1 Initialization
u_vec = zeros(N_dof,1);

% 5.) Prescribed and free DOFs
for p = 1:size(Up,1)
    Ip(p) = 6*(Up(p,2)-1) + Up(p,3);
    u_vec(Ip(p),1) = Up(p,1);
end %loop over rows in Up

If = setdiff(1:N_dof,Ip);
% Solve system

% Perform modal analysis
N_modes = 6; % Number of modes
[Phi, lambda, freq] = ModalAnalysis(K, M, If, N_modes, N_dof);

% 6) Solve system of equations (static case)
% 6.1 Solve system
u_vec(If,1) = [K(If,If)]\(f_vec(If,1)-[K(If,Ip)]*u_vec(Ip,1));
f_r = K*u_vec-f_vec;

% 7) Postprocess: Computing local strain and stress in shell elements:
sigma_VM = AssyStressShells(Tn,R,u_vec,Bb,Bmn,Bmt,Bs,nu,E,h,Tm);

%% POSTPROCESS

% Deflections and twist angle distribution
% Efficient extraction of DOFs
idxS_1 = (indSpar1 - 1) * 6;
idxS_2 = (indSpar2 - 1) * 6;

uS_1 = u_vec(idxS_1' + (1:6)');
uS_2 = u_vec(idxS_2' + (1:6)');

% Twist angle and displacement estimation
theta_x = (uS_2(3,:) - uS_1(3,:)) / (y2 - y1);
u_z     = uS_1(3,:) + theta_x .* (yc - y1);
u_y     = 0.5 * (uS_1(2,:) + uS_2(2,:));

% Mode shape extraction
id_xS1 = (indSpar1 - 1) * 6;
id_xS2 = (indSpar2 - 1) * 6;

uz1 = Phi(id_xS1 + 3, :);
uz2 = Phi(id_xS2 + 3, :);
uy1 = Phi(id_xS1 + 2, :);
uy2 = Phi(id_xS2 + 2, :);

theta_x_avg = (uz2-uz1)/(y2-y1);
u_z_avg     = uz1+theta_x_avg.*(yc-y1);
u_y_avg     = (uy1+uy2)/2;

xS = xn(indSpar1, 1);  

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
save('shell_results.mat');

% Include plot functions

figure
hold on
for mode = 1:N_modes
plot(xS,u_z_avg(:,mode),'LineWidth',2)
end
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6',Location='best')
xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Vertical deflection, $u_z$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
title('Vertical deflection distribution along the wingspan', 'FontWeight', 'bold', 'FontSize', 14);
grid on
grid minor
box on
axis padded
xlim([0,b]);
hold off


figure
hold on
for mode = 1:N_modes
plot(xS,u_y_avg(:,mode),'LineWidth',2)
end
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6',Location='best')
xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Horizontal deflection, $u_y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
title('Horizontal deflection distribution along the wingspan', 'FontWeight', 'bold', 'FontSize', 14);
grid on
grid minor
box on
axis padded
xlim([0,b]);
hold off




figure
hold on
for mode = 1:N_modes
plot(xS,theta_x_avg(:,mode),'LineWidth',2)
end
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6',Location='best')
xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Twist angle $\theta_x$ (rad)', 'Interpreter', 'latex', 'FontSize', 14)
title('Twist angle distribution along the wingspan', 'FontWeight', 'bold', 'FontSize', 14);
grid on
grid minor
box on
axis padded
xlim([0,b]);


figure
hold on
plot(xn(indSpar1,1),u_z,'LineWidth',2,'Color','r');
xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Vertical deflection, $u_z$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
title('Vertical deflection distribution along the wingspan', 'FontWeight', 'bold', 'FontSize', 14);
grid on
grid minor
box on
axis padded
xlim([0,b]);
hold off

figure
hold on
plot(xn(indSpar1,1),theta_x,'LineWidth', 2, 'Color','g');
xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Twist angle $\theta_x$ (rad)', 'Interpreter', 'latex', 'FontSize', 14)
title('Twist angle distribution along the wingspan', 'FontWeight', 'bold', 'FontSize', 14);
grid on
grid minor
box on
axis padded
xlim([0,b]);
hold off











% Additional plot functions useful to visualize 3D model and modes
scale=1e6;
plotDeformed('shell',xn,Tn,u_vec,scale);
% This function plots the deformed structure: 
% xn : Nodal coordinates matrix [Nnodes x 3]
% Tn : Nodal connectivities matrix [Nelem x 4]
% u : Displacements vector obtained from the system solution. It is expected
%     to be given as a column vector with dimensions [Ndof x 1].
% scale : Scale factor to amplify the displacements (set to appropriate 
%         number to visualize the deformed structure properly).

imodes = [1,2,3,4,5,6];
plotModes('shell',Phi,freq,imodes)
% This function plots the specified modes resulting from a modal analysis
% in sets of 9.
% Phi : Modal displacements matrix in which each column corresponds to the
%       mode shape of the corresponding mode. Expected dimensions [Ndof x Nmodes]
% freq : Natural frequencies array. Expected dimensions [Nmodes x 1]
% imodes : Array selecting which modes to plot. A maximum of 9 modes must
%          can be selected. Example: imodes = [1,2,3,4,5,6,7,8,9] will plot
%          the modes stored in the first 9 columns of Phi / imodes = [1,4,5,10] 
%          will plot modes in columns 1, 4, 5 and 10 of Phi. 