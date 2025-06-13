%% COMPUTATIONAL ENGINEERING
% VÍCTOR SAVALL DE RAMON
% CSM Project - Wing modellingsdfasdfasdf


%% PART I - BEAM MODELLING

% Useful commands to initialize script
clear
close all

%% 1) INPUT DATA

% Define the problem's data (e.g. dimensions, material properties, etc.)
% 1.1 Cross-section data and materials. For each different cross-section (m):

E   = 68.8e9;      % Young's modulus [Pa]
nu  = 0.33;        % Poisson ratio [-]
G   = E/(2*(1+nu)); % Shear modulus [Pa]
rho = 2700;        % Density [kg/m^3]
j_vec = [0;1;0];      % Orientation of y' axis
A   = 0.022;       % Area [m^2]
J   = 13.184e-4;   % Polar inertia J [m^4]
I_y = 1.097e-4;    % Inertia I_y' [m^4]
I_z = 12.087e-4;   % Inertia I_z' [m^4]
k_y = 0.5975;      % Shear correction factor k_y'
k_z = 0.1375;      % Shear correction factor k_z'
k_t = 0.2564;      % Torsion correction factor k_t

%Others
b = 5;              % Wingspan [m]

%% PREPROCESS

% 1.2 Mesh (discretization) data
load('beam.mat','xn','Tn','Tm');
N_elem = size(Tn,1);
N_nod = size(xn,1);
N_dof = 6*N_nod; % total number degrees of freedom
% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn : Nodal connectivities [N_elem x 2]
% Tm : Material connectivities [N_elem x 1]

% Define boundary conditions and forces data matrices: Up, Fe, Qe, Be

%Up: 1.3 Boundary conditions: - Fix nodes and DOFs matrix
Up = [0, 1, 1;  
      0, 1, 2;
      0, 1, 3;
      0, 1, 4;
      0, 1, 5;
      0, 1, 6]; %Restricts the first node in all DOFs

%Fe: 1.4 External forces: - Non-null point forces matrix (N)
%Qe: 1.4 External forces: - Distributed loads matrix (N/m)
%Be: 1.4 External forces: - Body forces matrix (N/kg)

% Load Type Setup
case_load = 'uForce';  % Options: 'uForce', 'uTorque'

Fe = zeros(0,3);
Qe = zeros(0,3);
Be = zeros(0,3); 

% Define Fe based on case
switch case_load
    case 'uForce'
        Fe = [1, N_elem, 3];

    case 'uTorque'
        Fe = [1, N_elem, 4];

    otherwise
        error('Unsupported load type: %s', case_load);
end


%% SOLVER
% TIP: To avoid recomputing the system matrices use a save/load structure:

if 1 % Set to 1 to (re)compute the system matrices and '0' to load them

    % 2) Assembly of global matrices
    [K,M,N,R,Me,Ba,Bs,Bt,Bb,Ka,Kb,Ks,Kt,le] = AssyMatrixBeams(Tm,Tn,xn,E,rho,I_y,I_z,J,G,A,k_y,k_z,k_t,j_vec);
    
    % Save the matrices to avoid recomputation
    save('beam_matrices.mat','K','M','N','R','Me','Ba','Bs','Bt','Bb','Ka','Kb','Ks','Kt','le');

else
    % Load previously computed results
    load('beam_matrices.mat','K','M','N','R','Me','Ba','Bs','Bt','Bb','Ka','Kb','Ks','Kt','le');

end

% 3) Compute global force vector:
% 3.1 Point loads:
f_vec = zeros(N_dof,1);
for q = 1:size(Fe,1)
    f_vec(6*(Fe(q,2)-1)+Fe(q,3),1) =  f_vec(6*(Fe(q,2)-1)+Fe(q,3),1) + Fe(q,1);
end


% 3.2 Nodal distributed forces:
Q = zeros(N_nod,6);
for r = 1:size(Qe,1)
    Q(Qe(r,2),Qe(r,3)) = Q(Qe(r,2),Qe(r,3)) + Qe(r,1);
end

% 3.3 Nodal body forces:
B = zeros(N_nod,6);
for s = 1:size(Be,1)
    B(Be(s,2),Be(s,3)) = B(Be(s,2),Be(s,3)) + Be(s,1);
end

%  3.4 Assembly process:
f_vec = AssyForceBeams(Tn,xn,R,N,Me,B,Q,f_vec);

% 4) Boundary conditions:
% 4.1 Initialization:
u_vec = zeros(N_dof,1);

% 4.2 Prescribed and free DOFs:
for p = 1:size(Up,1)
    Ip(p) = 6*(Up(p,2)-1) + Up(p,3);
    u_vec(Ip(p),1) = Up(p,1);
end
If = setdiff(1:N_dof,Ip);

% Perform modal analysis
N_modes = 6; % Number of modes
[Phi, lambda, freq] = ModalAnalysis(K, M, If, N_modes, N_dof);

% 5) Solve system of equations (static case):
% 5.1 Solve system:
u_vec(If,1) = [K(If,If)]\(f_vec(If,1)-[K(If,Ip)]*u_vec(Ip,1)); % displacements/rotations at free DOFs
fr = K*u_vec-f_vec; % reaction forces/moments at prescribed DOFs



%% 6) Postprocess: Computing strain and internal forces in beam elements

% 6.1 Get strain at each beam element:
for e = 1:N_elem
    % 6.1 a) Get element displacements:
    for j = 1:6
        Idof(j,1) = 6*(Tn(e,1)-1) + j;
        Idof(6+j,1) = 6*(Tn(e,2)-1) + j;
    end
    u_e = u_vec(Idof,1);

    % 6.1 b) Get each strain component:
    [strain_a, strain_s, strain_t, strain_b] = deal(zeros(1, N_elem));
    strain_a(1,e) = [Ba(1,:,e)]*[R(:,:,e)]*u_e;
    strain_s(:,e) = [Bs(1,:,e)]*[R(:,:,e)]*u_e;
    strain_t(1,e) = [Bt(1,:,e)]*[R(:,:,e)]*u_e;
    strain_b(:,e) = [Bb(1,:,e)]*[R(:,:,e)]*u_e;

    % 6.1.c) Get internal forces and moments at each element node:
    fint = zeros(12,N_elem);
    [Fx, Fy, Fz, Mx, My, Mz] = deal(zeros(2, N_elem));
    
    fint(:,e) = [R(:,:,e)]*(Ka(:,:,e) + Kb(:,:,e) + Ks(:,:,e) + Kt(:,:,e))*u_e;
    Fx(:,e) = [fint(1,e);-fint(7,e)];
    Fy(:,e) = [fint(2,e);-fint(8,e)];
    Fz(:,e) = [fint(3,e);-fint(9,e)];
    Mx(:,e) = [fint(4,e);-fint(10,e)];
    My(:,e) = [fint(5,e);-fint(11,e)];
    Mz(:,e) = [fint(6,e);-fint(12,e)];

end



% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)

% Extract u_y, u_z, and theta_x from vector u
idx = 1:6:6*N_nod;
u_y   = u_vec(idx + 1);
u_z   = u_vec(idx + 2);
theta_x = u_vec(idx + 3);


% Preallocate matrices
[phi_uy, phi_uz, phi_theta_x] = deal(zeros(N_nod, N_modes));

% Compute deflections and twist angle distributions for all modes
node_idx = 1:6:6*N_nod;
for i = 1:N_modes
    phi_uy(:, i) = Phi(node_idx + 1, i);
    phi_uz(:, i) = Phi(node_idx + 2, i);
    phi_theta_x(:, i) = Phi(node_idx + 3, i);
end
%save('beam_results.mat');

%Include plot functions

figure
hold on
plot(xn(:,1), u_z, 'LineWidth', 2, 'Color','r'); 
xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Vertical deflection, $u_z$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
title('Vertical deflection distribution along the wingspan', 'FontWeight', 'bold', 'FontSize', 14);
grid on;
grid minor;
box on;
xlim([0, b]);
hold off


figure
hold on
plot(xn(:,1),theta_x,'LineWidth', 2, 'Color','g');
xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Twist angle \theta_x (rad)')
title('Twist angle distribution along the wingspan', 'FontWeight', 'bold', 'FontSize', 14);

grid on
grid minor
box on
axis padded
xlim([0,b]);
hold off

figure
hold on
for i = 1:N_modes
    plot(xn(:,1),phi_uy(:,i),'LineWidth', 2)
end
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6',Location='best');
xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Horizontal deflection, $u_y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
grid on
grid minor
box on
axis padded
xlim([0,b]);
hold off

figure
hold on
for i = 1:N_modes
    plot(xn(:,1),phi_uz(:,i),'LineWidth', 2)
end
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6',Location='best')
xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Vertical deflection, $u_z$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
grid on
grid minor
box on
axis padded
xlim([0,b]);
hold off

figure
hold on
for i = 1:N_modes
    plot(xn(:,1),phi_theta_x(:,i),'LineWidth', 2)
end
legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6',Location='best')
xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Twist angle $\theta_x$ (rad)', 'Interpreter', 'latex', 'FontSize', 14)
grid on
grid minor
box on
axis padded
xlim([0,b]);
hold off


