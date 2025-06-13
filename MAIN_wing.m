%% PART III - WING MODELLING (TEMPLATE)

% Useful commands to initialize script
clear
close all

%% DATA

% Define the problem's data (e.g. dimensions, material properties, etc.)
E   = 68.8e9;       % Young's modulus [Pa]
nu  = 0.33;         % Poisson ratio [-]
G   = E/(2*(1+nu)); % Shear modulus [Pa]
rho = 2700;         % Density [kg/m^3]
j_vec = [0;1;0];    % Orientation of y' axis
g   = 9.8;

h1  = 30e-3;        % Thickness (front) [m]
h2  = 25e-3;        % Thickness (aft) [m]
h3  = 11e-3;        % Thickness (upper) [m]
h   = [h1,h2,h3];   % Thickness (vector) [m]
y1  = 0.345;        % Length from leading edge to front of wingbox [m]
y2  = 0.960;        % Length from leading edge to aft of wingbox [m]
yc  = 0.5791;       % Shear center position [m]


b   = 5;               % Wingspan [m]
c   = 1.5;             % Chord [m]
h_s = 5e-3;            % Skin thickness [m]
h_r = 6e-3;            % Ribs thickness [m]
D   = 15e-3;           % Diameter (stringers) [m]
A   = pi*D^2;          % Area (stringers) [m^2]
I_y = 0.25*pi*(D/2)^4; % Area inertia y [m^4]
I_z = I_y;             % Area inertia z [m^4] (circular)
J   = I_y+I_z;         % Polar inertia  [m^4]
k_y_s = 5/6;           % Shear correction factor k_y' (stringer)
k_z_s = 5/6;           % Shear correction factor k_z' (stringer) 
k_t_s = 1;             % Torsion correction factor k_t (stringer)


%% PREPROCESS

% Load mesh data
load('wing.mat','xn','Tn_st','Tm_st','Tn_wb','Tm_wb','Tn_rb','Tm_rb','Tn_sk','Tm_sk','indRoot','indPointA','indPointB','indSpar1','indSpar2','n_u','n_l');
N_nod = size(xn,1);
N_dof = 6*N_nod; 

% xn : Nodal coordinates       
%      |x1 y1 z1|
%      |x2 y2 z2|
%      |  ...   |
% Tn_st : Nodal connectivities for beam elements (st: stringers) [Nelem x 2]
% Tn_wb, Tn_rb, Tn_sk : Nodal connectivities for shell elements (wb: wingbox, rb: ribs, sk: skin) [Nelem x 4]
% Tm_st, Tm_wb, Tm_rb, Tm_sk : Material connectivities for the different elements [Nelem x 1]
% indRoot   : Array of indices for root section nodes.
% indPointA : Index for node at point A.
% indPointB : Index for node at point B.
% indSpar1  : Array of indices for front spar centerline nodes.
% indSpar2  : Array of indices for rear spar centerline nodes.
% n_u, n_l  : Matrices containing information about the unit normals in the upper 
%             and lower surfaces, respectively. 
%       |nx1 ny1 nz1 id1|   First three columns are components of normal vector 
%       |nx2 ny2 nz2 id2|   Last (fourth) column is the nodal index
%       |      ...      |

% Define boundary conditions and forces data matrices: Up, Fe, Pe, Be

% 1.3 Boundary conditions Up (fix nodes and DOFs matrix)

Up = [zeros(n*6,1), repelem(indRoot,6,1), repmat((1:6)',n,1)]; %always fix nodes at root
Qe = zeros(0,3);

case_load = 'windTunnel'; % Options 'uForce', 'uTorque', 'tunnel'.
switch case_load
    case 'uForce' %F_z=1
        F_A = (y2-yc)/(y2-y1); 
        F_B = (yc-y1)/(y2-y1);
        Fe =[F_A, indPointA, 3; F_B, indPointB, 3];
        [Pe, Be] = deal(zeros(0,3));

    case 'uTorque'
        FA = -1/(y2-y1);
        FB = +1/(y2-y1);
        Fe =[F_A, indPointA, 3; F_B, indPointB, 3];
        [Pe, Be] = deal(zeros(0,3));

    case 'tunnel'
        p_inf    = 2.5e5;     % [Pa]
        alpha   = 9*2*pi/360; % [rad]

        for i = 1:length(n_u)
            x_u(i,:) = xn(n_u(i,4),:);
        end
        for i = 1:length(n_l)
            x_l(i,:) = xn(n_l(i,4),:);
        end
        [p_u, p_l] = deal(zeros(length(x_u),1), zeros(length(x_l),1));

        for i = 1:length(x_u)
            if x_u(i,1) <= 0.25*b
                P_x = p_inf*(0.5+4.8*x_u(i,1)/b-11.2*(x_u(i,1)/b)^2);
            elseif x_u(i,1) > 0.25*b && x_u(i,1) <= 0.75*b
                P_x = p_inf*(1.2-0.8*x_u(i,1)/b);
            else
                P_x = p_inf*(-2.4+8.8*x_u(i,1)/b-6.4*(x_u(i,1)/b)^2);
            end

            for j = 1:3
                p_u(i,j) = alpha*P_x*((1-x_u(i,2)/c)^4 + sqrt(1-x_u(i,2)/c))*n_u(i,j);
            end
        end
        for i = 1:length(x_l)
            if x_l(i,1) <= 0.25*b
                P_x = p_inf*(0.5+4.8*x_l(i,1)/b-11.2*(x_l(i,1)/b)^2);
            elseif x_l(i,1) > 0.25*b && x_l(i,1) <= 0.75*b
                P_x = p_inf*(1.2-0.8*x_l(i,1)/b);
            else
                P_x = p_inf*(-2.4+8.8*x_l(i,1)/b-6.4*(x_l(i,1)/b)^2);
            end
            for j = 1:3
                p_l(i,j) = -alpha*P_x*((1-x_l(i,2)/c)^4 -(1/4)*sqrt(1-x_l(i,2)/c))*n_l(i,j);
            end
        end

       n1 = length(x_u);
       n2 = length(x_l);
       Pe_vec = [
         reshape(p_u.', [], 1), repelem(n_u(:,4), 3), repmat((1:3).', n1, 1);
         reshape(p_l.', [], 1), repelem(n_l(:,4), 3), repmat((1:3).', n2, 1)];

        Be = zeros(N_nod,3);
        Be(:,1) = -g;
        Be(:,3) = 3;
        Be(:,2) = 1:N_nod;
        Fe = zeros(0,3);
    otherwise
        error('Unsupported load type: %s', case_load);
end


%% SOLVER

% Obtain system matrices

% TIP: To avoid recomputing the system matrices use a save/load structure:
if 1 % Set to 1 to (re)compute the system matrices and '0' to load them
    
    % Compute system matrices (as long as parameters don't change there is 
    % no need to repeat the matrix assembly on every run)
    % ...

   
    
% Stringers (with beam elements)
[K_st,M_st,N_st,R_st,Me_st,Ba_st,Bs_st,Bt_st,Bb_st,Ka_st,Kb_st,Ks_st,Kt_st] = AssyMatrixBeams(Tm_st,Tn_st,xn,E,rho,Iy,Iz,J,G,A,ky_st,kz_st,kt_st);

% Wingbox (with shell elements)
[K_wb,M_wb,N_wb,R_wb,Me_wb,S4_wb,Bb_wb,Bmn_wb,Bmt_wb,Bs_wb] = AssyMatrixShells(Tm_wb,Tn_wb,Tn_st,xn,h_wb,E,nu,rho);

% Skin (with shell elements)
[K_sk,M_sk,N_sk,R_sk,Me_sk,S4_sk,Bb_sk,Bmn_sk,Bmt_sk,Bs_sk] = AssyMatrixShells(Tm_sk,Tn_sk,Tn_st,xn,h_sk,E,nu,rho);

% Ribs (with shell elements)
[K_rb,M_rb,N_rb,R_rb,Me_rb,S4_rb,Bb_rb,Bmn_rb,Bmt_rb,Bs_rb] = AssyMatrixShells(Tm_rb,Tn_rb,Tn_st,xn,h_rb,E,nu,rho);


% K,m matrices

onoff = struct('wb',1, 'sk',1, 'rb',1, 'st',1);
K = onoff.wb*K_wb + onoff.sk*K_sk + onoff.rb*K_rb + onoff.st*K_st;
M = onoff.wb*M_wb + onoff.sk*M_sk + onoff.rb*M_rb + onoff.st*M_st;

    % Once (re)computed, save them to a separate data file
    save('wing_matrices.mat','K','M'); 
    % TIP: Add other potential results that can be reused in other parts
    % (e.g. element's length 'l', elements rotations matrices 'R', etc.)
    
else
    
    % Load previously computed results
    load('wing_matrices.mat','K','M');
    
end

% Compute external forces vector

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
f_vec = AssyForceBeams(Tn_st,xn,R_st,N_st,Me_st,B,Q,f_vec);
f_vec = AssyForceShells(Tn_wb,R_wb,N_wb,Me_wb,S4_wb,B,P,f_vec);
f_vec = AssyForceShells(Tn_rb,R_rb,N_rb,Me_rb,S4_rb,B,P,f_vec);
f_vec = AssyForceShells(Tn_sk,R_sk,N_sk,Me_sk,S4_sk,B,P,f_vec);

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

% 7) Compute stresses
[sigma_VM_wb] = AssyStressShells(Tn_wb,R_wb,u,Bb_wb,Bmn_wb,Bmt_wb,Bs_wb,nu,E,h_wb,Tm_wb);
[sigma_VM_sk] = AssyStressShells(Tn_sk,R_sk,u,Bb_sk,Bmn_sk,Bmt_sk,Bs_sk,nu,E,h_sk,Tm_sk);
[sigma_VM_rb] = AssyStressShells(Tn_rb,R_rb,u,Bb_rb,Bmn_rb,Bmt_rb,Bs_rb,nu,E,h_rb,Tm_rb);

%%%%%%%%%%--------FALTA MODEL ORDER REDUCTION-------------%%%%%%%%%%%%%

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

%% POSTPROCESS

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
save('wing_results.mat');

% Include plot functions
% ...

% Additional plot functions useful to visualize 3D model and modes

plotDeformed('wing',xn,Tn_wb,u,scale,sigma_VM_wb); % For wingbox elements
plotDeformed('wing',xn,Tn_rb,u,scale,sigma_VM_rb); % For rib elements
plotDeformed('wing',xn,Tn_sk,u,scale,sigma_VM_sk); % For skin elements
% This function plots the deformed structure and Von Mises stress distribution: 
% xn : Nodal coordinates matrix [Nnodes x 3]
% Tn : Nodal connectivities matrix [Nelem x 4]
% u : Displacements vector obtained from the system solution. It is expected
%     to be given as a column vector with dimensions [Ndof x 1].
% scale : Scale factor to amplify the displacements (set to appropriate 
%         number to visualize the deformed structure properly).
% sigVM : Von Mises stress at each Gauss point. It is expected to be given as 
%         a matrix with dimensions [Nelem x Ngauss].

plotModes('wing',Phi,freq,imodes)
% This function plots the specified modes resulting from a modal analysis
% in sets of 9.
% Phi : Modal displacements matrix in which each column corresponds to the
%       mode shape of the corresponding mode. Expected dimensions [Ndof x Nmodes]
% freq : Natural frequencies array. Expected dimensions [Nmodes x 1]
% imodes : Array selecting which modes to plot. A maximum of 9 modes must
%          can be selected. Example: imodes = [1,2,3,4,5,6,7,8,9] will plot
%          the modes stored in the first 9 columns of Phi / imodes = [1,4,5,10] 
%          will plot modes in columns 1, 4, 5 and 10 of Phi. 