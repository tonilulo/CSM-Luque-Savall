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

case_load = 'uTorque';  % Options: 'uForce', 'uTorque'
switch case_load
    case 'uForce' %F_z=1
        F_A = -(y2-yc)/(y2-y1); 
        F_B = -(yc-y1)/(y2-y1);
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
  
    K = sparse(N_dof,N_dof);
M = sparse(N_dof,N_dof);

% 2.2) Assembly process:
for e = 1:N_elem

    % 2.2 a) Compute rotation matrix:
    S = cross((xn(Tn(e,3),:)'-xn(Tn(e,1),:)'),(xn(Tn(e,4),:)'-xn(Tn(e,2),:)'))/2;
    k_vec = S/norm(S); %mormal vector of the flat shell element
    d = (xn(Tn(e,2),:)' + xn(Tn(e,3),:)' - xn(Tn(e,4),:)' - xn(Tn(e,1),:)')/2;
    i_vec = d/norm(d);
    j_vec = cross(k_vec,i_vec);

    R_mat = [i_vec j_vec k_vec zeros(3,2); zeros(3,3) i_vec j_vec]';
    R(:,:,e) = blkdiag(R_mat, R_mat, R_mat, R_mat);

    % 2.2 b) Get nodal coefficients for the shape functions:
    a = [-1, 1, 1, -1];
    b = [-1, -1, 1, 1];

    % 2.2 c) Compute element matrices:
    % c1) 1 Gauss point quadrature matrices:
    N_1 = [1, 1, 1, 1]'/4;
    N_1ksi = a/4;
    N_1eta = b/4;
    J_1 = zeros (2,2);

    for i = 1:4
        J_1 = J_1 + [N_1ksi(i); N_1eta(i)]*xn(Tn(e,i),:)*[i_vec, j_vec];
    end %loop over nodes

    N1_xmat = J_1^(-1)*[N_1ksi; N_1eta];
    S_1 = 4*det(J_1); %area associated to Gauss point

    % c1.1) Shear component of stiffness matrix:
    Bs_i = zeros(2,5,4);
    for i = 1:4
        Bs_i(:,:,i) = [0, 0, N1_xmat(1,i), 0, N_1(i);
                       0, 0, N1_xmat(2,i), -N_1(i), 0];
    end %loop over nodes

    Cs = [1, 0; 
          0, 1]*5*h(Tm(e))*E/(12*(1 + nu)); %E,nu=ct.

    Bs(:,:,e) = [Bs_i(:,:,1),Bs_i(:,:,2),Bs_i(:,:,3),Bs_i(:,:,4)];

    Ks(:,:,e) = S_1*[R(:,:,e)]'*[Bs(:,:,e)]'*Cs*[Bs(:,:,e)]*[R(:,:,e)];

    % c1.2) Membrane transverse component of stiffness matrix:
    Bmt_i = zeros(1,5,4);
    for i = 1:4
        Bmt_i(:,:,i) = [N1_xmat(2,i), N1_xmat(1,i), 0, 0, 0];
    end %loop over nodes

    Cmt = h(Tm(e))*E/(2*(1+nu));
    Bmt(:,:,e) = [Bmt_i(:,:,1),Bmt_i(:,:,2),Bmt_i(:,:,3),Bmt_i(:,:,4)];
    Km(:,:,e) = S_1*[R(:,:,e)]'*[Bmt(:,:,e)]'*Cmt*[Bmt(:,:,e)]*[R(:,:,e)];
    
    % c2) 4 Gauss points quadrature matrices:
    Kb(:,:,e) = zeros(24,24);
    M_e(:,:,e) = zeros(24,24);
    ksi_4 = [-1, 1, 1, -1]/sqrt(3);
    eta_4 = [-1, -1, 1, 1]/sqrt(3);
    w_4 = [1, 1, 1, 1];

    for k = 1:4
        J_4 = zeros(2,2);
        for i = 1:4
            N4(i) = (1 + a(i)*ksi_4(k))*(1 + b(i)*eta_4(k))/4;
            N4_ksi(1,i) = a(i)*(1 + b(i)*eta_4(k))/4;
            N4_eta(1,i) = b(i)*(1 + a(i)*ksi_4(k))/4;
            J_4 = J_4 + [N4_ksi(i); 
                         N4_eta(i)]*xn(Tn(e,i),:)*[i_vec, j_vec];
        end

        N4x_mat = J_4^(-1)*[N4_ksi;N4_eta];
        S_4(e,k) = w_4(k)*det(J_4); %area associated to Gauss point

        % c2.1) Membrane normal component of stiffness matrix:
        Bmn_i = zeros(2,5,4);
        for i = 1:4
            Bmn_i(:,:,i) = [N4x_mat(1,i),    0,   0, 0, 0;
                               0,   N4x_mat(2,i), 0, 0, 0];
        end %loop over nodes
        Cmn = [1, nu; 
               nu, 1] *h(Tm(e))*E/(1-nu^2);
        Bmn(:,:,e,k) = [Bmn_i(:,:,1),Bmn_i(:,:,2),Bmn_i(:,:,3),Bmn_i(:,:,4)]; 
        Km(:,:,e) = Km(:,:,e) + S_4(e,k)*[R(:,:,e)]'*[Bmn(:,:,e,k)]'*Cmn*[Bmn(:,:,e,k)]*[R(:,:,e)];
        
        % c2.2) Bending component of stiffness matrix:
        Bb_i = zeros(3,5,4);
        for i = 1:4
            Bb_i(:,:,i) = [0, 0, 0,       0,       N4x_mat(1,i);
                           0, 0, 0, N4x_mat(2,i),      0;
                           0, 0, 0, -N4x_mat(1,i), N4x_mat(2,i)];
        end %loop over nodes

        Cb = [1, nu, 0; 
             nu,  1, 0; 
              0,  0, (1-nu)/2]   *h(Tm(e))^3*E/(12*(1-nu^2));
        Bb(:,:,e,k) = [Bb_i(:,:,1),Bb_i(:,:,2),Bb_i(:,:,3),Bb_i(:,:,4)];
        Kb(:,:,e) = Kb(:,:,e) + S_4(e,k)*[R(:,:,e)]'*[Bb(:,:,e,k)]'*Cb*[Bb(:,:,e,k)]*[R(:,:,e)];
        
        % c2.3) Mass matrix:
        for i = 1:4
            N_i(:,:,i) = N4(i)*eye(5);
        end

        rho_mat = [1, 0, 0,         0,         0;
                   0, 1, 0,         0,         0;
                   0, 0, 1,         0,         0;
                   0, 0, 0, (h(Tm(e))^2)/12,   0;
                   0, 0, 0,         0, (h(Tm(e))^2)/12]*rho*h(Tm(e));

        N(:,:,e,k) = [N_i(:,:,1),N_i(:,:,2),N_i(:,:,3),N_i(:,:,4)];
        M_e(:,:,e) = M_e(:,:,e) + S_4(e,k)*[R(:,:,e)]'*[N(:,:,e,k)]'*rho_mat*[N(:,:,e,k)]*[R(:,:,e)];
    end %loop over nodes

    % 2.2 d) Assembly to global matrices
    for j = 1:6
        Idof(j,1)    = 6*(Tn(e,1)-1) + j;
        Idof(6+j,1)  = 6*(Tn(e,2)-1) + j;
        Idof(12+j,1) = 6*(Tn(e,3)-1) + j;
        Idof(18+j,1) = 6*(Tn(e,4)-1) + j;
    end %loop over DOFs

    K(Idof,Idof) = K(Idof,Idof) + Km(:,:,e) + Kb(:,:,e) + Ks(:,:,e);
    M(Idof,Idof) = M(Idof,Idof) + M_e(:,:,e);
end

% 3) Compute artificial rotation stiffness matrix:
% 3.1 Find nodal normal to set criteria for finding coplanar nodes:
n = zeros(3,N_nod);

for e = 1:N_elem
    % a) Compute normal and surface
    s = cross((xn(Tn(e,3),:)'-xn(Tn(e,1),:)'),(xn(Tn(e,4),:)'-xn(Tn(e,2),:)'))/2;
    Se(e) = sqrt((s(1))^2+(s(2))^2+(s(3))^2);
    k_vec(:,e) = s/Se(e);

    % b) Assemble to get nodal normal
    for i = 1:4
      n(:,Tn(e,i)) = n(:,Tn(e,i)) + k_vec(:,e);
    end %loop over element nodes       
end % loop over elements

% 3.2) Compute artificial rotation matrix
Kr = sparse(N_dof,N_dof);
for e = 1:N_elem
    for i = 1:4
        % a) Determine whether it is or not a coplanar node
        Tnb=0; %revisar
        ind_beam = ismember(Tn(e,i),Tnb(:)); %only correct if node is not already a beam node
        alpha = acosd(dot(n(:,Tn(e,i)),k_vec(:,e))/norm(n(:,Tn(e,i))));
       if alpha<5  && ind_beam == false %we can consider node coplanar

            % b) Evaluate artificial rotation stiffness component
            Idof = 6*(Tn(e,i)-1) + [4, 5, 6]';
            Kr(Idof,Idof) = Kr(Idof,Idof) + E*h(Tm(e))*Se(e)*k_vec(:,e)*k_vec(:,e)';
        end %if
    end %loop over element nodes   
end  %loop over elements

% 3.3) Update stiffness matrix
K = K + Kr;

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
b = zeros(24, N_elem);
p = zeros(24, N_elem);
fe = zeros(24, N_elem);
for e = 1:N_elem
    % a) Compute element force vector:
    b(:,e) = [B(Tn(e,1),:),B(Tn(e,2),:),B(Tn(e,3),:),B(Tn(e,4),:)]';
    p(:,e) = [P(Tn(e,1),:),P(Tn(e,2),:),P(Tn(e,3),:),P(Tn(e,4),:)]';
    fe(:,e) = [M_e(:,:,e)]*b(:,e);

for k = 1:4
    fe(:,e) = fe(:,e) + S_4(e,k)*[R(:,:,e)]'*[N(:,:,e,k)]'*[N(:,:,e,k)]*[R(:,:,e)]*p(:,e);
end %loop over Gauss points

% b) Assembly to global force vector:
for j = 1:6
        Idof(j,1)    = 6*(Tn(e,1)-1) + j;
        Idof(6+j,1)  = 6*(Tn(e,2)-1) + j;
        Idof(12+j,1) = 6*(Tn(e,3)-1) + j;
        Idof(18+j,1) = 6*(Tn(e,4)-1) + j;
end %loop over DOFs

f_vec(Idof,1) = f_vec(Idof,1) + fe(:,e);
end
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
% Symmetrize and sparsify stiffness and mass matrices
    K = sparse(0.5 * (K + K.'));
    M = sparse(0.5 * (M + M.'));

    % Perform eigenvalue analysis
    [V, D] = eigs(K(If, If), M(If, If), N_modes, 'sm');

    % Preallocate outputs
    Phi = zeros(N_dof, N_modes);
    lambda = zeros(1,N_modes);                
    freq = zeros(N_modes,1);  

    % Normalize eigenvectors and populate global mode shapes
    Minv = M(If, If);                  
    for k = 1:size(V, 2)
        v = V(:, k);
        Phi(If, k) = v / sqrt(v' * Minv * v);
        lambda(k) = D(k, k);
        freq(k) = sqrt(lambda(k)) / (2 * pi);
    end

% 6) Solve system of equations (static case)
% 6.1 Solve system
u_vec(If,1) = [K(If,If)]\(f_vec(If,1)-[K(If,Ip)]*u_vec(Ip,1));
f_r = K*u_vec-f_vec;

% 7) Postprocess: Computing local strain and stress in shell elements
% 7.1) Get stress and strain at each Gauss point
for e = 1:N_elem
    % a) Get each strain component:
    for j = 1:6
        Idof(j,1)    = 6*(Tn(e,1)-1) + j;
        Idof(6+j,1)  = 6*(Tn(e,2)-1) + j;
        Idof(12+j,1) = 6*(Tn(e,3)-1) + j;
        Idof(18+j,1) = 6*(Tn(e,4)-1) + j;
    end %loop over DOFs

    for k = 1:4
        strain_b(:,e,k)   = [Bb(:,:,e,k)]*[R(:,:,e)]*u_vec(Idof,1);
        strain_m(1:2,e,k) = [Bmn(:,:,e,k)]*[R(:,:,e)]*u_vec(Idof,1);
        strain_m(3,e,k)   = [Bmt(:,:,e)]*[R(:,:,e)]*u_vec(Idof,1);
        strain_s(:,e,k)   = [Bs(:,:,e)]*[R(:,:,e)]*u_vec(Idof,1);
    end

    % b) Get stress
    Cp = [1, nu, 0;
          nu, 1, 0;
          0,  0, (1-nu)/2] *E/(1-nu^2);

    Cs = [1, 0; 
          0, 1] *E/(2*(1+nu));

    for k = 1:4
        sigma_m(:,e,k) = Cp*strain_m(:,e,k); %(constant membrane stress over the thickness)
        sigma_s(:,e,k) = Cs*strain_s(:,e,k); %(constant shear stress over the thickness assumed)
        sigma_b(:,e,k) = Cp*h(Tm(e))*strain_b(:,e,k)/2; %(bending stress on the top surface)
        sigma_plus = [sigma_m(:,e,k) + sigma_b(:,e,k); sigma_s(:,e,k)]'; %(stress on the top surface)
        sigma_VMplus = (sigma_plus(1)^2 + sigma_plus(2)^2 -sigma_plus(1)*sigma_plus(2) + 3*(sigma_plus(3)^2 + sigma_plus(4)^2 + sigma_plus(5)^2))^(1/2);
        sigma_minus = [sigma_m(:,e,k) - sigma_b(:,e,k); sigma_s(:,e,k)]'; %(stress on the bottom surface)
        sigma_VMminus = (sigma_minus(1)^2 + sigma_minus(2)^2 -sigma_minus(1)*sigma_minus(2) + 3*(sigma_minus(3)^2 + sigma_minus(4)^2 + sigma_minus(5)^2))^(1/2);
        sigma_VM(e,k) = max(sigma_VMplus,sigma_VMminus);
    end %loop over Gauss points
end %loop over elements

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

% figure
% hold on
% for mode = 1:N_modes
% plot(xS,u_z_avg(:,mode),'LineWidth',2)
% end
% legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6',Location='best')
% xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Vertical deflection, $u_z$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
% title('Vertical deflection distribution along the wingspan', 'FontWeight', 'bold', 'FontSize', 14);
% grid on
% grid minor
% box on
% axis padded
% %xlim([0,b]);
% hold off
% 
% 
% figure
% hold on
% for mode = 1:N_modes
% plot(xS,u_y_avg(:,mode),'LineWidth',2)
% end
% legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6',Location='best')
% xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Horizontal deflection, $u_y$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
% title('Horizontal deflection distribution along the wingspan', 'FontWeight', 'bold', 'FontSize', 14);
% grid on
% grid minor
% box on
% axis padded
% %xlim([0,b]);
% hold off




% figure
% hold on
% for mode = 1:N_modes
% plot(xS,theta_x_avg(:,mode),'LineWidth',2)
% end
% legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5','Mode 6',Location='best')
% xlabel('Span position, $x$ (m)', 'Interpreter', 'latex', 'FontSize', 14);
% ylabel('Twist angle $\theta_x$ (rad)', 'Interpreter', 'latex', 'FontSize', 14)
% title('Twist angle distribution along the wingspan', 'FontWeight', 'bold', 'FontSize', 14);
% grid on
% grid minor
% box on
% axis padded
% %xlim([0,b]);

% --- 3 subgráficas (u_y , u_z , theta_x) para cada modo ---
for iMode = 1:N_modes
    
    figure('Name',['Average fields – Mode ',num2str(iMode)],'Position',[200 200 560 420]);

    % ---------- u_y -----------------------------------------------------
    subplot(3,1,1)
    plot(xS, u_y_avg(:, iMode), 'LineWidth', 1);
    grid on;  box on
    title(['$u_y$ - Mode ', num2str(iMode)], 'Interpreter', 'latex', 'FontSize',12,'FontWeight','bold');  % ← sin LaTeX
    ylabel('u_y  [m]', 'FontSize',11);
    axis padded
    xlim([0 5])

    % ---------- u_z -----------------------------------------------------
    subplot(3,1,2)
    plot(xS, u_z_avg(:, iMode), 'LineWidth', 1);
    grid on;  box on
    title(['$u_z$ - Mode ', num2str(iMode)], 'Interpreter', 'latex', 'FontSize',12,'FontWeight','bold');
    ylabel('u_z  [m]', 'FontSize',11);
    axis padded
    xlim([0 5])
    
    % ---------- theta_x -------------------------------------------------
    subplot(3,1,3)
    plot(xS, theta_x_avg(:, iMode), 'LineWidth', 1);
    grid on;  box on
    title(['$\theta_x$ - Mode ', num2str(iMode)], 'Interpreter', 'latex', 'FontSize',12,'FontWeight','bold');
    ylabel('\theta_x  [rad]', 'FontSize',11);
    xlabel('Spanwise position  x  [m]', 'FontSize',11);
    axis padded
    xlim([0 5])
    
end


%% Vertical deflection  (u_z)  +  twist angle  (theta_x)  in a single plot
xSpan = xn(indSpar1,1);     % posición along spanwise
figure;hold on
yyaxis left
h0 = plot(xSpan, zeros(size(xSpan)), 'k--', 'LineWidth',1);    % undeformed
h1 = plot(xSpan, u_z,  'b', 'LineWidth',2);                    % u_z deformed
ylabel('Vertical deflection $u_z$ [m]', 'Interpreter','latex','FontSize',13);
ylim padded
ax = gca;                      
ax.YColor = 'b';
yyaxis right
h2 = plot(xSpan, theta_x, 'r', 'LineWidth',2);                 % theta_x
ylabel('Twist angle $\theta_x$ [rad]','Interpreter','latex','FontSize',13);
ylim padded
ax.YColor = 'r';
xlabel('Spanwise position $x$ [m]', 'Interpreter','latex','FontSize',13);
title('$u_z$ and $\theta_{x}$ along the span', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
legend([h0 h1 h2], {'Original', '$u_z$ (deformed)', '$\theta_x$ (twist)'}, 'Interpreter','latex','Location','best');
grid on;  grid minor;  box on;  axis padded;
xlim([0 5]);


scale=3e5;
plotDeformed('shell',xn,Tn,u_vec,scale);

imodes = [1,2,3,4,5,6];
plotModes('shell',Phi,freq,imodes)