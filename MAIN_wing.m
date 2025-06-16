%% PART III - WING MODELLING

clear
close all

%% DATA
% Define the problem's data (e.g. dimensions, material properties, etc.)
E   = 68.8e9;       % Young's modulus [Pa]
nu  = 0.33;         % Poisson ratio [-]
G   = E/(2*(1+nu)); % Shear modulus [Pa]
rho = 2700;         % Density [kg/m^3]
j_vec = [0;1;0];    % Orientation of y' axis
g   = 9.81;

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
Iy = 0.25*pi*(D/2)^4; % Area inertia y [m^4]
Iz = Iy;             % Area inertia z [m^4] (circular)
J   = Iy+Iz;         % Polar inertia  [m^4]
ky_st = 5/6;           % Shear correction factor k_y' (stringer)
kz_st = 5/6;           % Shear correction factor k_z' (stringer) 
kt_st = 1;             % Torsion correction factor k_t (stringer)


%% PREPROCESS

% Load mesh data
load('wing.mat','xn','Tn_st','Tm_st','Tn_wb','Tm_wb','Tn_rb','Tm_rb','Tn_sk','Tm_sk','indRoot','indPointA','indPointB','indSpar1','indSpar2','n_u','n_l');
N_nod = size(xn,1);
N_dof = 6*N_nod; 

onoff = struct('wb',1, 'sk',1, 'rb',1, 'st',1); %Turn off some elements
case_load = 'tunnel'; % Options 'uForce', 'uTorque', 'tunnel'.

% 1.3 Boundary conditions Up (fix nodes and DOFs matrix)
n = size(indRoot,1);
Up = [zeros(n*6,1), repelem(indRoot,6,1), repmat((1:6)',n,1)]; %always fix nodes at root
Qe = zeros(0,3);

switch case_load
    case 'uForce' %F_z=-1
        F_A = -(y2-yc)/(y2-y1); 
        F_B = -(yc-y1)/(y2-y1);
        Fe =[F_A, indPointA, 3; F_B, indPointB, 3];
        [Pe, Be] = deal(zeros(0,3));

    case 'uTorque'
        F_A = -1/(y2-y1);
        F_B = +1/(y2-y1);
        Fe =[F_A, indPointA, 3; F_B, indPointB, 3];
        [Pe, Be] = deal(zeros(0,3));

    case 'tunnel'
        p_inf    = 2.5e5;     % [Pa]
        alpha   = deg2rad(9); % [rad]

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
       Pe = [
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
% Initialization
N_elem_st       = numel(Tm_st);         
K_st       = sparse(N_dof,N_dof); %using sparse
M_st       = sparse(N_dof,N_dof);

%Matrix definitions
R_mat_1   = zeros(6,6);
R_st       = zeros(12,12,N_elem_st);
Ba_st      = zeros(1,12,N_elem_st);
Bb_st      = zeros(2,12,N_elem_st);
Bs_st      = zeros(2,12,N_elem_st);
Bt_st      = zeros(1,12,N_elem_st);

%2.2 Assembly process:
for e = 1:N_elem_st
  
 % a) Compute rotation matrix:
    le_st = norm(xn(Tn_st(e,2),:)-xn(Tn_st(e,1),:)); % (element size = beam length) [m]
    i_vec = (xn(Tn_st(e,2),:)'-xn(Tn_st(e,1),:)')/le_st;
    %j_vec defined already [0, 1, 0]
    k_vec = cross(i_vec,j_vec);
    R_mat_1(1:3,1:3) = [i_vec, j_vec, k_vec];
    R_mat_1(4:6,4:6) = [i_vec, j_vec, k_vec];
    R_mat_1 = R_mat_1';
    R_st(:,:,e) = [R_mat_1, zeros(6,6); zeros(6,6), R_mat_1];

  % b) Compute shape function derivatives:
    N_x_st(1) = -1/le_st;
    N_x_st(2) = 1/le_st;

  % c) Compute each element matrix:
    % c1) Axial component of stiffness matrix:
    Ba_st(1,:,e) = [N_x_st(1), 0, 0, 0, 0, 0, N_x_st(2), 0, 0, 0, 0, 0];
    Ca_st = E*A;
    Ka_st(:,:,e) = le_st*[R_st(:,:,e)]'*[Ba_st(1,:,e)]'*Ca_st*[Ba_st(1,:,e)]*[R_st(:,:,e)];

    % c2) Bending component of stiffness matrix:
    Bb_st(:,:,e) = [0, 0, 0, 0, N_x_st(1), 0, 0, 0, 0, 0, N_x_st(2), 0;
                 0, 0, 0, 0, 0, N_x_st(1), 0, 0, 0, 0, 0, N_x_st(2)];
    Cb_st = E* [Iy, 0; 0, Iz];
    Kb_st(:,:,e) = le_st*[R_st(:,:,e)]'*[Bb_st(:,:,e)]'*Cb_st*[Bb_st(:,:,e)]*[R_st(:,:,e)];

    % c3) Shear component of stiffness matrix
    N_i = 1/2; %(shape functions assuming only 1 Gauss point)
    Bs_st(:,:,e) = [0, N_x_st(1), 0, 0, 0, -N_i, 0, N_x_st(2), 0, 0, 0, -N_i;
                 0, 0, N_x_st(1), 0, N_i, 0, 0, 0, N_x_st(2), 0, N_i, 0 ];
    Cs_st = G*A*[ky_st, 0; 0, kz_st];
    Ks_st(:,:,e) = le_st*[R_st(:,:,e)]'*[Bs_st(:,:,e)]'*Cs_st*[Bs_st(:,:,e)]*[R_st(:,:,e)];

    % c4) Torsion component of stiffness matrix:
    Bt_st(1,:,e) = [0, 0, 0, N_x_st(1), 0, 0, 0, 0, 0, N_x_st(2), 0, 0];
    Ct_st = G*J*kt_st;
    Kt_st(:,:,e) = le_st*[R_st(:,:,e)]'*[Bt_st(1,:,e)]'*Ct_st*[Bt_st(1,:,e)]*[R_st(:,:,e)];

    % c5) Mass matrix:
    ksi = [-1/sqrt(3);1/sqrt(3)]; %Gauss point coordinates
    w = [1;1];
    rho_mat = rho*diag([A, A, A, J, Iy, Iz]);
    Me_st(:,:,e) = zeros(12,12);
    for k = 1:2
        N_i_st(1) = (1-ksi(k))/2;
        N_i_st(2) = (1+ksi(k))/2;
        N_st(:,:,e,k) = [N_i_st(1)*diag([ones(1,6)]), N_i_st(2)*diag([ones(1,6)])];
        Me_st(:,:,e) = Me_st(:,:,e) + w(k)*le_st*[R_st(:,:,e)]'*[N_st(:,:,e,k)]'*rho_mat*[N_st(:,:,e,k)]*[R_st(:,:,e)]/2;
    end

    % d) Assembly to global matrices
    Idof_st = zeros(12,1);
    for j = 1:6
        Idof_st(j,1) = 6*(Tn_st(e,1)-1) + j;
        Idof_st(6+j,1) = 6*(Tn_st(e,2)-1) + j;
    end

    K_st(Idof_st,Idof_st) = K_st(Idof_st,Idof_st) + Ka_st(:,:,e) + Kb_st(:,:,e) + Ks_st(:,:,e) + Kt_st(:,:,e);
    M_st(Idof_st,Idof_st) = M_st(Idof_st,Idof_st) + Me_st(:,:,e);

end %loop over elements


% Wingbox (with shell elements)

N_elem_wb = numel(Tm_wb);

K_wb = sparse(N_dof,N_dof);
M_wb = sparse(N_dof,N_dof);

% 2.2) Assembly process:
for e = 1:N_elem_wb

    % 2.2 a) Compute rotation matrix:
    S_wb = cross((xn(Tn_wb(e,3),:)'-xn(Tn_wb(e,1),:)'),(xn(Tn_wb(e,4),:)'-xn(Tn_wb(e,2),:)'))/2;
    k_vec_wb = S_wb/norm(S_wb); %mormal vector of the flat shell element
    d_wb = (xn(Tn_wb(e,2),:)' + xn(Tn_wb(e,3),:)' - xn(Tn_wb(e,4),:)' - xn(Tn_wb(e,1),:)')/2;
    i_vec_wb = d_wb/norm(d_wb);
    j_vec_wb = cross(k_vec_wb,i_vec_wb);

    R_mat_wb = [i_vec_wb j_vec_wb k_vec_wb zeros(3,2); zeros(3,3) i_vec_wb j_vec_wb]';
    R_wb(:,:,e) = blkdiag(R_mat_wb, R_mat_wb, R_mat_wb, R_mat_wb);

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
        J_1 = J_1 + [N_1ksi(i); N_1eta(i)]*xn(Tn_wb(e,i),:)*[i_vec_wb, j_vec_wb];
    end %loop over nodes

    N1_xmat = J_1^(-1)*[N_1ksi; N_1eta];
    S_1 = 4*det(J_1); %area associated to Gauss point

    % c1.1) Shear component of stiffness matrix:
    Bs_i_wb = zeros(2,5,4);
    for i = 1:4
        Bs_i_wb(:,:,i) = [0, 0, N1_xmat(1,i), 0, N_1(i);
                       0, 0, N1_xmat(2,i), -N_1(i), 0];
    end %loop over nodes

    Cs_wb = [1, 0; 
          0, 1]*5*h(Tm_wb(e))*E/(12*(1 + nu)); %E,nu=ct.

    Bs_wb(:,:,e) = [Bs_i_wb(:,:,1),Bs_i_wb(:,:,2),Bs_i_wb(:,:,3),Bs_i_wb(:,:,4)];

    Ks_wb(:,:,e) = S_1*[R_wb(:,:,e)]'*[Bs_wb(:,:,e)]'*Cs_wb*[Bs_wb(:,:,e)]*[R_wb(:,:,e)];

    % c1.2) Membrane transverse component of stiffness matrix:
    Bmt_i_wb = zeros(1,5,4);
    for i = 1:4
        Bmt_i_wb(:,:,i) = [N1_xmat(2,i), N1_xmat(1,i), 0, 0, 0];
    end %loop over nodes

    Cmt_wb = h(Tm_wb(e))*E/(2*(1+nu));
    Bmt_wb(:,:,e) = [Bmt_i_wb(:,:,1),Bmt_i_wb(:,:,2),Bmt_i_wb(:,:,3),Bmt_i_wb(:,:,4)];
    Km_wb(:,:,e) = S_1*[R_wb(:,:,e)]'*[Bmt_wb(:,:,e)]'*Cmt_wb*[Bmt_wb(:,:,e)]*[R_wb(:,:,e)];
    
    % c2) 4 Gauss points quadrature matrices:
    Kb_wb(:,:,e) = zeros(24,24);
    Me_wb(:,:,e) = zeros(24,24);
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
                         N4_eta(i)]*xn(Tn_wb(e,i),:)*[i_vec_wb, j_vec_wb];
        end

        N4x_mat = J_4^(-1)*[N4_ksi;N4_eta];
        S4_wb(e,k) = w_4(k)*det(J_4); %area associated to Gauss point

        % c2.1) Membrane normal component of stiffness matrix:
        Bmn_i = zeros(2,5,4);
        for i = 1:4
            Bmn_i(:,:,i) = [N4x_mat(1,i),    0,   0, 0, 0;
                               0,   N4x_mat(2,i), 0, 0, 0];
        end %loop over nodes
        Cmn = [1, nu; 
               nu, 1] *h(Tm_wb(e))*E/(1-nu^2);
        Bmn_wb(:,:,e,k) = [Bmn_i(:,:,1),Bmn_i(:,:,2),Bmn_i(:,:,3),Bmn_i(:,:,4)]; 
        Km_wb(:,:,e) = Km_wb(:,:,e) + S4_wb(e,k)*[R_wb(:,:,e)]'*[Bmn_wb(:,:,e,k)]'*Cmn*[Bmn_wb(:,:,e,k)]*[R_wb(:,:,e)];
        
        % c2.2) Bending component of stiffness matrix:
        Bb_i = zeros(3,5,4);
        for i = 1:4
            Bb_i(:,:,i) = [0, 0, 0,       0,       N4x_mat(1,i);
                           0, 0, 0, N4x_mat(2,i),      0;
                           0, 0, 0, -N4x_mat(1,i), N4x_mat(2,i)];
        end %loop over nodes

        Cb_wb = [1, nu, 0; 
             nu,  1, 0; 
              0,  0, (1-nu)/2]   *h(Tm_wb(e))^3*E/(12*(1-nu^2));
        Bb_wb(:,:,e,k) = [Bb_i(:,:,1),Bb_i(:,:,2),Bb_i(:,:,3),Bb_i(:,:,4)];
        Kb_wb(:,:,e) = Kb_wb(:,:,e) + S4_wb(e,k)*[R_wb(:,:,e)]'*[Bb_wb(:,:,e,k)]'*Cb_wb*[Bb_wb(:,:,e,k)]*[R_wb(:,:,e)];
        
        % c2.3) Mass matrix:
        for i = 1:4
            N_i_wb(:,:,i) = N4(i)*eye(5);
        end

        rho_mat = [1, 0, 0,         0,         0;
                   0, 1, 0,         0,         0;
                   0, 0, 1,         0,         0;
                   0, 0, 0, (h(Tm_wb(e))^2)/12,   0;
                   0, 0, 0,         0, (h(Tm_wb(e))^2)/12]*rho*h(Tm_wb(e));

        N_wb(:,:,e,k) = [N_i_wb(:,:,1),N_i_wb(:,:,2),N_i_wb(:,:,3),N_i_wb(:,:,4)];
        Me_wb(:,:,e) = Me_wb(:,:,e) + S4_wb(e,k)*[R_wb(:,:,e)]'*[N_wb(:,:,e,k)]'*rho_mat*[N_wb(:,:,e,k)]*[R_wb(:,:,e)];
    end %loop over nodes

    % 2.2 d) Assembly to global matrices
    for j = 1:6
        Idof_wb(j,1)    = 6*(Tn_wb(e,1)-1) + j;
        Idof_wb(6+j,1)  = 6*(Tn_wb(e,2)-1) + j;
        Idof_wb(12+j,1) = 6*(Tn_wb(e,3)-1) + j;
        Idof_wb(18+j,1) = 6*(Tn_wb(e,4)-1) + j;
    end %loop over DOFs

    K_wb(Idof_wb,Idof_wb) = K_wb(Idof_wb,Idof_wb) + Km_wb(:,:,e) + Kb_wb(:,:,e) + Ks_wb(:,:,e);
    M_wb(Idof_wb,Idof_wb) = M_wb(Idof_wb,Idof_wb) + Me_wb(:,:,e);
end

% 3) Compute artificial rotation stiffness matrix:
% 3.1 Find nodal normal to set criteria for finding coplanar nodes:
n = zeros(3,N_nod);

for e = 1:N_elem_wb
    % a) Compute normal and surface
    s = cross((xn(Tn_wb(e,3),:)'-xn(Tn_wb(e,1),:)'),(xn(Tn_wb(e,4),:)'-xn(Tn_wb(e,2),:)'))/2;
    Se(e) = sqrt((s(1))^2+(s(2))^2+(s(3))^2);
    k_vec_wb(:,e) = s/Se(e);

    % b) Assemble to get nodal normal
    for i = 1:4
      n(:,Tn_wb(e,i)) = n(:,Tn_wb(e,i)) + k_vec_wb(:,e);
    end %loop over element nodes       
end % loop over elements

% 3.2) Compute artificial rotation matrix
Kr_wb = sparse(N_dof,N_dof);
for e = 1:N_elem_wb
    for i = 1:4
        % a) Determine whether it is or not a coplanar node
        ind_beam_wb = ismember(Tn_wb(e,i),Tn_st(:)); %only correct if node is not already a beam node
        alpha_wb = acosd(dot(n(:,Tn_wb(e,i)),k_vec_wb(:,e))/norm(n(:,Tn_wb(e,i))));
       if alpha_wb<5  && ind_beam_wb == false %we can consider node coplanar

            % b) Evaluate artificial rotation stiffness component
            Idof_wb = 6*(Tn_wb(e,i)-1) + [4, 5, 6]';
            Kr_wb(Idof_wb,Idof_wb) = Kr_wb(Idof_wb,Idof_wb) + E*h(Tm_wb(e))*Se(e)*k_vec_wb(:,e)*k_vec_wb(:,e)';
        end %if
    end %loop over element nodes   
end  %loop over elements

% 3.3) Update stiffness matrix
K_wb = K_wb + Kr_wb;


% Skin (with shell elements)

% 2.1) Initialization:
N_elem_sk = numel(Tm_sk);

K_sk = sparse(N_dof,N_dof);
M_sk = sparse(N_dof,N_dof);

% 2.2) Assembly process:
for e = 1:N_elem_sk

    % 2.2 a) Compute rotation matrix:
    S_sk = cross((xn(Tn_sk(e,3),:)'-xn(Tn_sk(e,1),:)'),(xn(Tn_sk(e,4),:)'-xn(Tn_sk(e,2),:)'))/2;
    k_vec_sk = S_sk/norm(S_sk); %mormal vector of the flat shell element
    d_sk = (xn(Tn_sk(e,2),:)' + xn(Tn_sk(e,3),:)' - xn(Tn_sk(e,4),:)' - xn(Tn_sk(e,1),:)')/2;
    i_vec_sk = d_sk/norm(d_sk);
    j_vec_sk = cross(k_vec_sk,i_vec_sk);

    R_mat_sk = [i_vec_sk j_vec_sk k_vec_sk zeros(3,2); zeros(3,3) i_vec_sk j_vec_sk]';
    R_sk(:,:,e) = blkdiag(R_mat_sk, R_mat_sk, R_mat_sk, R_mat_sk);

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
        J_1 = J_1 + [N_1ksi(i); N_1eta(i)]*xn(Tn_sk(e,i),:)*[i_vec_sk, j_vec_sk];
    end %loop over nodes

    N1_xmat = J_1^(-1)*[N_1ksi; N_1eta];
    S_1 = 4*det(J_1); %area associated to Gauss point

    % c1.1) Shear component of stiffness matrix:
    Bs_i_sk = zeros(2,5,4);
    for i = 1:4
        Bs_i_sk(:,:,i) = [0, 0, N1_xmat(1,i), 0, N_1(i);
                       0, 0, N1_xmat(2,i), -N_1(i), 0];
    end %loop over nodes

    Cs_sk = [1, 0; 
          0, 1]*5*h(Tm_sk(e))*E/(12*(1 + nu)); %E,nu=ct.

    Bs_sk(:,:,e) = [Bs_i_sk(:,:,1),Bs_i_sk(:,:,2),Bs_i_sk(:,:,3),Bs_i_sk(:,:,4)];

    Ks_sk(:,:,e) = S_1*[R_sk(:,:,e)]'*[Bs_sk(:,:,e)]'*Cs_sk*[Bs_sk(:,:,e)]*[R_sk(:,:,e)];

    % c1.2) Membrane transverse component of stiffness matrix:
    Bmt_i_sk = zeros(1,5,4);
    for i = 1:4
        Bmt_i_sk(:,:,i) = [N1_xmat(2,i), N1_xmat(1,i), 0, 0, 0];
    end %loop over nodes

    Cmt_sk = h(Tm_sk(e))*E/(2*(1+nu));
    Bmt_sk(:,:,e) = [Bmt_i_sk(:,:,1),Bmt_i_sk(:,:,2),Bmt_i_sk(:,:,3),Bmt_i_sk(:,:,4)];
    Km_sk(:,:,e) = S_1*[R_sk(:,:,e)]'*[Bmt_sk(:,:,e)]'*Cmt_sk*[Bmt_sk(:,:,e)]*[R_sk(:,:,e)];
    
    % c2) 4 Gauss points quadrature matrices:
    Kb_sk(:,:,e) = zeros(24,24);
    Me_sk(:,:,e) = zeros(24,24);
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
                         N4_eta(i)]*xn(Tn_sk(e,i),:)*[i_vec_sk, j_vec_sk];
        end

        N4x_mat = J_4^(-1)*[N4_ksi;N4_eta];
        S4_sk(e,k) = w_4(k)*det(J_4); %area associated to Gauss point

        % c2.1) Membrane normal component of stiffness matrix:
        Bmn_i = zeros(2,5,4);
        for i = 1:4
            Bmn_i(:,:,i) = [N4x_mat(1,i),    0,   0, 0, 0;
                               0,   N4x_mat(2,i), 0, 0, 0];
        end %loop over nodes
        Cmn = [1, nu; 
               nu, 1] *h(Tm_sk(e))*E/(1-nu^2);
        Bmn_sk(:,:,e,k) = [Bmn_i(:,:,1),Bmn_i(:,:,2),Bmn_i(:,:,3),Bmn_i(:,:,4)]; 
        Km_sk(:,:,e) = Km_sk(:,:,e) + S4_sk(e,k)*[R_sk(:,:,e)]'*[Bmn_sk(:,:,e,k)]'*Cmn*[Bmn_sk(:,:,e,k)]*[R_sk(:,:,e)];
        
        % c2.2) Bending component of stiffness matrix:
        Bb_i = zeros(3,5,4);
        for i = 1:4
            Bb_i(:,:,i) = [0, 0, 0,       0,       N4x_mat(1,i);
                           0, 0, 0, N4x_mat(2,i),      0;
                           0, 0, 0, -N4x_mat(1,i), N4x_mat(2,i)];
        end %loop over nodes

        Cb_sk = [1, nu, 0; 
             nu,  1, 0; 
              0,  0, (1-nu)/2]   *h(Tm_sk(e))^3*E/(12*(1-nu^2));
        Bb_sk(:,:,e,k) = [Bb_i(:,:,1),Bb_i(:,:,2),Bb_i(:,:,3),Bb_i(:,:,4)];
        Kb_sk(:,:,e) = Kb_sk(:,:,e) + S4_sk(e,k)*[R_sk(:,:,e)]'*[Bb_sk(:,:,e,k)]'*Cb_sk*[Bb_sk(:,:,e,k)]*[R_sk(:,:,e)];
        
        % c2.3) Mass matrix:
        for i = 1:4
            N_i_sk(:,:,i) = N4(i)*eye(5);
        end

        rho_mat = [1, 0, 0,         0,         0;
                   0, 1, 0,         0,         0;
                   0, 0, 1,         0,         0;
                   0, 0, 0, (h(Tm_sk(e))^2)/12,   0;
                   0, 0, 0,         0, (h(Tm_sk(e))^2)/12]*rho*h(Tm_sk(e));

        N_sk(:,:,e,k) = [N_i_sk(:,:,1),N_i_sk(:,:,2),N_i_sk(:,:,3),N_i_sk(:,:,4)];
        Me_sk(:,:,e) = Me_sk(:,:,e) + S4_sk(e,k)*[R_sk(:,:,e)]'*[N_sk(:,:,e,k)]'*rho_mat*[N_sk(:,:,e,k)]*[R_sk(:,:,e)];
    end %loop over nodes

    % 2.2 d) Assembly to global matrices
    for j = 1:6
        Idof_sk(j,1)    = 6*(Tn_sk(e,1)-1) + j;
        Idof_sk(6+j,1)  = 6*(Tn_sk(e,2)-1) + j;
        Idof_sk(12+j,1) = 6*(Tn_sk(e,3)-1) + j;
        Idof_sk(18+j,1) = 6*(Tn_sk(e,4)-1) + j;
    end %loop over DOFs

    K_sk(Idof_sk,Idof_sk) = K_sk(Idof_sk,Idof_sk) + Km_sk(:,:,e) + Kb_sk(:,:,e) + Ks_sk(:,:,e);
    M_sk(Idof_sk,Idof_sk) = M_sk(Idof_sk,Idof_sk) + Me_sk(:,:,e);
end

% 3) Compute artificial rotation stiffness matrix:
% 3.1 Find nodal normal to set criteria for finding coplanar nodes:
n = zeros(3,N_nod);

for e = 1:N_elem_sk
    % a) Compute normal and surface
    s = cross((xn(Tn_sk(e,3),:)'-xn(Tn_sk(e,1),:)'),(xn(Tn_sk(e,4),:)'-xn(Tn_sk(e,2),:)'))/2;
    Se(e) = sqrt((s(1))^2+(s(2))^2+(s(3))^2);
    k_vec_sk(:,e) = s/Se(e);

    % b) Assemble to get nodal normal
    for i = 1:4
      n(:,Tn_sk(e,i)) = n(:,Tn_sk(e,i)) + k_vec_sk(:,e);
    end %loop over element nodes       
end % loop over elements

% 3.2) Compute artificial rotation matrix
Kr_sk = sparse(N_dof,N_dof);
for e = 1:N_elem_sk
    for i = 1:4
        % a) Determine whether it is or not a coplanar node
        ind_beam_sk = ismember(Tn_sk(e,i),Tn_st(:)); %only correct if node is not already a beam node
        alpha_sk = acosd(dot(n(:,Tn_sk(e,i)),k_vec_sk(:,e))/norm(n(:,Tn_sk(e,i))));
       if alpha_sk<5  && ind_beam_sk == false %we can consider node coplanar

            % b) Evaluate artificial rotation stiffness component
            Idof_sk = 6*(Tn_sk(e,i)-1) + [4, 5, 6]';
            Kr_sk(Idof_sk,Idof_sk) = Kr_sk(Idof_sk,Idof_sk) + E*h(Tm_sk(e))*Se(e)*k_vec_sk(:,e)*k_vec_sk(:,e)';
        end %if
    end %loop over element nodes   
end  %loop over elements

% 3.3) Update stiffness matrix
K_sk = K_sk + Kr_sk;

% Ribs (with shell elements)

% 2.1) Initialization:
N_elem_rb = numel(Tm_rb);

K_rb = sparse(N_dof,N_dof);
M_rb = sparse(N_dof,N_dof);

% 2.2) Assembly process:
for e = 1:N_elem_rb

    % 2.2 a) Compute rotation matrix:
    S_rb = cross((xn(Tn_rb(e,3),:)'-xn(Tn_rb(e,1),:)'),(xn(Tn_rb(e,4),:)'-xn(Tn_rb(e,2),:)'))/2;
    k_vec_rb = S_rb/norm(S_rb); %mormal vector of the flat shell element
    d_rb = (xn(Tn_rb(e,2),:)' + xn(Tn_rb(e,3),:)' - xn(Tn_rb(e,4),:)' - xn(Tn_rb(e,1),:)')/2;
    i_vec_rb = d_rb/norm(d_rb);
    j_vec_rb = cross(k_vec_rb,i_vec_rb);

    R_mat_rb = [i_vec_rb j_vec_rb k_vec_rb zeros(3,2); zeros(3,3) i_vec_rb j_vec_rb]';
    R_rb(:,:,e) = blkdiag(R_mat_rb, R_mat_rb, R_mat_rb, R_mat_rb);

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
        J_1 = J_1 + [N_1ksi(i); N_1eta(i)]*xn(Tn_rb(e,i),:)*[i_vec_rb, j_vec_rb];
    end %loop over nodes

    N1_xmat = J_1^(-1)*[N_1ksi; N_1eta];
    S_1_rb = 4*det(J_1); %area associated to Gauss point

    % c1.1) Shear component of stiffness matrix:
    Bs_i_rb = zeros(2,5,4);
    for i = 1:4
        Bs_i_rb(:,:,i) = [0, 0, N1_xmat(1,i), 0, N_1(i);
                       0, 0, N1_xmat(2,i), -N_1(i), 0];
    end %loop over nodes

    Cs_rb = [1, 0; 
          0, 1]*5*h(Tm_rb(e))*E/(12*(1 + nu)); %E,nu=ct.

    Bs_rb(:,:,e) = [Bs_i_rb(:,:,1),Bs_i_rb(:,:,2),Bs_i_rb(:,:,3),Bs_i_rb(:,:,4)];

    Ks_rb(:,:,e) = S_1_rb*[R_rb(:,:,e)]'*[Bs_rb(:,:,e)]'*Cs_rb*[Bs_rb(:,:,e)]*[R_rb(:,:,e)];

    % c1.2) Membrane transverse component of stiffness matrix:
    Bmt_i_rb = zeros(1,5,4);
    for i = 1:4
        Bmt_i_rb(:,:,i) = [N1_xmat(2,i), N1_xmat(1,i), 0, 0, 0];
    end %loop over nodes

    Cmt_rb = h(Tm_rb(e))*E/(2*(1+nu));
    Bmt_rb(:,:,e) = [Bmt_i_rb(:,:,1),Bmt_i_rb(:,:,2),Bmt_i_rb(:,:,3),Bmt_i_rb(:,:,4)];
    Km_rb(:,:,e) = S_1_rb*[R_rb(:,:,e)]'*[Bmt_rb(:,:,e)]'*Cmt_rb*[Bmt_rb(:,:,e)]*[R_rb(:,:,e)];
    
    % c2) 4 Gauss points quadrature matrices:
    Kb_rb(:,:,e) = zeros(24,24);
    Me_rb(:,:,e) = zeros(24,24);
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
                         N4_eta(i)]*xn(Tn_rb(e,i),:)*[i_vec_rb, j_vec_rb];
        end

        N4x_mat = J_4^(-1)*[N4_ksi;N4_eta];
        S4_rb(e,k) = w_4(k)*det(J_4); %area associated to Gauss point

        % c2.1) Membrane normal component of stiffness matrix:
        Bmn_i = zeros(2,5,4);
        for i = 1:4
            Bmn_i(:,:,i) = [N4x_mat(1,i),    0,   0, 0, 0;
                               0,   N4x_mat(2,i), 0, 0, 0];
        end %loop over nodes
        Cmn = [1, nu; 
               nu, 1] *h(Tm_rb(e))*E/(1-nu^2);
        Bmn_rb(:,:,e,k) = [Bmn_i(:,:,1),Bmn_i(:,:,2),Bmn_i(:,:,3),Bmn_i(:,:,4)]; 
        Km_rb(:,:,e) = Km_rb(:,:,e) + S4_rb(e,k)*[R_rb(:,:,e)]'*[Bmn_rb(:,:,e,k)]'*Cmn*[Bmn_rb(:,:,e,k)]*[R_rb(:,:,e)];
        
        % c2.2) Bending component of stiffness matrix:
        Bb_i = zeros(3,5,4);
        for i = 1:4
            Bb_i(:,:,i) = [0, 0, 0,       0,       N4x_mat(1,i);
                           0, 0, 0, N4x_mat(2,i),      0;
                           0, 0, 0, -N4x_mat(1,i), N4x_mat(2,i)];
        end %loop over nodes

        Cb_rb = [1, nu, 0; 
             nu,  1, 0; 
              0,  0, (1-nu)/2]   *h(Tm_rb(e))^3*E/(12*(1-nu^2));
        Bb_rb(:,:,e,k) = [Bb_i(:,:,1),Bb_i(:,:,2),Bb_i(:,:,3),Bb_i(:,:,4)];
        Kb_rb(:,:,e) = Kb_rb(:,:,e) + S4_rb(e,k)*[R_rb(:,:,e)]'*[Bb_rb(:,:,e,k)]'*Cb_rb*[Bb_rb(:,:,e,k)]*[R_rb(:,:,e)];
        
        % c2.3) Mass matrix:
        for i = 1:4
            N_i_rb(:,:,i) = N4(i)*eye(5);
        end

        rho_mat = [1, 0, 0,         0,         0;
                   0, 1, 0,         0,         0;
                   0, 0, 1,         0,         0;
                   0, 0, 0, (h(Tm_rb(e))^2)/12,   0;
                   0, 0, 0,         0, (h(Tm_rb(e))^2)/12]*rho*h(Tm_rb(e));

        N_rb(:,:,e,k) = [N_i_rb(:,:,1),N_i_rb(:,:,2),N_i_rb(:,:,3),N_i_rb(:,:,4)];
        Me_rb(:,:,e) = Me_rb(:,:,e) + S4_rb(e,k)*[R_rb(:,:,e)]'*[N_rb(:,:,e,k)]'*rho_mat*[N_rb(:,:,e,k)]*[R_rb(:,:,e)];
    end %loop over nodes

    % 2.2 d) Assembly to global matrices
    for j = 1:6
        Idof_rb(j,1)    = 6*(Tn_rb(e,1)-1) + j;
        Idof_rb(6+j,1)  = 6*(Tn_rb(e,2)-1) + j;
        Idof_rb(12+j,1) = 6*(Tn_rb(e,3)-1) + j;
        Idof_rb(18+j,1) = 6*(Tn_rb(e,4)-1) + j;
    end %loop over DOFs

    K_rb(Idof_rb,Idof_rb) = K_rb(Idof_rb,Idof_rb) + Km_rb(:,:,e) + Kb_rb(:,:,e) + Ks_rb(:,:,e);
    M_rb(Idof_rb,Idof_rb) = M_rb(Idof_rb,Idof_rb) + Me_rb(:,:,e);
end

% 3) Compute artificial rotation stiffness matrix:
% 3.1 Find nodal normal to set criteria for finding coplanar nodes:
n = zeros(3,N_nod);

for e = 1:N_elem_rb
    % a) Compute normal and surface
    s = cross((xn(Tn_rb(e,3),:)'-xn(Tn_rb(e,1),:)'),(xn(Tn_rb(e,4),:)'-xn(Tn_rb(e,2),:)'))/2;
    Se(e) = sqrt((s(1))^2+(s(2))^2+(s(3))^2);
    k_vec_rb(:,e) = s/Se(e);

    % b) Assemble to get nodal normal
    for i = 1:4
      n(:,Tn_rb(e,i)) = n(:,Tn_rb(e,i)) + k_vec_rb(:,e);
    end %loop over element nodes       
end % loop over elements

% 3.2) Compute artificial rotation matrix
Kr_rb = sparse(N_dof,N_dof);
for e = 1:N_elem_rb
    for i = 1:4
        % a) Determine whether it is or not a coplanar node
        ind_beam_rb = ismember(Tn_rb(e,i),Tn_st(:)); %only correct if node is not already a beam node
        alpha_rb = acosd(dot(n(:,Tn_rb(e,i)),k_vec_rb(:,e))/norm(n(:,Tn_rb(e,i))));
       if alpha_rb<5  && ind_beam_rb == false %we can consider node coplanar

            % b) Evaluate artificial rotation stiffness component
            Idof_rb = 6*(Tn_rb(e,i)-1) + [4, 5, 6]';
            Kr_rb(Idof_rb,Idof_rb) = Kr_rb(Idof_rb,Idof_rb) + E*h(Tm_rb(e))*Se(e)*k_vec_rb(:,e)*k_vec_rb(:,e)';
        end %if
    end %loop over element nodes   
end  %loop over elements

% 3.3) Update stiffness matrix
K_rb = K_rb + Kr_rb;

% K, M matrices

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

P = zeros(N_nod,6);
for r = 1:size(Pe,1)
    P(Pe(r,2),Pe(r,3)) = P(Pe(r,2),Pe(r,3)) + Pe(r,1);
end

% 3.3 Nodal body forces:
B = zeros(N_nod,6);
for s = 1:size(Be,1)
    B(Be(s,2),Be(s,3)) = B(Be(s,2),Be(s,3)) + Be(s,1);
end

%Assembly process:
N_st_elem = size(Tn_st,1); 
%  3.4 Assembly process:
b_st = zeros(12,N_st_elem);
q_st = zeros(12,N_st_elem);
fe_st = zeros(12,N_st_elem);
for e = 1:N_st_elem
    le_st = norm(xn(Tn_st(e,2),:)-xn(Tn_st(e,1),:)); %from previously
    % a) Compute element force vector:
    b_st(:,e) = [B(Tn_st(e,1),:), B(Tn_st(e,2),:)].';
    q_st(:,e) = [Q(Tn_st(e,1),:), Q(Tn_st(e,2),:)].';
    fe_st(:,e) = [Me_st(:,:,e)]*b_st(:,e);
    for k = 1:2
        fe_st(:,e) = fe_st(:,e) + w(k)*le_st*[R_st(:,:,e)].'*[N_st(:,:,e,k)].'*[N_st(:,:,e,k)]*[R_st(:,:,e)]*q_st(:,e)/2;
    end
    % b) Assembly to global force vector:
    for j = 1:6
            Idof_st(j,1) = 6*(Tn_st(e,1)-1) + j;
            Idof_st(6+j,1) = 6*(Tn_st(e,2)-1) + j;
    end
    f_vec(Idof_st,1) = f_vec(Idof_st,1) + fe_st(:,e);
end

N_elem_wb = size(Tn_wb,1);
% 4.4) Assembly process
for e = 1:N_elem_wb
    % a) Compute element force vector:
    b_wb(:,e) = [B(Tn_wb(e,1),:),B(Tn_wb(e,2),:),B(Tn_wb(e,3),:),B(Tn_wb(e,4),:)]';
    p_wb(:,e) = [P(Tn_wb(e,1),:),P(Tn_wb(e,2),:),P(Tn_wb(e,3),:),P(Tn_wb(e,4),:)]';
    fe_wb(:,e) = [Me_wb(:,:,e)]*b_wb(:,e);

for k = 1:4
    fe_wb(:,e) = fe_wb(:,e) + S4_wb(e,k)*[R_wb(:,:,e)]'*[N_wb(:,:,e,k)]'*[N_wb(:,:,e,k)]*[R_wb(:,:,e)]*p_wb(:,e);
end %loop over Gauss points

% b) Assembly to global force vector:
for j = 1:6
        Idof_wb(j,1)    = 6*(Tn_wb(e,1)-1) + j;
        Idof_wb(6+j,1)  = 6*(Tn_wb(e,2)-1) + j;
        Idof_wb(12+j,1) = 6*(Tn_wb(e,3)-1) + j;
        Idof_wb(18+j,1) = 6*(Tn_wb(e,4)-1) + j;
end %loop over DOFs

f_vec(Idof_wb,1) = f_vec(Idof_wb,1) + fe_wb(:,e);
end %loop over elements

N_elem_rb = size(Tn_rb,1);
% 4.4) Assembly process
for e = 1:N_elem_rb
    % a) Compute element force vector:
    b_rb(:,e) = [B(Tn_rb(e,1),:),B(Tn_rb(e,2),:),B(Tn_rb(e,3),:),B(Tn_rb(e,4),:)]';
    p_rb(:,e) = [P(Tn_rb(e,1),:),P(Tn_rb(e,2),:),P(Tn_rb(e,3),:),P(Tn_rb(e,4),:)]';
    fe_rb(:,e) = [Me_rb(:,:,e)]*b_rb(:,e);

for k = 1:4
    fe_rb(:,e) = fe_rb(:,e) + S4_rb(e,k)*[R_rb(:,:,e)]'*[N_rb(:,:,e,k)]'*[N_rb(:,:,e,k)]*[R_rb(:,:,e)]*p_rb(:,e);
end %loop over Gauss points

% b) Assembly to global force vector:
for j = 1:6
        Idof_rb(j,1)    = 6*(Tn_rb(e,1)-1) + j;
        Idof_rb(6+j,1)  = 6*(Tn_rb(e,2)-1) + j;
        Idof_rb(12+j,1) = 6*(Tn_rb(e,3)-1) + j;
        Idof_rb(18+j,1) = 6*(Tn_rb(e,4)-1) + j;
end %loop over DOFs

f_vec(Idof_rb,1) = f_vec(Idof_rb,1) + fe_rb(:,e);
end %loop over elements

N_elem_sk = size(Tn_sk,1);
% 4.4) Assembly process
for e = 1:N_elem_sk
    % a) Compute element force vector:
    b_sk(:,e) = [B(Tn_sk(e,1),:),B(Tn_sk(e,2),:),B(Tn_sk(e,3),:),B(Tn_sk(e,4),:)]';
    p_sk(:,e) = [P(Tn_sk(e,1),:),P(Tn_sk(e,2),:),P(Tn_sk(e,3),:),P(Tn_sk(e,4),:)]';
    fe_sk(:,e) = [Me_sk(:,:,e)]*b_sk(:,e);

for k = 1:4
    fe_sk(:,e) = fe_sk(:,e) + S4_sk(e,k)*[R_sk(:,:,e)]'*[N_sk(:,:,e,k)]'*[N_sk(:,:,e,k)]*[R_sk(:,:,e)]*p_sk(:,e);
end %loop over Gauss points

% b) Assembly to global force vector:
for j = 1:6
        Idof_sk(j,1)    = 6*(Tn_sk(e,1)-1) + j;
        Idof_sk(6+j,1)  = 6*(Tn_sk(e,2)-1) + j;
        Idof_sk(12+j,1) = 6*(Tn_sk(e,3)-1) + j;
        Idof_sk(18+j,1) = 6*(Tn_sk(e,4)-1) + j;
end %loop over DOFs

f_vec(Idof_sk,1) = f_vec(Idof_sk,1) + fe_sk(:,e);
end %loop over elements

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

% 7) Compute stresses
N_elem_wb = size(Tn_wb,1);
% 7) Postprocess: Computing local strain and stress in shell elements
% 7.1) Get stress and strain at each Gauss point
for e = 1:N_elem_wb
    % a) Get each strain component:
    for j = 1:6
        Idof_wb(j,1)    = 6*(Tn_wb(e,1)-1) + j;
        Idof_wb(6+j,1)  = 6*(Tn_wb(e,2)-1) + j;
        Idof_wb(12+j,1) = 6*(Tn_wb(e,3)-1) + j;
        Idof_wb(18+j,1) = 6*(Tn_wb(e,4)-1) + j;
    end %loop over DOFs

    for k = 1:4
        strain_b_wb(:,e,k)   = [Bb_wb(:,:,e,k)]*[R_wb(:,:,e)]*u_vec(Idof_wb,1);
        strain_m_wb(1:2,e,k) = [Bmn_wb(:,:,e,k)]*[R_wb(:,:,e)]*u_vec(Idof_wb,1);
        strain_m_wb(3,e,k)   = [Bmt_wb(:,:,e)]*[R_wb(:,:,e)]*u_vec(Idof_wb,1);
        strain_s_wb(:,e,k)   = [Bs_wb(:,:,e)]*[R_wb(:,:,e)]*u_vec(Idof_wb,1);
    end

    % b) Get stress
    Cp = [1, nu, 0;
          nu, 1, 0;
          0,  0, (1-nu)/2] *E/(1-nu^2);

    Cs = [1, 0; 
          0, 1] *E/(2*(1+nu));

    for k = 1:4
        sigma_m_wb(:,e,k) = Cp*strain_m_wb(:,e,k); %(constant membrane stress over the thickness)
        sigma_s_wb(:,e,k) = Cs*strain_s_wb(:,e,k); %(constant shear stress over the thickness assumed)
        sigma_b_wb(:,e,k) = Cp*h(Tm_wb(e))*strain_b_wb(:,e,k)/2; %(bending stress on the top surface)
        sigma_plus_wb = [sigma_m_wb(:,e,k) + sigma_b_wb(:,e,k); sigma_s_wb(:,e,k)]'; %(stress on the top surface)
        sigma_VMplus_wb = (sigma_plus_wb(1)^2 + sigma_plus_wb(2)^2 -sigma_plus_wb(1)*sigma_plus_wb(2) + 3*(sigma_plus_wb(3)^2 + sigma_plus_wb(4)^2 + sigma_plus_wb(5)^2))^(1/2);
        sigma_minus_wb = [sigma_m_wb(:,e,k) - sigma_b_wb(:,e,k); sigma_s_wb(:,e,k)]'; %(stress on the bottom surface)
        sigma_VMminus_wb = (sigma_minus_wb(1)^2 + sigma_minus_wb(2)^2 -sigma_minus_wb(1)*sigma_minus_wb(2) + 3*(sigma_minus_wb(3)^2 + sigma_minus_wb(4)^2 + sigma_minus_wb(5)^2))^(1/2);
        sigma_VM_wb(e,k) = max(sigma_VMplus_wb,sigma_VMminus_wb);
    end %loop over Gauss points
end %loop over elements
% 

N_elem_sk = size(Tn_sk,1);
% 7) Postprocess: Computing local strain and stress in shell elements
% 7.1) Get stress and strain at each Gauss point
for e = 1:N_elem_sk
    % a) Get each strain component:
    for j = 1:6
        Idof_sk(j,1)    = 6*(Tn_sk(e,1)-1) + j;
        Idof_sk(6+j,1)  = 6*(Tn_sk(e,2)-1) + j;
        Idof_sk(12+j,1) = 6*(Tn_sk(e,3)-1) + j;
        Idof_sk(18+j,1) = 6*(Tn_sk(e,4)-1) + j;
    end %loop over DOFs

    for k = 1:4
        strain_b_sk(:,e,k)   = [Bb_sk(:,:,e,k)]*[R_sk(:,:,e)]*u_vec(Idof_sk,1);
        strain_m_sk(1:2,e,k) = [Bmn_sk(:,:,e,k)]*[R_sk(:,:,e)]*u_vec(Idof_sk,1);
        strain_m_sk(3,e,k)   = [Bmt_sk(:,:,e)]*[R_sk(:,:,e)]*u_vec(Idof_sk,1);
        strain_s_sk(:,e,k)   = [Bs_sk(:,:,e)]*[R_sk(:,:,e)]*u_vec(Idof_sk,1);
    end

    % b) Get stress
    Cp = [1, nu, 0;
          nu, 1, 0;
          0,  0, (1-nu)/2] *E/(1-nu^2);

    Cs = [1, 0; 
          0, 1] *E/(2*(1+nu));

    for k = 1:4
        sigma_m_sk(:,e,k) = Cp*strain_m_sk(:,e,k); %(constant membrane stress over the thickness)
        sigma_s_sk(:,e,k) = Cs*strain_s_sk(:,e,k); %(constant shear stress over the thickness assumed)
        sigma_b_sk(:,e,k) = Cp*h(Tm_sk(e))*strain_b_sk(:,e,k)/2; %(bending stress on the top surface)
        sigma_plus_sk = [sigma_m_sk(:,e,k) + sigma_b_sk(:,e,k); sigma_s_sk(:,e,k)]'; %(stress on the top surface)
        sigma_VMplus_sk = (sigma_plus_sk(1)^2 + sigma_plus_sk(2)^2 -sigma_plus_sk(1)*sigma_plus_sk(2) + 3*(sigma_plus_sk(3)^2 + sigma_plus_sk(4)^2 + sigma_plus_sk(5)^2))^(1/2);
        sigma_minus_sk = [sigma_m_sk(:,e,k) - sigma_b_sk(:,e,k); sigma_s_sk(:,e,k)]'; %(stress on the bottom surface)
        sigma_VMminus_sk = (sigma_minus_sk(1)^2 + sigma_minus_sk(2)^2 -sigma_minus_sk(1)*sigma_minus_sk(2) + 3*(sigma_minus_sk(3)^2 + sigma_minus_sk(4)^2 + sigma_minus_sk(5)^2))^(1/2);
        sigma_VM_sk(e,k) = max(sigma_VMplus_sk,sigma_VMminus_sk);
    end %loop over Gauss points
end %loop over elements

N_elem_rb = size(Tn_rb,1);
% 7) Postprocess: Computing local strain and stress in shell elements
% 7.1) Get stress and strain at each Gauss point
for e = 1:N_elem_rb
    % a) Get each strain component:
    for j = 1:6
        Idof_rb(j,1)    = 6*(Tn_rb(e,1)-1) + j;
        Idof_rb(6+j,1)  = 6*(Tn_rb(e,2)-1) + j;
        Idof_rb(12+j,1) = 6*(Tn_rb(e,3)-1) + j;
        Idof_rb(18+j,1) = 6*(Tn_rb(e,4)-1) + j;
    end %loop over DOFs

    for k = 1:4
        strain_b_rb(:,e,k)   = [Bb_rb(:,:,e,k)]*[R_rb(:,:,e)]*u_vec(Idof_rb,1);
        strain_m_rb(1:2,e,k) = [Bmn_rb(:,:,e,k)]*[R_rb(:,:,e)]*u_vec(Idof_rb,1);
        strain_m_rb(3,e,k)   = [Bmt_rb(:,:,e)]*[R_rb(:,:,e)]*u_vec(Idof_rb,1);
        strain_s_rb(:,e,k)   = [Bs_rb(:,:,e)]*[R_rb(:,:,e)]*u_vec(Idof_rb,1);
    end

    % b) Get stress
    Cp = [1, nu, 0;
          nu, 1, 0;
          0,  0, (1-nu)/2] *E/(1-nu^2);

    Cs = [1, 0; 
          0, 1] *E/(2*(1+nu));

    for k = 1:4
        sigma_m_rb(:,e,k) = Cp*strain_m_rb(:,e,k); %(constant membrane stress over the thickness)
        sigma_s_rb(:,e,k) = Cs*strain_s_rb(:,e,k); %(constant shear stress over the thickness assumed)
        sigma_b_rb(:,e,k) = Cp*h(Tm_rb(e))*strain_b_rb(:,e,k)/2; %(bending stress on the top surface)
        sigma_plus_rb = [sigma_m_rb(:,e,k) + sigma_b_rb(:,e,k); sigma_s_rb(:,e,k)]'; %(stress on the top surface)
        sigma_VMplus_rb = (sigma_plus_rb(1)^2 + sigma_plus_rb(2)^2 -sigma_plus_rb(1)*sigma_plus_rb(2) + 3*(sigma_plus_rb(3)^2 + sigma_plus_rb(4)^2 + sigma_plus_rb(5)^2))^(1/2);
        sigma_minus_rb = [sigma_m_rb(:,e,k) - sigma_b_rb(:,e,k); sigma_s_rb(:,e,k)]'; %(stress on the bottom surface)
        sigma_VMminus_rb = (sigma_minus_rb(1)^2 + sigma_minus_rb(2)^2 -sigma_minus_rb(1)*sigma_minus_rb(2) + 3*(sigma_minus_rb(3)^2 + sigma_minus_rb(4)^2 + sigma_minus_rb(5)^2))^(1/2);
        sigma_VM_rb(e,k) = max(sigma_VMplus_rb,sigma_VMminus_rb);
    end %loop over Gauss points
end %loop over elements
Imodes=[1,2];

N_mode_reduc = numel(Imodes);            % Number of selected modes
N_cases = size(f_vec, 2);        % Number of forcing vectors

Ureduced = zeros(N_dof, N_cases);    % Initialize output

% --- Reduced-order modeling loop ---
for k = 1:N_cases
    f = f_vec(:, k);                          % Current force vector
    alpha = zeros(N_mode_reduc, 1);               % Modal participation coefficients

    for j = 1:N_mode_reduc
        modeIdx = Imodes(j);                    % Index of the selected mode
        phi_j = Phi(:, modeIdx);            % Mode shape
        lambda_j = lambda(modeIdx);         % Corresponding eigenvalue
        alpha(j) = (phi_j' * f) / lambda_j; % Modal coefficient
        Ureduced(:, k) = Ureduced(:, k) + phi_j * alpha(j);
    end
end
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

save('wing_results.mat');

scale=5;
plotDeformed('wing',xn,Tn_wb,u_vec,scale,sigma_VM_wb); % For wingbox elements
plotDeformed('wing',xn,Tn_rb,u_vec,scale,sigma_VM_rb); % For rib elements
plotDeformed('wing',xn,Tn_sk,u_vec,scale,sigma_VM_sk); % For skin elements


% xSpan = xn(indSpar1,1);     
% figure;hold on
% yyaxis left
% h0 = plot(xSpan, zeros(size(xSpan)), 'k--', 'LineWidth',1);   
% h1 = plot(xSpan, u_z,  'b', 'LineWidth',2);                   
% ylabel('Vertical deflection $u_z$ [m]', 'Interpreter','latex','FontSize',13);
% ylim padded
% ax = gca;                      
% ax.YColor = 'b';
% yyaxis right
% h2 = plot(xSpan, theta_x, 'r', 'LineWidth',2);                
% ylabel('Twist angle $\theta_x$ [rad]','Interpreter','latex','FontSize',13);
% ylim padded
% ax.YColor = 'r';
% xlabel('Spanwise position $x$ [m]', 'Interpreter','latex','FontSize',13);
% title('$u_z$ and $\theta_{x}$ along the span', 'Interpreter', 'latex', 'FontSize', 13, 'FontWeight', 'bold');
% legend([h0 h1 h2], {'Original', '$u_z$ (deformed)', '$\theta_x$ (twist)'}, 'Interpreter','latex','Location','best');
% grid on;  grid minor;  box on;  axis padded;
% xlim([0 5]);
% 
% %Modes
% plotModes('wing',Phi,freq,Imodes)
