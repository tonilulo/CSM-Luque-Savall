function [K,M,N,R,Me,Ba,Bs,Bt,Bb,Ka,Kb,Ks,Kt,le] = AssyMatrixBeams(Tm,Tn,xn,E,rho,I_y,I_z,J,G,A,k_y,k_z,k_t,j_vec)

% Initialization
N_elem   = numel(Tm);
N_nod    = size(xn,1);
N_dof    = 6*N_nod;           % total number of degrees of freedom
K       = sparse(N_dof,N_dof); %using sparse
M       = sparse(N_dof,N_dof);

%Matrix definitions
R_mat   = zeros(6,6);
R       = zeros(12,12,N_elem);
Ba      = zeros(1,12,N_elem);
Bb      = zeros(2,12,N_elem);
Bs      = zeros(2,12,N_elem);
Bt      = zeros(1,12,N_elem);

%2.2 Assembly process:
for e = 1:N_elem
  
 % a) Compute rotation matrix:
    le = norm(xn(Tn(e,2),:)-xn(Tn(e,1),:)); % (element size = beam length) [m]
    i_vec = (xn(Tn(e,2),:)'-xn(Tn(e,1),:)')/le;
    %j_vec defined already [0, 1, 0]
    k_vec = cross(i_vec,j_vec);
    R_mat(1:3,1:3) = [i_vec, j_vec, k_vec];
    R_mat(4:6,4:6) = [i_vec, j_vec, k_vec];
    R_mat = R_mat';
    R(:,:,e) = [R_mat, zeros(6,6); zeros(6,6), R_mat];

  % b) Compute shape function derivatives:
    N_x(1) = -1/le;
    N_x(2) = 1/le;

  % c) Compute each element matrix:
    % c1) Axial component of stiffness matrix:
    Ba(1,:,e) = [N_x(1), 0, 0, 0, 0, 0, N_x(2), 0, 0, 0, 0, 0];
    Ca = E*A;
    Ka(:,:,e) = le*[R(:,:,e)]'*[Ba(1,:,e)]'*Ca*[Ba(1,:,e)]*[R(:,:,e)];

    % c2) Bending component of stiffness matrix:
    Bb(:,:,e) = [0, 0, 0, 0, N_x(1), 0, 0, 0, 0, 0, N_x(2), 0;
                 0, 0, 0, 0, 0, N_x(1), 0, 0, 0, 0, 0, N_x(2)];
    Cb = E* [I_y, 0; 0, I_z];
    Kb(:,:,e) = le*[R(:,:,e)]'*[Bb(:,:,e)]'*Cb*[Bb(:,:,e)]*[R(:,:,e)];

    % c3) Shear component of stiffness matrix
    N_i = 1/2; %(shape functions assuming only 1 Gauss point)
    Bs(:,:,e) = [0, N_x(1), 0, 0, 0, -N_i, 0, N_x(2), 0, 0, 0, -N_i;
                 0, 0, N_x(1), 0, N_i, 0, 0, 0, N_x(2), 0, N_i, 0 ];
    Cs = G*A*[k_y, 0; 0, k_z];
    Ks(:,:,e) = le*[R(:,:,e)]'*[Bs(:,:,e)]'*Cs*[Bs(:,:,e)]*[R(:,:,e)];

    % c4) Torsion component of stiffness matrix:
    Bt(1,:,e) = [0, 0, 0, N_x(1), 0, 0, 0, 0, 0, N_x(2), 0, 0];
    Ct = G*J*k_t;
    Kt(:,:,e) = le*[R(:,:,e)]'*[Bt(1,:,e)]'*Ct*[Bt(1,:,e)]*[R(:,:,e)];

    % c5) Mass matrix:
    ksi = [-1/sqrt(3);1/sqrt(3)]; %Gauss point coordinates
    w = [1;1];
    rho_mat = rho*diag([A, A, A, J, I_y, I_z]);
    Me(:,:,e) = zeros(12,12);
    for k = 1:2
        N_i(1) = (1-ksi(k))/2;
        N_i(2) = (1+ksi(k))/2;
        N(:,:,e,k) = [N_i(1)*diag([ones(1,6)]), N_i(2)*diag([ones(1,6)])];
        Me(:,:,e) = Me(:,:,e) + w(k)*le*[R(:,:,e)]'*[N(:,:,e,k)]'*rho_mat*[N(:,:,e,k)]*[R(:,:,e)]/2;
    end

    % d) Assembly to global matrices
    Idof = zeros(12,1);
    for j = 1:6
        Idof(j,1) = 6*(Tn(e,1)-1) + j;
        Idof(6+j,1) = 6*(Tn(e,2)-1) + j;
    end

    K(Idof,Idof) = K(Idof,Idof) + Ka(:,:,e) + Kb(:,:,e) + Ks(:,:,e) + Kt(:,:,e);
    M(Idof,Idof) = M(Idof,Idof) + Me(:,:,e);

end %loop over elements

end %function