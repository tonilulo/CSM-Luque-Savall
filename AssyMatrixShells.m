function [K_up,M,N,R,M_e,S_4,Bb,Bmn,Bmt,Bs] = AssyMatrixShells(Tm,Tn,xn,h,E,nu,rho,Tnb)
% 2.1) Initialization:
N_elem = numel(Tm);
N_nod = size(xn,1);
N_dof = 6*N_nod; %(total number of degrees of freedom)
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
K_up = K + Kr;
end