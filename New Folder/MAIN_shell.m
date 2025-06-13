%% PART II - SHELL MODELLING (TEMPLATE)

clear
close all

%% DATA

% Material (same for all beam elements in the wing‑box)
E    = [68.8e9, 68.8e9, 68.8e9];        % Young's modulus   [Pa]
nu   = [0.33, 0.33, 0.33];          % Poisson ratio     [‑]
G    = E/(2*(1+nu));  % Shear modulus     [Pa]
rho  = [2700, 2700, 2700];          % Density           [kg/m^3]

% Cross‑sectional properties of the wing‑box (NACA0012 – Table 1 in project)
A  = 0.022;           % Area                [m^2]
Iy_prima = 1.097e-4;        % I_y' (about local y') [m^4]
Iz_prima = 12.087e-4;       % I_z' (about local z') [m^4]
J  = 13.184e-4;       % Polar inertia        [m^4]
ky_prima = 0.5975;          % Shear correction factor along y'
kz_prima = 0.1375;          % Shear correction factor along z'
kt_prima = 0.2564;          % Torsion correction factor

% Geometric data (needed to locate the shear centre & apply loads)
y0 = 0.6226;    %Inertia I_z
y1 = 0.345;    % front spar           [m] (see Table 1)
y2 = 0.960;    % rear  spar           [m]
yc = 0.5791;   % shear centre y‑coord [m]
c=1.5;
b=5;
h1=0.030; %h1 [m]
h2=0.025;    %h2 [m]
h3=0.011;    %h3 [m]
h = [h1 , h2 , h3];

%% PREPROCESS

% Load mesh data
load('shell.mat','xn','Tn','Tm','indRoot','indPointA','indPointB','indSpar1','indSpar2');

% Define boundary conditions and forces data matrices: Up, Fe, Pe, Be

Up = zeros(6*numel(indRoot) , 3);        % pre-asigna Up con 6 filas por cada nodo
rowUp = 0;
for k = 1:numel(indRoot)                 % recorre los nodos de la raíz
    np = indRoot(k);         % nodo actual
    for jp = 1:6             % recorre los 6 DOF de ese nodo
        rowUp = rowUp + 1;      % avanza una fila en Up
        Up(rowUp,:) = [ 0 , np , jp ];   
    end
end

Fz=-1;
FA=Fz*(y2-yc)/(y2-y1);
FB=Fz*(yc-y1)/(y2-y1);

Fe=[
   -FA, indPointA, 3 %FA
   -FB, indPointB, 3 %FB
];

Pe=[];
Be=[];


%2)Assembly of global matrices:

    %2.1)Inicialization
    Nnodes=length(xn);
    Nelem=numel(Tm);
    Ndof=6*Nnodes;

    K=sparse(Ndof,Ndof);
    M=sparse(Ndof,Ndof);

    %2.2)Assembly process

    for e=1:Nelem

        %a)Compute rotation matrix
        S=cross(xn(Tn(e,3),:)'-xn(Tn(e,1),:)',0.5*(xn(Tn(e,4),:)'-xn(Tn(e,2),:)'));
        k_prima = S / norm(S);
        d=0.5*(xn(Tn(e,2),:)'+ xn(Tn(e,3),:)'-xn(Tn(e,4),:)'-xn(Tn(e,1),:)');
        i_prima=d / norm(d);
        j_prima=cross(k_prima, i_prima);
        R3 = [i_prima, j_prima, k_prima ];
        R2 = [i_prima, j_prima];
        R_prima = blkdiag(R3, R2)';  
        R(:,:,e)=blkdiag(R_prima, R_prima,R_prima,R_prima);
    
       
        %b)Get nodal coefficients for the shape functions
        
        a=[-1,1,1,-1];
        b=[-1,-1,1,1];

        %c)Compute element matrices
            
            %c1)1 Gauss point quadrature matrices

            N1=0.25*[1,1,1,1]';
            dN1dxi=a/4;
            dN1deta=b/4;
            J1=zeros(2,2);

            for i=1:4
                J1=J1+[dN1dxi(i);dN1deta(i)]*xn(Tn(e,i),:)*[i_prima, j_prima];
            end

            dN1dx_prima = J1 \ [ dN1dxi ; dN1deta ];
            S1=4*det(J1);

                %c1.1 Shear component of stiffness matrix:

                for i=1:4
                    Bs_priman(:,:,i)=[0, 0, dN1dx_prima(1,i), 0, N1(i);
                                     0, 0, dN1dx_prima(2,i), -N1(i), 0];
                end
    
                Cs_prima=(eye(2)*5*h(Tm(e))*E(Tm(e)))/(12*(1+nu(Tm(e))));
                Bs_primae(:,:,e)=[Bs_priman(:,:,1), Bs_priman(:,:,2), Bs_priman(:,:,3), Bs_priman(:,:,4)];
                Ks(:,:,e) = S1 * R(:,:,e)' * Bs_primae(:,:,e)' * Cs_prima * Bs_primae(:,:,e) * R(:,:,e);

                %c1.2) Membrane transverse component of stiffness matrix:

                for i = 1:4

                    B_mt_prima(:,:,i) = [dN1dx_prima(2,i) dN1dx_prima(1,i) 0 0 0 ]; 

                end
                
                h_e  = h(Tm(e));                    
                E_e  = E(Tm(e));                   
                nu_e = nu(Tm(e));                 
                C_mt = (h_e * E_e) / (2*(1 + nu_e)); 
                
                B_mte(:,:,e) = [ B_mt_prima(:,:,1), B_mt_prima(:,:,2), B_mt_prima(:,:,3), B_mt_prima(:,:,4) ];  
                Km(:,:,e) = S1 * R(:,:,e)' * B_mte(:,:,e)' * C_mt * B_mte(:,:,e) * R(:,:,e);

                %c2) 4 Gauss points quadrature matrices:
              Kb(:,:,e)=zeros(24,24);
              Me(:,:,e)=zeros(24,24);
              xi4=[-1,1,1,-1]/sqrt(3);
              eta4=[-1,-1,1,1]/sqrt(3);
              w4=[1,1,1,1];

              for k=1:4
                  J4=zeros(2,2);

                  for i=1:4
                    N4=(1+a(i)*xi4(i))*(1+b(i)*eta4)/4;
                    N4xi=a(i)*(1+b(i)*eta4)/4;
                    N4eta=b(i)*(1+a(i)*xi4)/4;

                    J4=J4+[N4xi(i);N4eta(i)]*(xn(Tn(e,i),:))*[i_prima j_prima];
                  end
                N4xprima=J4^-1*[N4xi;N4eta];
                S4(e,k)=w4(k)*det(J4);

                %c2.1) Membrane normal component of stiffness matrix:

                for i=1:4
                        Bmn_prima1(:,:,i) = [dN1dx_prima(1,i), 0,  0,  0,  0 ;
                            0,   dN1dx_prima(2,i),   0,  0,  0 ];
                end
                Cmn_prima=[1 nu(Tm(e));nu(Tm(e)) 1]*h(Tm(e))*E(Tm(e))/(1-nu(Tm(e))^2);
                Bmn_prima(:,:,e,k)=[Bmn_prima1(:,:,1),Bmn_prima1(:,:,2),Bmn_prima1(:,:,3),Bmn_prima1(:,:,4)];
                Km(:,:,e)= Km(:,:,e)+S4(e,k)*R(:,:,e)'*Bmn_prima(:,:,e,k)'*Cmn_prima*Bmn_prima(:,:,e,k)*R(:,:,e);
              

              %c2.2) Bending component of stiffness matrix: 

              for i=1:4 
                Bb_prima1(:,:,i)=[0 0 0 0 N4xprima(1,i);
                                  0 0 0 N4xprima(2,i) 0;
                                  0 0 0 -N4xprima(1,i) N4xprima(2,i)
                                  ];
              end

              Cb_prima1=[1 nu(Tm(e)) 0;
                         nu(Tm(e)) 1 0;
                         0 0 0.5*(1-nu(Tm(e)))];
              Cb_prima=(Cb_prima1*h(Tm(e))^3*E(Tm(e)))/(12*(1-nu(Tm(e))^2));
              
              Bb_prima(:,:,e,k)  =[Bb_prima1(:,:,1), Bb_prima1(:,:,2), Bb_prima1(:,:,3), Bb_prima1(:,:,4)];
              Kb(:,:,e)=Kb(:,:,e)+S4(e,k)*R(:,:,e)'*Bb_prima(:,:,e,k)'*Cb_prima*Bb_prima(:,:,e,k)*R(:,:,e);
              %c2.3 Mass matrix
              for i=1:4 
                  Ni(:,:,i)=N4(i)*eye(5);
              end
        
              rho_prima1=[1 0 0 0 0;
                          0 1 0 0 0;
                          0 0 1 0 0;
                          0 0 0 (1/12)*h(Tm(e))^2 0;
                          0 0 0 0 (1/12)*h(Tm(e))^2];
              rho_prima=rho(Tm(e))*h(Tm(e))*rho_prima1;
                
              N(:,:,e,k)=[Ni(:,:,1),Ni(:,:,2),Ni(:,:,3),Ni(:,:,4)];
              Me(:,:,e)= Me(:,:,e)+S4(e,k)*R(:,:,e)'* N(:,:,e,k)'*rho_prima* N(:,:,e,k)*R(:,:,e);
              end            
%d)Assembly to global matrices:

    for j=1:6
        Idof(j,1)=6*(Tn(e,1)-1)+j;
        Idof(6+j,1)=6*(Tn(e,2)-1)+j;
        Idof(12+j,1)=6*(Tn(e,3)-1)+j;
        Idof(18+j,1)=6*(Tn(e,4)-1)+j;
    end

K(Idof,Idof)=K(Idof,Idof)+Km(:,:,e)+Kb(:,:,e)+Ks(:,:,e);
M(Idof,Idof)=M(Idof,Idof)+Me(:,:,e);
end

%3)Compute artificial rotation stiffness matrix
%3.1 Find nodal normal to set criteria for finding coplanar nodes:
n=zeros(3,Nnodes);

for e=1:Nelem
    %a) Compute normal and surface
    Sf1 = 0.5 * cross( ...
       xn(Tn(e,3),:).' - xn(Tn(e,1),:).', ...  % primer vector (3×1)
       xn(Tn(e,4),:).' - xn(Tn(e,2),:).'  ...  % segundo vector (3×1)
     );
    Sf(e)=sqrt(Sf1(1)^2+Sf1(2)^2+Sf1(3)^2);
    k_prima(:,e)=Sf1/Sf(e);

    %b)Assemble to get nodal normal:

    for i=1:4 %esto es 4??
        n(:,Tn(e,i))=n(:,Tn(e,i))+k_prima(:,e);
    end    
end

%3.2)Compute artificial rotation matrix:

Kr=zeros(Ndof,Ndof);

for e=1:Nelem
    
    for i=1:4
        
        %a)Determine whether it is or not a coplanar node
        alpha = acos(dot(n(:,Tn(e,i)), k_prima(:,e))/norm(n(:,Tn(e,i))));
        if alpha < deg2rad(5)
        end
        %b) Evaluate artificial rotation stiffness component
        Idof=6*(Tn(e,i)-1)+[4,5,6]';
        Kr(Idof,Idof)= Kr(Idof,Idof)+E(Tm(e))*h(Tm(e))*Sf(e)*k_prima(:,e)*k_prima(:,e)';
        
    end

end

%3.3)Update stiffness matrix:

K=K+Kr;

%4)Compute global force vector:

    f=zeros(Ndof,1);

    %4.1) Point loads
    for q = 1:size(Fe,1)    % recorre cada fila de la matriz Fe
        f(6*(Fe(q,2)-1) + Fe(q,3),1) = f(6*(Fe(q,2)-1)+Fe(q,3),1)+Fe(q,1);
    end

    %4.2) Nodal distributed forces
    P=zeros(Nnodes,6);

    for r = 1:size(Pe,1)
        P(Pe(r,2),Pe(r,3))=P(Pe(r,2),Pe(r,3))+Pe(r,1);
    end

    %4.3) Nodal body forces

    B=zeros(Nnodes,6);

    for s=1:size(Be,1)
        B(Be(s,2),Be(s,3))=B(Be(s,2),Be(s,3))+Be(s,1);
    end
    
    %4.4) Assembly process

    b = zeros(24 , Nelem);
    p = zeros(24 , Nelem);
    fe= zeros(24 , Nelem);

    for e=1:Nelem

        %a)Compute element force vector
        b(:,e)=[B(Tn(e,1),:), B(Tn(e,2),:), B(Tn(e,3),:), B(Tn(e,4),:)]';
        p(:,e)=[P(Tn(e,1),:), P(Tn(e,2),:), P(Tn(e,3),:), P(Tn(e,4),:)]';
        fe(:,e)=Me(:,:,e)*b(:,e);
        
        for k=1:4
            fe(:,e)=fe(:,e)+S4(e,k)*R(:,:,e)'*N(:,:,e,k)'*N(:,:,e,k)*R(:,:,e)*p(:,e);
        end

        %b)Assembly to global force vector

        for j=1:6
            Idof(j,1)=6*(Tn(e,1)-1)+j;
            Idof(6+j,1)=6*(Tn(e,2)-1)+j;
            Idof(12+j,1)=6*(Tn(e,3)-1)+j;
            Idof(18+j,1)=6*(Tn(e,4)-1)+j;
        end
        f(Idof,1)=f(Idof,1)+fe(:,e);
    end
 %5) Boundary conditions
 %5.1) Initialization
    u=zeros(Ndof,1);

    %4.2) Prescribed and free DOF's

    for p = 1:size(Up,1)
        Ip(p)=6*(Up(p,2)-1)+Up(p,3);
        u(Ip(p),1)=Up(p,1);
    end
    
    If = setdiff(1:Ndof,Ip);

  %6) Solve system of equation (static case)
    u(If,1)=K(If,If) \ (f(If,1)-K(If,Ip)*u(Ip,1)); %\Calcula la inversa
    fR= K*u-f;
 %7) Postprocess: Computing local strain and stress in shell elements

% epsilon_m_prima=zeros(3,5,4); 

for e=1:Nelem
    
    for j=1:6
         Idof(j,1)=6*(Tn(e,1)-1)+j;
         Idof(6+j,1)=6*(Tn(e,2)-1)+j;
         Idof(12+j,1)=6*(Tn(e,3)-1)+j;
         Idof(18+j,1)=6*(Tn(e,4)-1)+j; 
    end
    for k=1:4
        epsilon_b_prima(:,e,k)=Bb_prima(:,:,e,k)*R(:,:,e)*u(Idof,1);
        epsilon_m_prima(1:2,e,k)=Bmn_prima(:,:,e,k)*R(:,:,e)*u(Idof,1);
        epsilon_m_prima(3,e,k)=B_mte(:,:,e)*R(:,:,e)*u(Idof,1); 
        epsilon_s_prima(:,e,k)=Bs_primae(:,:,e)*R(:,:,e)*u(Idof,1);
    end

%b)Get stress:

C_p1=[1 nu(Tm(e)) 0;
     nu(Tm(e)) 1 0;
     0 0 0.5*(1-nu(Tm(e)))];
C_p=(C_p1*E(Tm(e)))/(1-nu(Tm(e))^2);

C_s=(eye(2)*E(Tm(e)))/(2*(1+nu(Tm(e))));
 
for k=1:4

    sigma_m_prima(:,e,k)=C_p*epsilon_m_prima(:,e,k);
    sigma_s_prima(:,e,k)=C_s*epsilon_s_prima(:,e,k);
    sigma_b_prima(:,e,k)=C_p*0.5*(h(Tm(e))*epsilon_b_prima(:,e,k));
    sigma_plus_prima=[sigma_m_prima(:,e,k)+sigma_b_prima(:,e,k); sigma_s_prima(:,e,k)]';
    sigma_plus_VM=(sigma_plus_prima(1)^2-sigma_plus_prima(2)^2-sigma_plus_prima(1)*sigma_plus_prima(2)+3*(sigma_plus_prima(3)^2+sigma_plus_prima(4)^2+sigma_plus_prima(5)^2))^(0.5);
    sigma_minus_prima=[sigma_m_prima(:,e,k)-sigma_b_prima(:,e,k);sigma_s_prima(:,e,k)]';
    sigma_minus_VM=(sigma_minus_prima(1)^2-sigma_minus_prima(2)^2-sigma_minus_prima(1)*sigma_minus_prima(2)+3*(sigma_minus_prima(3)^2+sigma_minus_prima(4)^2+sigma_minus_prima(5)^2))^(0.5);
    sigma_VM(e,k)=max(sigma_plus_VM,sigma_minus_VM);
end 

end

scale=10^6;
plotDeformed('shell',xn,Tn,u,scale);

%plotModes('shell',Phi,freq,imodes)
% % This function plots the specified modes resulting from a modal analysis
% % in sets of 9.
% % Phi : Modal displacements matrix in which each column corresponds to the
% %       mode shape of the corresponding mode. Expected dimensions [Ndof x Nmodes]
% % freq : Natural frequencies array. Expected dimensions [Nmodes x 1]
% % imodes : Array selecting which modes to plot. A maximum of 9 modes must
% %          can be selected. Example: imodes = [1,2,3,4,5,6,7,8,9] will plot
% %          the modes stored in the first 9 columns of Phi / imodes = [1,4,5,10] 
% %          will plot modes in columns 1, 4, 5 and 10 of Phi. 