clc; clear; close all;

%% PART I - BEAM MODELLING (TEMPLATE)

% Useful commands to initialize script
clear
close all

%% DATA

% Material (same for all beam elements in the wing‑box)
E    = 68.8e9;        % Young's modulus   [Pa]
nu   = 0.33;          % Poisson ratio     [‑]
G    = E/(2*(1+nu));  % Shear modulus     [Pa]
rho  = 2700;          % Density           [kg/m^3]

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

%% PREPROCESS
zeros(0,3);
% Load mesh data
load('beam.mat','xn','Tn','Tm');

% Define boundary conditions and forces data matrices: Up, Fe, Qe, Be
%Se restringe totalmente el nodo 1??
Up=[
   0, 1, 1
   0, 1, 2
   0, 1, 3
   0, 1, 4
   0, 1, 5
   0, 1, 6
];

Fe=[
%-1, length(xn), 3 %1.2   
1, length(xn), 4 %1.3
];

Qe=[];
Be=[];

%%2.1 Inicialization
Nnodes=length(xn); %Cantidad de nodos
Nelem=Nnodes-1;
Ndof=6*Nnodes; %Cantidad de DOF totales(en el sistema)
K=sparse(Ndof, Ndof);
M=sparse(Ndof, Ndof);

%%2.2 Assembly process

l=zeros(1,Nelem);
R=zeros(12, 12);

Ba_prima=zeros(1,12);
Ca_prima=zeros(1, Nelem);
Ka=zeros(12,12);

Bb_prima=zeros(2,12);
Kb=zeros(12,12);

Bs_prima=zeros(2,12);
Ks=zeros(12,12);

Bt_prima=zeros(1,12);
Kt=zeros(12,12);

Idof=zeros(12, 1);

for e=1:Nelem
    %a) Compute rotation matrix
    l(e) = norm(xn(Tn(e,2),:)-xn(Tn(e,1),:)); %Calculo longitud elemento 
    i_prima=(xn(Tn(e,2),:)'-xn(Tn(e,1),:)')/l(e);
    j_prima= [0, 1, 0]';
    k_prima = cross(i_prima, j_prima);
    R3 = [ i_prima,j_prima, k_prima ];
    R_prima = blkdiag(R3, R3)';    
    R(:,:,e)=blkdiag(R_prima, R_prima);

    %b) Compute shape function derivatives
    Nx_prima(1)=-1/l(e);
    Nx_prima(2)=1/l(e);

    %c) Compute each element matrix
        %c1) Axial component of stiffness matrix:
        Ba_prima(1, 1, e)=Nx_prima(1);
        Ba_prima(1, 7, e)=Nx_prima(2);
        Ca_prima(e)=E(Tm(e))*A(Tm(e));%%Es un escalar o un vector??
        Ka(:, :, e)=l(e)* (R(:,:,e))'*(Ba_prima(1, :, e))'*Ca_prima(e)*(Ba_prima(1, :, e))*(R(:,:,e));
        
        %c2) Bending component of stiffness matrix
        Bb_prima(1, 5, e)=Nx_prima(1); Bb_prima(1, 11, e)=Nx_prima(2);
        Bb_prima(2, 6, e)=Nx_prima(1); Bb_prima(2, 12, e)=Nx_prima(2);
        
        Cb_prima=E(Tm(e))*diag([ Iy_prima(Tm(e)) , Iz_prima(Tm(e)) ]);
        Kb(:,:,e)=l(e)*(R(:,:,e))'*(Bb_prima(:, :, e))'*Cb_prima*(Bb_prima(:, :, e))*(R(:,:,e));
        
        %c3) Shear component of stiffness matrix
        N=1/2;
        
        Bs_prima(1, 2, e)=Nx_prima(1); Bs_prima(1, 6, e)=-N; Bs_prima(1, 8, e)=Nx_prima(2); Bs_prima(1, 12, e)=-N;
        Bs_prima(2, 3, e)=Nx_prima(1); Bs_prima(2, 5, e)=+N; Bs_prima(2, 9, e)=Nx_prima(2); Bs_prima(2, 11, e)=+N;

        Cs_prima=G(Tm(e))*A(Tm(e))*diag([ ky_prima(Tm(e)) , kz_prima(Tm(e)) ]);
        Ks(:,:,e)=l(e)*(R(:,:,e))'*(Bs_prima(:, :, e))'*Cs_prima*(Bs_prima(:, :, e))*(R(:,:,e));

        %c4) Torsion component of stiffness matrix
        
        Bt_prima(1, 4, e)=Nx_prima(1);
        Bt_prima(1, 10, e)=Nx_prima(2);

        Ct_prima=G(Tm(e)) * J(Tm(e)) * kt_prima(Tm(e));%%es un escalar o un vector
        Kt(:,:,e)=l(e)*(R(:,:,e))'*(Bt_prima(1, :, e))'*Ct_prima*(Bt_prima(1, :, e))*(R(:,:,e));

        %c5) Mass matrix
        xi=[-(1/sqrt(3)),+(1/sqrt(3))]; %%diferencia entre [] y {}??
        w=[1,1];
        rho_prima=rho(Tm(e))* (diag([A(Tm(e)), A(Tm(e)), A(Tm(e)), J(Tm(e)), Iy_prima(Tm(e)), Iz_prima(Tm(e))]));
        Me(:, :, e)=zeros(12,12);

        for k=1:2
            N(1) = (1 - xi(k))/2;
            N(2) = (1 + xi(k))/2; 

            N_mat(:,:,e,k)=[N(1)*eye(6), N(2)*eye(6)];
            Me(:,:,e)=(Me(:,:,e)+w(k)*l(e)*R(:,:,e)'*N_mat(:,:,e,k)'*rho_prima*N_mat(:,:,e,k)*R(:,:,e))/2;
        end
    %d) Assembly to global matrices
    for j=1:6
        Idof(j,1)=6*(Tn(e,1)-1)+j;
        Idof(6+j,1)=6*(Tn(e,2)-1)+j;
    end

    K(Idof,Idof) = K(Idof,Idof) + Ka(:,:,e)+Kb(:,:,e)+Ks(:,:,e)+Kt(:,:,e);
    M(Idof,Idof) = M(Idof,Idof) + Me(:,:,e);
end

%% SOLVER

% Obtain system matrices

%% ------------------- GUARDAR / CARGAR MATRICES ------------------------
recompute = true;   % ← pon TRUE si cambias la malla o propiedades

matFile = 'beam_matrices.mat';

if ~recompute && isfile(matFile)
    % ----------  Cargar matrices ya calculadas  ------------------------
    load(matFile , 'K' , 'M');
    fprintf('[INFO] Matrices K y M cargadas desde %s\n', matFile);

else
    % ----------  Recalcular y guardar  ---------------------------------
    % (en este punto K, M, ell y R ya los has creado en el bucle anterior)

    save(matFile , 'K' , 'M');
    fprintf('[INFO] Matrices K y M guardadas en %s\n', matFile);
end

%% ---------------------------------------------------------------------

%3) Compute global force vector
    %3.1) Point loads
    f=zeros(Ndof,1);
    
    for q = 1:size(Fe,1)    % recorre cada fila de la matriz Fe
        f(6*(Fe(q,2)-1) + Fe(q,3),1) = f(6*(Fe(q,2)-1)+Fe(q,3),1)+Fe(q,1);
    end

    %3.2) Nodal distributed forces
    Q=zeros(Nnodes,6);

    for r = 1:size(Qe,1)
        Q(Qe(r,2),Qe(r,3))=Q(Qe(r,2),Qe(r,3))+Qe(r,1);
    end

    %3.3) Nodal body forces

    B=zeros(Nnodes,6);

    for s=1:size(Be,1)
        B(Be(s,2),Be(s,3))=B(Be(s,2),Be(s,3))+Be(s,1);
    end
    
    %3.4) Assembly process

    b = zeros(12 , Nelem);
    q = zeros(12 , Nelem);
    fe= zeros(12 , Nelem);

    for e=1:Nelem

        %a)Compute element force vector
        b(:,e)=[B(Tn(e,1),:), B(Tn(e,2),:)]';
        q(:,e)=[Q(Tn(e,1),:), Q(Tn(e,2),:)]';
        fe(:,e)=Me(:,:,e)*b(:,e);
        
        for k=1:2
            fe(:,e)=fe(:,e)+w(k)*l(e)*R(:,:,e)'*N_mat(:,:,e,k)'*N_mat(:,:,e,k)*R(:,:,e)*q(:,e)/2;
        end

        %b)Assembly to global force vector

        for j=1:6
            Idof(j,1)=6*(Tn(e,1)-1)+j;
            Idof(6+j,1)=6*(Tn(e,2)-1)+j;
        end
        f(Idof,1)=f(Idof,1)+fe(:,e);
    end

%4) Boundary conditions

    %4.1) Initialization
    u=zeros(Ndof,1);
    Ip=zeros(1,6);

    %4.2) Prescribed and free DOF's

    for p = 1:size(Up,1)
        Ip(p)=6*(Up(p,2)-1)+Up(p,3);
        u(Ip(p),1)=Up(p,1);
    end
    
    If = setdiff(1:Ndof,Ip);
   
%5)Solve system of equations(static case)

    %5.1)Solve system

    u(If,1)=K(If,If) \ (f(If,1)-K(If,Ip)*u(Ip,1)); %\Calcula la inversa
    fR= K*u-f;

%6)Postprocess: Computing strain and internal forces in beam elements

epsilon_a=zeros(1,Nelem);
epsilon_s=zeros(2,Nelem);
epsilon_t=zeros(1,Nelem);
epsilon_b=zeros(2,Nelem);

fint_prima=zeros(12,Nelem);
Fx_prima=zeros(1,Nelem);
Fy_prima=zeros(1,Nelem);
Fz_prima=zeros(1,Nelem);
Mx_prima=zeros(1,Nelem);
My_prima=zeros(1,Nelem);
Mz_prima=zeros(1,Nelem);

    %6.1)Get strain at each beam element

    for e=1:Nelem

        %a)Get element displacements

        for j=1:6
            
            Idof(j,1)=6*(Tn(e,1)-1)+j;
            Idof(6+j,1)=6*(Tn(e,2)-1)+j;
        end
        ue=u(Idof,1);

        %b) Get each strain component

        epsilon_a(1,e)=Ba_prima(1,:,e)*R(:,:,e)*ue;
        epsilon_s(:,e)=Bs_prima(:,:,e)*R(:,:,e)*ue;
        epsilon_t(1,e)=Bt_prima(1,:,e)*R(:,:,e)*ue;
        epsilon_b(:,e)=Bb_prima(:,:,e)*R(:,:,e)*ue;

        %b) Get each strain component

        fint_prima(:,e)=R(:,:,e)*(Ka(:,:,e)+Kb(:,:,e)+Ks(:,:,e)+Kt(:,:,e))*ue;
        Fx_prima(:,e)=fint_prima(1,e);-fint_prima(7,e);
        Fy_prima(:,e)=fint_prima(2,e);-fint_prima(8,e);
        Fz_prima(:,e)=fint_prima(3,e);-fint_prima(9,e);
        Mx_prima(:,e)=fint_prima(4,e);-fint_prima(10,e);
        My_prima(:,e)=fint_prima(5,e);-fint_prima(11,e);
        Mz_prima(:,e)=fint_prima(6,e);-fint_prima(12,e);
    end

% Perform modal analysis
% ...

% Compute external forces vector
% ...

% Solve system
% ...

%% POSTPROCESS

% Save data for postprocessing in separate script file (useful when results
% from different runs need to be compared)
save('beam_results.mat');

%% --- 1 POST-PROCESS : dibujo de la viga original y deformada -------------
scale = 100000;          % factor de ampliación de la deformada  (ajusta a ojo)
Nnodes = size(xn,1);
% Coordenadas globales
x = xn(:,1);         % eje x global (longitud de la viga)
y = xn(:,2);         % eje y  (¡casi 0 en tu malla!)
z = xn(:,3);         % eje z  (corda)         
% Desplazamientos verticales (u_z local = DOF 3    de cada nodo)
uz = u(3:6:end);     % toma 3,9,15,…
% Coordenadas deformadas
z_def = z + scale*uz;   % solo despl. vertical; añade otros si quieres

x = xn(:,1);               % Coordenada x (envergadura)
uz = u(3:6:end);           % Desplazamiento vertical (DOF 3)
theta_x = u(4:6:end);      % Rotación en torno a x (torsión, DOF 4)

figure; grid on; hold on;

yyaxis left
hh2 = plot(x, 1e5 * uz, 'b-', 'LineWidth', 2);                % deformada
hh1 = plot(x , z ,'--','Color',[0.6 0.6 0.6],'LineWidth',1.5);  % original
ylabel('Vertical deflection u_z [m]')
ylim([-0.02 0.08])

yyaxis right
hh3 = plot(x, theta_x * 1e5, 'r-', 'LineWidth', 2);           % torsión
ylabel('Twist angle \theta_x [rad]')

xlabel('Spanwise position x [m]')
title('Vertical deflection and twist angle along the span')
legend([hh1 hh2 hh3], 'Original', 'u_z (deformed)', '\theta_x (twist)')  % <- CORRECTO

xlim([0 5])
ylim([-0.02 0.08])

%% --- 2 POST-PROCESS : modes-------------

F = zeros(Ndof,1);              
F(6*(Nnodes-1)+3) = -1;          

[U, lambda, Phi] = frequency(F, Up);



