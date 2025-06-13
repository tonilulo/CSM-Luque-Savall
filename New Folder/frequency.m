function [U, lambda, Phi] = frequency(F, Up)

%1.1 Input data

Nm=6;
load('beam_matrices.mat','K','M');
load('beam.mat','xn');

% (Opcional) garantizar que siguen siendo matrices dispersas
K = sparse(K);
M = sparse(M);

Ndof=sqrt(numel(K));
Nnodes=Ndof/6;


% F=[
% -1, length(Ndof/6), 3 %1.2   
% % 1, length(Ndof/6), 4 %1.3
% ];

% F = zeros(Ndof,1);          % una sola columna = un solo caso de carga
% F(6*(Nnodes-1)+3) = -1;     % DOF 3 del último nodo


w=0;

% Up=[
%    0, 1, 1
%    0, 1, 2
%    0, 1, 3
%    0, 1, 4
%    0, 1, 5
%    0, 1, 6
% ];

%2) Boundary conditions
Nw=numel(w);

U=zeros(Ndof,Nw);

for p = 1:size(Up,1)
    Ip(p)=6*(Up(p,2)-1)+Up(p,3);
    u(Ip(p),:)=Up(p,1);
end

If = setdiff(1:Ndof,Ip);

%3)Solve system of equations

for k = 1:size(F,2)
U(If,k) = (K(If,If) - w(k)^2 * M(If,If)) \(F(If,k) - (K(If,Ip) - w(k)^2 * M(If,Ip)) * U(Ip,k));
end

%4)Modal analisys

K=0.5*(K+K');
M=0.5*(M+M');
[V,D]=eigs(K(If, If), M(If, If), Nm,'sm');

Phi=zeros(Ndof,Nm);
lambda=zeros(1,Nm);

for k = 1:size(V,2)
    Phi(If,k) = V(:,k) / sqrt( V(:,k).' * M(If,If) * V(:,k) );
    lambda(k) = D(k,k);
end

%5)Model order reduction

Im = 1:Nm;

Ustar=zeros(Ndof,Nw);

for k = 1:size(F,2)                     

    Ustar(:,k) = zeros(size(F,1),1);      

    for j = 1:length(Im)                             
        alpha(j,k) = ( Phi(:,Im(j)).' * F(:,k) ) / ( lambda(j) - w(k)^2 );
        Ustar(:,k) = Ustar(:,k) + Phi(:,j) * alpha(j,k);
    end
end
%Note: The model-order reduction procedure is also applicable to static problems, just by
%setting w(k) = 0 and making F(:,k) the corresponding (static) force vector.

x = xn(:,1);

uy = Phi(2:6:end, 1);       % DOF 2 (desplazamiento en y)
uz = Phi(3:6:end, 1);       % DOF 3 (desplazamiento en z)
thetax = Phi(4:6:end, 1);   % DOF 4 (rotación torsional)

for i = 1:Nm
    
    figure; 
    subplot(3,1,1)
    plot(x, Phi(2:6:end, i)); title(['u_y - Mode ', num2str(i)]); ylabel('u_y'); grid on

    subplot(3,1,2)
    plot(x, Phi(3:6:end, i)); title(['u_z - Mode ', num2str(i)]); ylabel('u_z'); grid on

    subplot(3,1,3)
    plot(x, Phi(4:6:end, i)); title(['\theta_x - Mode ', num2str(i)]); ylabel('\theta_x'); xlabel('x [m]'); grid on
end

fprintf('\nPrimeros %d modos naturales:\n', Nm);
fprintf('Modo\t lambda (rad^2/s^2)\t omega (rad/s)\t f (Hz)\n');
for i = 1:Nm
    omega_i = sqrt(lambda(i));
    f_i = omega_i / (2*pi);
    fprintf('%d\t %.4e\t %.4f\t %.4f\n', i, lambda(i), omega_i, f_i);
end

end

