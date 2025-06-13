function [sigma_VM] = AssyStressShells(Tn,R,u_vec,Bb,Bmn,Bmt,Bs,nu,E,h,Tm)
N_elem = size(Tn,1);
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
end %function