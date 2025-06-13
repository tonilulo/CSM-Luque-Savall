function [f_vec] = AssyForceShells(Tn,B,P,R,N,Me,S_4,f_vec)
N_elem = size(Tn,1);
% 4.4) Assembly process
for e = 1:N_elem
    % a) Compute element force vector:
    b(:,e) = [B(Tn(e,1),:),B(Tn(e,2),:),B(Tn(e,3),:),B(Tn(e,4),:)]';
    p(:,e) = [P(Tn(e,1),:),P(Tn(e,2),:),P(Tn(e,3),:),P(Tn(e,4),:)]';
    fe(:,e) = [Me(:,:,e)]*b(:,e);

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
end %loop over elements
end %function