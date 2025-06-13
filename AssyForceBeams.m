function [f_vec] = AssyForceBeams(Tn,xn,R,N,Me,B,Q,f_vec)
N_elem = size(Tn,1); %from previously
w = [1;1]; %from previously
%  3.4 Assembly process:
b = zeros(12,N_elem);
q = zeros(12,N_elem);
fe = zeros(12,N_elem);
for e = 1:N_elem
    le = norm(xn(Tn(e,2),:)-xn(Tn(e,1),:)); %from previously
    % a) Compute element force vector:
    b(:,e) = [B(Tn(e,1),:), B(Tn(e,2),:)].';
    q(:,e) = [Q(Tn(e,1),:), Q(Tn(e,2),:)].';
    fe(:,e) = [Me(:,:,e)]*b(:,e);
    for k = 1:2
        fe(:,e) = fe(:,e) + w(k)*le*[R(:,:,e)].'*[N(:,:,e,k)].'*[N(:,:,e,k)]*[R(:,:,e)]*q(:,e)/2;
    end
    % b) Assembly to global force vector:
    for j = 1:6
            Idof(j,1) = 6*(Tn(e,1)-1) + j;
            Idof(6+j,1) = 6*(Tn(e,2)-1) + j;
    end
    f_vec(Idof,1) = f_vec(Idof,1) + fe(:,e);
end
end