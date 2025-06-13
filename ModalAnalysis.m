function [Phi, lambda, freq] = ModalAnalysis(K, M, I_f, N_m, N_dof)
    % Symmetrize and sparsify stiffness and mass matrices
    K = sparse(0.5 * (K + K.'));
    M = sparse(0.5 * (M + M.'));

    % Perform eigenvalue analysis
    [V, D] = eigs(K(I_f, I_f), M(I_f, I_f), N_m, 'sm');

    % Preallocate outputs
    Phi = zeros(N_dof, N_m);
    lambda = zeros(1,N_m);                
    freq = zeros(N_m,1);  

    % Normalize eigenvectors and populate global mode shapes
    Minv = M(I_f, I_f);                  
    for k = 1:size(V, 2)
        v = V(:, k);
        Phi(I_f, k) = v / sqrt(v' * Minv * v);
        lambda(k) = D(k, k);
        freq(k) = sqrt(lambda(k)) / (2 * pi);
    end
end
