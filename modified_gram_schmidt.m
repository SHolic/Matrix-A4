function [Q, R] = modified_gram_schmidt(V)
    % Input: V is an m x n matrix (m vectors, n dimensions)
    % Output: Q is an m x n matrix containing orthonormal vectors
    %         R is an n x n upper triangular matrix
    
    [m, n] = size(V);
    Q = zeros(m, n); % Initialize the matrix for orthonormal vectors
    R = zeros(n, n); % Initialize the upper triangular matrix
    
    for k = 1:n
        % Step 1: Orthogonalize v_k against all previous q_i vectors
        for i = 1:k-1
            R(i, k) = Q(:, i)' * V(:, k);
            V(:, k) = V(:, k) - Q(:, i) * R(i, k);
        end
        
        % Step 2: Normalize the vector to produce q_k
        R(k, k) = norm(V(:, k), 2);
        if R(k, k) == 0
            error('Vectors are linearly dependent');
        end
        Q(:, k) = V(:, k) / R(k, k);
    end
end