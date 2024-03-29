function [Q, R] = reorthogonalization(V)
    % Input: V is an m x n matrix (m vectors, n dimensions)
    % Output: Q is an m x n matrix containing orthonormal vectors
    %         R is an n x n upper triangular matrix
    
    [m, n] = size(V);
    Q = zeros(m, n); % Initialize the matrix for orthonormal vectors
    R = zeros(n, n); % Initialize the upper triangular matrix
    
    for k = 1:n
        % Orthogonalization
        for i = 1:k-1
            R(i, k) = Q(:, i)' * V(:, k);
        end
        for i = 1:k-1
            V(:, k) = V(:, k) - Q(:, i) * R(i, k);
        end
        
        % Reorthogonalization
        S = zeros(1, k-1);
        for i = 1:k-1
            S(i) = Q(:, i)' * V(:, k);
        end
        for i = 1:k-1
            V(:, k) = V(:, k) - Q(:, i) * S(i);
        end

        % If R is needed, update it
        R(1:k-1, k) = R(1:k-1, k) + S';

        % Compute R(k, k) and check for linear dependence
        R(k, k) = norm(V(:, k), 2);
        if R(k, k) == 0
            error('Vectors are linearly dependent');
        end
        
        % Normalize the k-th vector
        Q(:, k) = V(:, k) / R(k, k);
    end
end