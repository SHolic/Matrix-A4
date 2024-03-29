function [Q, R] = classical_gram_schmidt(V)
    % Input: V is an m x n matrix (m vectors, n dimensions)
    % Output: Q is an m x n matrix containing orthonormal vectors
    %         R is an n x n upper triangular matrix

    [m, n] = size(V);
    Q = zeros(m, n); % Initialize the matrix for orthonormal vectors
    R = zeros(n, n); % Initialize the upper triangular matrix

    for k = 1:n
        v_k = V(:, k);
        for i = 1:k-1
            R(i, k) = Q(:, i)' * v_k;
        end
        for i = 1:k-1
            v_k = v_k - Q(:, i) * R(i, k);
        end

        R(k, k) = norm(v_k, 2);
        if R(k, k) == 0
            error('Vectors are linearly dependent');
        end
        Q(:, k) = v_k / R(k, k);
    end
end