function [sol, realisable] = heuristique_regret(m, n, c, A, b)
    % get A and c as matrix
    A_mat = vector_to_matrix(m, n, A);
    c_mat = vector_to_matrix(m, n, c);
    % compute preferences matrix
    f = -c_mat./A_mat;
    
    realisable = 1;
    beta = b;
    I = 1:m;
    J = 1:n;
    x = zeros(n, 1);
    zeta = 0;
    while ~isempty(J) && realisable == 1
        regmax = -Inf;
        for j = J
            Rj = I(A_mat(:,j) < beta);
            if isempty(Rj)
                realisable = 0;
            else
                [~, k] = max(f(Rj,j));
                if isempty(Rj(Rj ~= k))
                    reg = Inf;
                else
                    reg = f(k,j) - max(f(Rj(Rj ~= k),j));
                end
                if reg > regmax
                    regmax = reg;
                    imax = k;
                    jmax = j;
                end
            end
        end
        if realisable == 1
            x(jmax) = imax;
            zeta = zeta - c_mat(imax, jmax);
            beta(imax) = beta(imax) - A_mat(imax, jmax);
            J = J(J ~= jmax);
        end
    end
    
    sol = zeros(m,n);
    if realisable == 1
        for j = 1:n
            k = x(j);
            I_bis = I(I~=k);
            K = -c_mat(I_bis,:);
            K = K(A_mat(I_bis,j) <= beta(I_bis));
            if ~isempty(K)
                [maxval, q] = max(K);
                if maxval > -c_mat(k,j)
                    x(j) = q;
                    zeta = zeta + c_mat(k,j) - c_mat(q,j);
                    beta(k) = beta(k) + A_mat(k,j);
                    beta(q) = beta(q) - A_mat(q,j);
                end
            end
        end
    
        % create matrix solution
        for j = 1:n
            sol(x(j),j) = 1;
        end
    end
end

function M = vector_to_matrix(m, n, v)
    M = zeros(m,n);
    for i = 1:m
        M(i,:) = v((i-1)*n+1:(i-1)*n+n);
    end
end
