function [sol, realisable, zeta] = heuristique_regret(m, n, c, A, b)
    % get A and c as matrix
    A_mat = vector_to_matrix(m, n, A);
    c_mat = vector_to_matrix(m, n, c);
    % compute preferences matrix
    %f = -c_mat;
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
            Rj = I(A_mat(:,j) <= beta);
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
%        disp(['reg=',num2str(reg),', imax=',num2str(imax),', jmax=',num2str(jmax)])
        if realisable == 1
            x(jmax) = imax;
            zeta = zeta - c_mat(imax, jmax);
            beta(imax) = beta(imax) - A_mat(imax, jmax);
            J = J(J ~= jmax);
        end
    end
    
    sol = zeros(m,n);
%     for j = 1:n
%         sol(x(j),j) = 1;
%     end
%     realisable
%     x0_realisable = sol_check(sol, A, b)
    if realisable == 1
        for j = 1:n
            % get index of the best agent for task j found yet
            k = x(j);
            % compute indices i matching the two conditions:
            % - i != k
            % - agent i can perform the task j: a_ij <= beta_i
            I_bis = [];
            last_i = 1;
            for i = 1:m
                if i ~= k && A_mat(i,j) <= beta(i)
                    I_bis(last_i) = i;
                    last_i = last_i + 1;
                end
            end
            % get corresponding utilities
            K = -c_mat(I_bis,j);
            % if at least one other agent can perform the task j,
            if ~isempty(K)
                % get the best utility from all the possible replacements
                [minval, qtmp] = min(K);
                % if it is better than the current one, reassign task j
                if minval < -c_mat(k,j)
                    q = I_bis(qtmp);
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
    if m == 1
        if isrow(v)
            M = v;
            return;
        else
            M = v';
            return;
        end
    elseif n == 1
        if isrow(v)
            M = v';
            return;
        else
            M = v;
            return;
        end
    end

    M = zeros(m,n);
    for i = 1:m
        M(i,:) = v((i-1)*n+1:(i-1)*n+n);
    end
end

function ok = sol_check(x, A_vec, b)
    disp('Checking solution for feasibility.');
    ok = 1;
    [m, n] = size(x);
    A = vector_to_matrix(m, n, A_vec);
    tol = 0.1;
    % check weight constraints (3)
    for i = 1:m
        if sum(A(i,:).*x(i,:)) > b(i)
            disp(['Solution is infeasible: weight constraints (3) not ok for agent i=',num2str(i),':']);
            disp(['A*x = ',num2str(sum(A(i,:).*x(i,:))),' > b = ',num2str(b(i))])
            ok = 0;
            return;
        end
    end
    % check equality constraints (4)
    for j = 1:n
        if abs(sum(x(:,j)) - 1) > tol
            disp(['Solution is infeasible: equality constraints (4) not ok for task j=',num2str(j),'.']);
            disp(['Sum of column j = ',num2str(sum(x(:,j)))])
            ok = 0;
            return;
        elseif abs(sum(x(:,j)) - 1) > 0
            x(:,j) = round(x(:,j));
        end
    end
end
