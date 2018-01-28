function [m, n, A_vec, A, b, c, Aeq, beq, x0, r_values] = read_problem(file)
    % load data
    inf = fopen(file, 'r');
    r_dimensions = textscan(inf, '%f', 2);
    r_values = textscan(inf, '%f');
    fclose(inf);

    m = r_dimensions{1}(1);
    n = r_dimensions{1}(2);
    nvars = m*n;

    tmp = r_values{1};
    c = tmp(1:nvars);
    A_vec = tmp(nvars+1:2*nvars);
    b = tmp(2*nvars+1:end);
    
    % compose A matrix from A vector (for constraints (3))
    A = zeros(m, nvars);
    for i = 1:m
        for j = 1:n
            A(i, (i-1)*n+j) = A_vec((i-1)*n+j);
        end
    end

    % compose Aeq matrix (for constraints (4) on column sum)
    Aeq = repmat(eye(n), 1, m);
    beq = ones(n, 1);

    x0 = sol_initial(m, n, c, A_vec, b, A, Aeq, beq);
end

function x0 = sol_initial(m, n, c, A_vec, b, A, Aeq, beq)
    [x0_mat, x0_realisable] = heuristique_regret(m, n, c, A_vec, b);
    x0 = matrix_to_vector(x0_mat);

    % if x0 is not feasible:
    % find initial feasible solution with interrupted intlinprog
%    if x0_realisable == 0
    if sol_check(x0_mat, A_vec, b) == 0
        opts = optimoptions('intlinprog', 'MaxFeasiblePoints', 1);
        nvars = m*n;
        lb = zeros(nvars,1);
        ub = ones(nvars,1);
        [x0, ~] = intlinprog(c, 1:nvars, A, b, Aeq, beq, lb, ub, [], opts);
    end
end

function v = matrix_to_vector(M)
     v = reshape(M', [], 1);
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
            disp(['Sum of column j = ',num2str(sum(x(:,j))),': ',num2str(sum(x(:,j)))])
            ok = 0;
            return;
        elseif abs(sum(x(:,j)) - 1) > 0
            x(:,j) = round(x(:,j));
        end
    end
end
