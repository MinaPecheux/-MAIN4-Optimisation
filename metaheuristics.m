function metaheuristics()
    % load data
    inf = fopen('data-1/1/a10100', 'r');
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

%    x0 = zeros(nvars, 1);
    [x0_mat, x0_realisable] = heuristique_regret(m, n, c, A_vec, b);
    x0 = matrix_to_vector(x0_mat);

    % solve (P) problem
    % -----------------
    % if initial solution is feasible:
    if x0_realisable == 1
        run_problem(m, n, x0, c, A, b, Aeq, beq, A_vec);
    else
        % else find initial feasible solution with interrupted intlinprog
        opts = optimoptions('intlinprog', 'MaxTime', 10);
        [x0, ~] = intlinprog(c, 1:nvars, A, b, Aeq, beq, lb, ub, opts);
        run_problem(m, n, x0, c, A, b, Aeq, beq, A_vec);
    end

    %%
    % solve (L) problem without constraints (10)
    % ------------------------------------------
    % load data
    inf = fopen('data-1/1/a0505', 'r');
    r_dimensions = textscan(inf, '%f', 2);
    r_values = textscan(inf, '%f');
    fclose(inf);

    m = r_dimensions{1}(1);
    n = r_dimensions{1}(2);

    tmp = r_values{1};
    c = [tmp(1:m*n);zeros(m*n*m,1)];    % add zeros for no-cost z variables
    A_vec = tmp(m*n+1:2*m*n);
    b = tmp(2*m*n+1:end);

    % objective function
    f = @(x) obj_func(m, n, x, c);
    
    nvars = m*n * (1+m); % (x_11,...,x_mn,z_11^1,...,z_1n^1,...,z_n1^m,z_nm^m)

    % compose A matrix (for constraints (11) + (13) + (14) + (15))
    A = zeros(m*m + 3*m*m*n, nvars);
                            % SIZE OF THIS MATRIX:
                            % (m*m rows because there are constraints for i,k in {1,...m}
                            % + m*m*n for constraints (13) + m*m*n for constraints(14))
                            % + m*m*n for constraints (15))
    row = 1;
    for i = 1:m
        for k = 1:m
            disp([num2str(i), ', ', num2str(k), ' // ', num2str(row), ' // ', num2str((i-1)*m+k)]);
            A(row, row) = -(b(i) - A_vec((i-1)*m+k));
            for j = 1:n
                A(row, m*m + (i-1)*m*m + k + (j-1)*m) = A_vec((i-1)*m + j);
            end
            row = row + 1;
        end
    end
    % complete A matrix (for constraints (13) + (14))
    row = 1;
    for i = 1:m
        for j = 1:n
            for k = 1:m
                % constraints (13)
                A(m*m + row, j + (i-1)*m) = -1;
                A(m*m + row, m*m + row) = 1;
                % constraints (14)
                A(m*m + m*m*n + row, k + (i-1)*m) = -1;
                A(m*m + m*m*n + row, m*m + row) = 1;
                row = row + 1;
            end
        end
    end
    % complete A matrix (for constraints (15))
    row = 1;
    for i = 1:m
        for j = 1:n
            for k = 1:m
                % constraints (15)
                if j == k
                    A(m*m + 2*m*m*n + row, j + (i-1)*m) = 2;
                else
                    A(m*m + 2*m*m*n + row, j + (i-1)*m) = 1;
                    A(m*m + 2*m*m*n + row, k + (i-1)*m) = 1;
                end
                A(m*m + 2*m*m*n + row, m*m + row) = -1;
                row = row + 1;
            end
        end
    end
    % remove overlapping (doubled) constraints
    A = unique(A, 'stable', 'rows');

    % compose B matrix (for constraints (10) + (13) + (14) + (15))
    B = [zeros(size(A,1) - m*m*n,1);ones(m*m*n,1)];

    % compose Aeq matrix (for constraint (4) on column sum)
    Aeq = zeros(n, nvars);
    Aeq(:,1:m*n) = repmat(eye(n), 1, m);
    beq = ones(n, 1);

    x0 = zeros(nvars, 1);

    lb = zeros(nvars, 1);
    ub = [ones(m*n, 1);Inf*ones(nvars - m*n, 1)];

    % solve problem
    opts = optimoptions('patternsearch', 'PollMethod', 'GSSPositiveBasisNp1', 'TolBind',2);
    tic
    %[x_ps, obj_ps] = patternsearch(f, x0, A, B, Aeq, beq, lb, ub, [], opts);
    [x_ps, obj_ps] = ga(f, nvars, A, B, Aeq, beq, lb, ub);
    toc

    % display results
    disp_x_ps = sol_display(m, n, x_ps(1:m*n));
    disp_z_ps = sol_display(m*m, n, x_ps(m*n+1:end));

    disp('PATTERN SEARCH');
    disp('Solution x =')
    disp(disp_x_ps);
    disp('Solution z =')
    disp(disp_z_ps);
    disp(['Objective value (min): ', num2str(obj_ps)]);

    % re-solve (P) problem
    % --------------------
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

    run_problem(m, n, x_ps(1:m*n), c, A, b, Aeq, beq);

    %%
    % solve (L) problem without constraints (11)
    % ------------------------------------------
    % load data
    inf = fopen('data-1/1/a05100', 'r');
    r_dimensions = textscan(inf, '%f', 2);
    r_values = textscan(inf, '%f');
    fclose(inf);

    m = r_dimensions{1}(1);
    n = r_dimensions{1}(2);

    tmp = r_values{1};
    c = [tmp(1:m*n);zeros(m*n*m,1)];    % add zeros for no-cost z variables
    A_vec = tmp(m*n+1:2*m*n);
    b = tmp(2*m*n+1:end);

    % objective function
    f = @(x) obj_func(m, n, x, c);

    nvars = m*n * (1+m); % (x_11,...,x_mn,z_11^1,...,z_1n^1,...,z_n1^m,z_nm^m)

    % compose A matrix (for constraints (10) + (13) + (14) + (15))
    A = zeros(m*m + 3*m*m*n, nvars);
                            % SIZE OF THIS MATRIX:
                            % (m*m rows because there are constraints for i,k in {1,...m}
                            % + m*m*n for constraints (13) + m*m*n for constraints(14))
                            % + m*m*n for constraints (15))
    row = 1;
    for i = 1:m
        for k = 1:m
            A(row, row) = b(i);
            for j = 1:n
                if j ~= k
                    A(row, (i-1)*m + j) = A_vec((i-1)*m + j);
                    A(row, m*m + (i-1)*m + (j-1)*m + k) = -A_vec((i-1)*m + j);
                end
            end
            row = row + 1;
        end
    end
    % complete A matrix (for constraints (13) + (14))
    row = 1;
    for i = 1:m
        for j = 1:n
            for k = 1:m
                % constraints (13)
                A(m*m + row, j + (i-1)*m) = -1;
                A(m*m + row, m*m + row) = 1;
                % constraints (14)
                A(m*m + m*m*n + row, k + (i-1)*m) = -1;
                A(m*m + m*m*n + row, m*m + row) = 1;
                row = row + 1;
            end
        end
    end
    % complete A matrix (for constraints (15))
    row = 1;
    for i = 1:m
        for j = 1:n
            for k = 1:m
                % constraints (15)
                if j == k
                    A(m*m + 2*m*m*n + row, j + (i-1)*m) = 2;
                else
                    A(m*m + 2*m*m*n + row, j + (i-1)*m) = 1;
                    A(m*m + 2*m*m*n + row, k + (i-1)*m) = 1;
                end
                A(m*m + 2*m*m*n + row, m*m + row) = -1;
                row = row + 1;
            end
        end
    end
    % remove overlapping (doubled) constraints
    A = unique(A, 'stable', 'rows');

    % compose B matrix (for constraints (10) + (13) + (14) + (15))
    B = [zeros(size(A,1) - m*m*n,1);ones(m*m*n,1)];
    row = 1;
    for i = 1:m
        for k = 1:m
            B(row) = b(i);
            row = row + 1;
        end
    end

    % compose Aeq matrix (for constraint (4) on column sum)
    Aeq = zeros(n, nvars);
    Aeq(:,1:m*n) = repmat(eye(n), 1, m);
    beq = ones(n, 1);

    x0 = zeros(nvars, 1);

    lb = zeros(nvars, 1);
    ub = [ones(m*n, 1);Inf*ones(nvars - m*n, 1)];

    % solve problem
    opts = optimoptions('patternsearch','TolBind',0.2);
    tic
    [x_ps, obj_ps] = patternsearch(f, x0, A, B, Aeq, beq, lb, ub, [], opts);
    toc

    % display results
    disp_x_ps = sol_display(m, n, x_ps(1:m*n));
    disp_z_ps = sol_display(m*m, n, x_ps(m*n+1:end));

    disp('PATTERN SEARCH');
    disp('Solution x =')
    disp(disp_x_ps);
    disp('Solution z =')
    disp(disp_z_ps);
    disp(['Objective value (min): ', num2str(obj_ps)]);

    % re-solve (P) problem
    % --------------------
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

    run_problem(m, n, x_ps(1:m*n), c, A, b, Aeq, beq);
end

function run_problem(m, n, x0, c, A, b, Aeq, beq, A_vec)
    % prepare vars
    nvars = length(x0);
    lb = zeros(nvars, 1);
    ub = ones(nvars, 1);
    f = @(x) obj_func(m, n, x, c);
    
    % solve (P) problem with the 3 heuristics:
    % ----------------------------------------
    % PATTERNSEARCH
    % {'GPSPositiveBasis2N'} | 'GPSPositiveBasisNp1' | 'GSSPositiveBasis2N' | 'GSSPositiveBasisNp1' | 'MADSPositiveBasis2N' | 'MADSPositiveBasisNp1'
    opts = optimoptions('patternsearch', 'ScaleMesh', false, 'MeshTolerance', 0.99); % force integer solutions
    tic
    [x_ps, obj_ps] = patternsearch(f, x0, A, b, Aeq, beq, lb, ub, [], opts);
    t_ps = toc;
    
    % GA: ga with integer vars takes no equality constraints: workaround is to add
    % two inequality constraints equivalent to the equality constraint:
    A2 = [A;Aeq;-Aeq];
    b2 = [b;beq;-beq];
    %opts = optimoptions('ga','MaxStallGenerations',50,'FunctionTolerance',1e-10,'MaxGenerations',300,'ConstraintTolerance',1e-10);
    opts = optimoptions('ga', 'InitialPopulation', x0');
    tic
    [x_ga, ~] = ga(f, nvars, A2, b2, [], [], lb, ub, [], 1:nvars, opts);
    %x_ga = 0; obj_ga = 0;
    t_ga = toc;
    
    % INTLINPROG
    tic
    [x_intlinprog, obj_intlinprog] = intlinprog(c, 1:nvars, A, b, Aeq, beq, lb, ub, x0);
    t_ilp = toc;

    % clean and display results
    disp_x_ps = sol_cleanup_ps(m, n, x_ps, A_vec, c);
    obj_ps = obj_func(m, n, matrix_to_vector(disp_x_ps), c);
    %disp_x_ps = vector_to_matrix(m, n, x_ps);
    disp_x_ga = sol_cleanup_ga(m, n, x_ga, A_vec, c);
    obj_ga = obj_func(m, n, matrix_to_vector(disp_x_ga), c);
    %disp_x_ga = sol_display(m, n, x_ga);
    %disp_x_ga = x_ga;
    disp_x_int = vector_to_matrix(m, n, x_intlinprog);

    disp(['PATTERN SEARCH - Objective value (min): ', num2str(obj_ps)]);
    disp(['Time: ', num2str(t_ps), ' sec']);
    disp(['Feasible: ', num2str(sol_check(disp_x_ps, A_vec, b))]);
    disp('Solution x =');
    disp(disp_x_ps);

    disp(['GENETIC ALGORITHM - Objective value (min): ', num2str(obj_ga)]);
    disp(['Time: ', num2str(t_ga), ' sec']);
    disp(['Feasible: ', num2str(sol_check(disp_x_ga, A_vec, b))]);
    disp('Solution x =')
    disp(disp_x_ga);

    disp(['INTLINPROG - Objective value (min): ', num2str(obj_intlinprog)]);
    disp(['Time: ', num2str(t_ilp), ' sec']);
    disp(['Feasible: ', num2str(sol_check(disp_x_int, A_vec, b))]);
    disp('Solution x =');
    disp(disp_x_int);
end

function obj = obj_func(m, n, x, c)
    obj = 0;
    for i = 1:m
        for j = 1:n
            obj = obj + c((i-1)*n+j)*x((i-1)*n+j);
        end
    end
end

function v = matrix_to_vector(M)
     v = reshape(M', [], 1);
end

function M = vector_to_matrix(m, n, v)
    M = zeros(m,n);
    for i = 1:m
        M(i,:) = v((i-1)*n+1:(i-1)*n+n);
    end
end

function x = sol_cleanup_ps(m, n, dirty_x, A, c)
    % get A and c as matrix
    A_mat = vector_to_matrix(m, n, A);
    c_mat = vector_to_matrix(m, n, c);
    % compute ratio matrix
    ratios = -c_mat./A_mat;
    
    % new x is only zero, except for 1 coeff at the max efficacity pos
    x = zeros(m, n);
    % get dirty_x as matrix
    dirty_x_mat = vector_to_matrix(m, n, dirty_x);
    % foreach column of the matrix
    for j = 1:n
        col = dirty_x_mat(:,j);
        % check if column has non-integer values
        if floor(col) ~= col
            % get matching ratio column
            r = ratios(:, j);
            % apply 1 coeff to x vector at the pos corresponding to r max
            % value
            x(r == max(r), j) = 1;
        else
            x(:,j) = dirty_x_mat(:,j);
        end
    end
end

function x = sol_cleanup_ga(m, n, dirty_x, A, c)
    % get A and c as matrix
    A_mat = vector_to_matrix(m, n, A);
    c_mat = vector_to_matrix(m, n, c);
    % compute ratio matrix
    ratios = -c_mat./A_mat;
    
    % new x is only zero, except for 1 coeff at the max(r) pos
    x = zeros(m, n);
    % get dirty_x as matrix
    dirty_x_mat = vector_to_matrix(m, n, dirty_x);
    % foreach column of the matrix
    for j = 1:n
        col = dirty_x_mat(:,j);
        % get all non-zero values
        col_nonzeros_ind = find(col);
        % if there are more than one, clean
        if length(col_nonzeros_ind) ~= 1
            % get matching ratio columns at non-zero values indexes
            r = ratios(col_nonzeros_ind, j);
            % apply 1 coeff to x vector at the pos corresponding to r max
            % value
            x(r == max(r), j) = 1;
        else
            x(:,j) = dirty_x_mat(:,j);
        end
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
