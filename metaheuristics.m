function metaheuristics()
    % load data
    inf = fopen('data-1/1/a20200', 'r');
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
    Aeq = zeros(n, nvars);
    for i = 0:m-1
        Aeq(:,i*n+1:i*n+n) = eye(n);
    end
    beq = ones(n, 1);

    x0 = zeros(nvars, 1);

    % solve (P) problem
    % -----------------
    run_problem(m, n, x0, c, A, b, Aeq, beq);

    %%
    % solve (L) problem without constraints (10)
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

    % compose A matrix (for constraints (11) + (13) + (14) + (15))
    A = zeros(m*m + 3*m*m*n, nvars);
                            % SIZE OF THIS MATRIX:
                            % (m*m rows because there are constraints for i,k in {1,...m}
                            % + m*m*n for constraints (13) + m*m*n for constraints(14))
                            % + m*m*n for constraints (15))
    row = 1;
    for i = 1:m
        for k = 1:m
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
    Aeq = zeros(m, nvars);
    for i = 1:m
        Aeq(:,(i-1)*m+1:(i-1)*m+m) = eye(m);
    end
    beq = ones(m, 1);

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
    Aeq = zeros(n, nvars);
    for i = 0:m-1
        Aeq(:,i*n+1:i*n+n) = eye(n);
    end
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
    Aeq = zeros(m, nvars);
    for i = 1:m
        Aeq(:,(i-1)*m+1:(i-1)*m+m) = eye(m);
    end
    beq = ones(m, 1);

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
    Aeq = zeros(n, nvars);
    for i = 0:m-1
        Aeq(:,i*n+1:i*n+n) = eye(n);
    end
    beq = ones(n, 1);

    run_problem(m, n, x_ps(1:m*n), c, A, b, Aeq, beq);
end

function run_problem(m, n, x0, c, A, b, Aeq, beq)
    % prepare vars
    nvars = length(x0);
    lb = zeros(nvars, 1);
    ub = ones(nvars, 1);
    f = @(x) obj_func(m, n, x, c);

    % solve (P) problem with the 3 heuristics:
    % ----------------------------------------
    % PATTERNSEARCH
    opts = optimoptions('patternsearch', 'ScaleMesh', false, 'MeshTolerance', 0.99); % force integer solutions
    tic
    [x_ps, obj_ps] = patternsearch(f, x0, A, b, Aeq, beq, lb, ub, [], opts);
    toc
    
    % GA: ga with integer vars takes no equality constraints: workaround is to add
    % two inequality constraints equivalent to the equality constraint:
    A2 = [A;Aeq;-Aeq];
    b2 = [b;beq;-beq];
    opts = optimoptions('ga','MaxStallGenerations',50,'FunctionTolerance',1e-10,'MaxGenerations',300,'ConstraintTolerance',1e-10);
    tic
    [x_ga, obj_ga] = ga(f, nvars, A2, b2, [], [], lb, ub, [], 1:nvars);
    %x_ga = 0; obj_ga = 0;
    toc
    
    % INTLINPROG
    tic
    [x_intlinprog, obj_intlinprog] = intlinprog(c, 1:nvars, A, b, Aeq, beq, lb, ub);
    toc

    % display results
    disp_x_ps = sol_display(m, n, x_ps);
    disp_x_ga = sol_display(m, n, x_ga);
    %disp_x_ga = x_ga;
    disp_x_int = sol_display(m, n, x_intlinprog);

    disp('PATTERN SEARCH');
    disp('Solution x =');
    disp(disp_x_ps);
    disp(['Objective value (min): ', num2str(obj_ps)]);

    disp('GENETIC ALGORITHM');
    disp('Solution x =')
    disp(disp_x_ga);
    disp(['Objective value (min): ', num2str(obj_ga)]);

    disp('INTLINPROG');
    disp('Solution x =');
    disp(disp_x_int);
    disp(['Objective value (min): ', num2str(obj_intlinprog)]);
end

function obj = obj_func(m, n, x, c)
    obj = 0;
    for i = 1:m
        for j = 1:n
            obj = obj + c((i-1)*n+j)*x((i-1)*n+j);
        end
    end
end

function d = sol_display(m, n, x)
    d = zeros(m,n);
    for i = 1:m
        d(i,:) = x((i-1)*n+1:(i-1)*n+n);
    end
end
