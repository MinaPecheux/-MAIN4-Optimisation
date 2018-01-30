% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% ---------------------------------------------------
% Script comparant les métaheuristiques patternsearch
% et ga selon la solution initiale x0 choisie :
% - solution x0 donnée
% - solution de la relaxation continue de (L) sans
% . les contraintes (10)
% - solution de la relaxation continue de (L) sans
% . les contraintes (11)
%
% Chaque cas est dans une section différente qui peut
% être évaluée de manière autonome.
%
% ATTENTION ! Les résolutions du problème (L) sont
% relativement longues (~1h30) même pour la plus
% petite instance a05100.
% ---------------------------------------------------

function metaheuristics()
    % solve (P) problem directly
    % --------------------------
    
    % load file data + get initial feasible solution
    [m, n, A_vec, A, b, c, Aeq, beq, x0] = read_problem('../dat/1/a20100');
    % solve problem directly with patternsearch, ga, intlinprog
    %   -> option 'true' at the end to display solutions (turn to
    %   'false' to only get objective value, feasibility and exec time)
    run_problem(m, n, x0, c, A, b, Aeq, beq, A_vec, false);

    %%
    % solve (L) problem without constraints (10)
    % ------------------------------------------
    % load file data + get initial feasible solution
    [m, n, A_vec, ~, b, c, ~, ~, x0, r_values] = read_problem('../dat/1/a05100');

    % useful matrix
    a = vector_to_matrix(m, n, A_vec);
    % add zeros for no-cost z variables
    c = [c;zeros(m*n*m,1)];
    x0 = [x0;zeros(m*n*m,1)];
    % objective function
    f = @(x) obj_func(m, n, x, c);
    
    nvars = m*n * (1+m); % (x_11,...,x_mn,z_11^1,...,z_1n^1,...,z_n1^m,z_nm^m)

    % compose A matrix (for constraints (11) + (13) + (14) + (15))
    A = zeros(m*m + 3*m*m*n, nvars);
                            % SIZE OF THIS MATRIX:
                            % (m*m rows because there are constraints for i,k in {1,...m}
                            % + m*m*n for constraints (13) + m*m*n for constraints(14))
                            % + m*m*n for constraints (15))
    % start A matrix (constraints (11))
    row = 1;
    for i = 1:m
        for k = 1:m
            A(row, (i-1)*n + k) = b(i);
            for j = 1:n
                if j ~= k
                    A(row, (i-1)*n + j) = a(i,j);
                    A(row, m*n + (i-1)*m*n + k + (j-1)*m) = -a(i,j);
                end
            end
            row = row + 1;
        end
    end
    % complete A matrix (constraints (13))
    A(row:row + m*m*n - 1, m*n+1:end) = eye(m*m*n);
    pattern = zeros(m*n, n);
    for j = 1:n
        for i = 1:m
            t = (j-1)*m;
            pattern(t+1:t+m, j) = -1;
        end
    end
    for i = 1:m
        A(row:row + m*n - 1, (i-1)*n+1:(i-1)*n+n) = pattern;
        row = row + m*n;
    end
    % complete A matrix (constraints (14))
    A(row:row + m*m*n - 1, m*n+1:end) = eye(m*m*n);
    pattern = zeros(m*n, n);
    for j = 1:n
        for i = 1:m
            t = (j-1)*m;
            pattern(t+1:t+m, 1:m) = -eye(m);
        end
    end
    for i = 1:m
        A(row:row + m*n - 1, (i-1)*n+1:(i-1)*n+n) = pattern;
        row = row + m*n;
    end
    % complete A matrix (constraints (15))
    A(row:row + m*m*n - 1, m*n+1:end) = -eye(m*m*n);
    pattern = zeros(m*n, n);
    r = 1;
    for j = 1:n
        for k = 1:m
            if j == k
                pattern(r,k) = 2;
            else
                pattern(r,j) = 1;
                pattern(r,k) = 1;
            end
            r = r + 1;
        end
    end
    for i = 1:m
        A(row:row + m*n - 1, (i-1)*n+1:(i-1)*n+n) = pattern;
        row = row + m*n;
    end
    % remove overlapping (doubled) constraints
    A = unique(A, 'stable', 'rows');

    % compose B matrix (for constraints (11) + (13) + (14) + (15))
    B = [zeros(size(A,1) - m*m*n,1);ones(m*m*n,1)];

    % compose Aeq matrix (for constraint (4) on column sum)
    Aeq = zeros(n, nvars);
    Aeq(:,1:m*n) = repmat(eye(n), 1, m);
    beq = ones(n, 1);

    lb = zeros(nvars, 1);
    ub = [ones(m*n, 1);Inf*ones(nvars - m*n, 1)];
    
    % solve (L) problem with patternsearch
    opts = optimoptions('patternsearch', 'PollMethod', 'GSSPositiveBasisNp1', 'UseCompletePoll', true, 'Display', 'iter');
    disp('PATTERNSEARCH: Starting (L) resolution...');
    tic
    [x_ps, obj_ps] = patternsearch(f, x0, A, B, Aeq, beq, lb, ub, [], opts);
    tps_l = toc;

    % display results
    disp_x_ps = vector_to_matrix(m, n, x_ps(1:m*n));
    disp_z_ps = vector_to_matrix(m*m, n, x_ps(m*n+1:end));

    disp('PATTERNSEARCH');
    disp(['Solved (L): ',num2str(tps_l),' sec']);
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

    x0 = sol_cleanup_ps(m, n, x_ps(1:m*n), A_vec, b, c);
    disp(['x0 feasible? Feasibility: ',num2str(sol_check(x0, A_vec, b))])
    run_problem(m, n, matrix_to_vector(x0), c, A, b, Aeq, beq, A_vec, false);

    %%
    % solve (L) problem without constraints (11)
    % ------------------------------------------
    % load file data + get initial feasible solution
    [m, n, A_vec, ~, b, c, ~, ~, x0, r_values] = read_problem('../dat/1/a05100');

    % useful matrix
    a = vector_to_matrix(m, n, A_vec);
    % add zeros for no-cost z variables
    c = [c;zeros(m*n*m,1)];
    x0 = [x0;zeros(m*n*m,1)];
    % objective function
    f = @(x) obj_func(m, n, x, c);
    
    nvars = m*n * (1+m); % (x_11,...,x_mn,z_11^1,...,z_1n^1,...,z_n1^m,z_nm^m)

    % compose A matrix (for constraints (10) + (13) + (14) + (15))
    A = zeros(m*m + 3*m*m*n, nvars);
                            % SIZE OF THIS MATRIX:
                            % (m*m rows because there are constraints for i,k in {1,...m}
                            % + m*m*n for constraints (13) + m*m*n for constraints(14))
                            % + m*m*n for constraints (15))
    % start A matrix (constraints (10))
    row = 1;
    for i = 1:m
        for k = 1:m
            A(row, (i-1)*n + k) = -(b(i) - a(i,k));
            for j = 1:n
                A(row, m*n + (i-1)*m*n + k + (j-1)*m) = a(i,j);
            end
            row = row + 1;
        end
    end
    % complete A matrix (constraints (13))
    A(row:row + m*m*n - 1, m*n+1:end) = eye(m*m*n);
    pattern = zeros(m*n, n);
    for j = 1:n
        for i = 1:m
            t = (j-1)*m;
            pattern(t+1:t+m, j) = -1;
        end
    end
    for i = 1:m
        A(row:row + m*n - 1, (i-1)*n+1:(i-1)*n+n) = pattern;
        row = row + m*n;
    end
    % complete A matrix (constraints (14))
    A(row:row + m*m*n - 1, m*n+1:end) = eye(m*m*n);
    pattern = zeros(m*n, n);
    for j = 1:n
        for i = 1:m
            t = (j-1)*m;
            pattern(t+1:t+m, 1:m) = -eye(m);
        end
    end
    for i = 1:m
        A(row:row + m*n - 1, (i-1)*n+1:(i-1)*n+n) = pattern;
        row = row + m*n;
    end
    % complete A matrix (constraints (15))
    A(row:row + m*m*n - 1, m*n+1:end) = -eye(m*m*n);
    pattern = zeros(m*n, n);
    r = 1;
    for j = 1:n
        for k = 1:m
            if j == k
                pattern(r,k) = 2;
            else
                pattern(r,j) = 1;
                pattern(r,k) = 1;
            end
            r = r + 1;
        end
    end
    for i = 1:m
        A(row:row + m*n - 1, (i-1)*n+1:(i-1)*n+n) = pattern;
        row = row + m*n;
    end
    % remove overlapping (doubled) constraints
    A = unique(A, 'stable', 'rows');

    % compose B matrix (for constraints (10) + (13) + (14) + (15))
    B = [zeros(size(A,1) - m*m*n,1);ones(m*m*n,1)];
    B(1:m*m) = repelem(b,m,1);

    % compose Aeq matrix (for constraint (4) on column sum)
    Aeq = zeros(n, nvars);
    Aeq(:,1:m*n) = repmat(eye(n), 1, m);
    beq = ones(n, 1);

    lb = zeros(nvars, 1);
    ub = [ones(m*n, 1);Inf*ones(nvars - m*n, 1)];

    % solve (L) problem with patternsearch
    opts = optimoptions('patternsearch', 'PollMethod', 'GSSPositiveBasisNp1', 'UseCompletePoll', true, 'Display', 'iter');
    disp('PATTERNSEARCH: Starting (L) resolution...');
    tic
    [x_ps, obj_ps] = patternsearch(f, x0, A, B, Aeq, beq, lb, ub, [], opts);
    tps_l = toc;

    % display results
    disp_x_ps = vector_to_matrix(m, n, x_ps(1:m*n));
    disp_z_ps = vector_to_matrix(m*m, n, x_ps(m*n+1:end));

    disp('PATTERNSEARCH');
    disp(['Solved (L): ',num2str(tps_l),' sec']);
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

    x0 = sol_cleanup_ps(m, n, x_ps(1:m*n), A_vec, b, c);
    disp(['x0 feasible? Feasibility: ',num2str(sol_check(x0, A_vec, b))])
    run_problem(m, n, matrix_to_vector(x0), c, A, b, Aeq, beq, A_vec, false);
end

% -----------------------------------------------------------
% Fonction résolvant le problème (P) à partir des
% données lues dans le fichier et d'une solution
% initiale x0 choisie.
% Input(s):
% - m, n [entiers] : dimensions du problème (resp. nombre
% . d'agents et nombre de tâches)
% - x0 [vecteur de réels] : solution initiale réalisable
% . pour (P)
% - c [vecteur de réels] : vecteur des coûts
% - A [matrice de réels] : sous-matrice de productivité
% . correspondant aux contraintes d'inégalité
% - b [vecteur de réels] : sous-vecteur des disponibilités
% . correspondant aux contraintes d'inégalité
% - Aeq [matrice de réels] : sous-matrice de productivité
% . correspondant aux contraintes d'égalité
% - beq [vecteur de réels] : sous-vecteur des disponibilités
% . correspondant aux contraintes d'égalité
% - A_vec [vecteur de réels] : vecteur de productivité
% - show_solutions [bool] : option pour montrer les solutions
% . trouvées ou afficher seulement la valeur objectif, le
% . temps d'exécution et la réalisabilité de la solution
% -----------------------------------------------------------
function run_problem(m, n, x0, c, A, b, Aeq, beq, A_vec, show_solutions)
    % prepare vars
    nvars = length(x0);
    lb = zeros(nvars, 1);
    ub = ones(nvars, 1);
    f = @(x) obj_func(m, n, x, c);
    
    % solve (P) problem with the 3 heuristics:
    % ----------------------------------------
    % PATTERNSEARCH
    opts = optimoptions('patternsearch', 'PollMethod', 'GSSPositiveBasisNp1', 'UseCompletePoll', true, 'ScaleMesh', false, 'MeshTolerance', 0.99); % force integer solutions
    tic
    [x_ps, ~] = patternsearch(f, x0, A, b, Aeq, beq, lb, ub, [], opts);
    t_ps = toc;
    
    % GA: ga with integer vars takes no equality constraints: workaround is to add
    % two inequality constraints equivalent to the equality constraint:
    A2 = [A;Aeq;-Aeq];
    b2 = [b;beq;-beq];
    opts = optimoptions('ga', 'InitialPopulation', x0');
    tic
    [x_ga, ~] = ga(f, nvars, A2, b2, [], [], lb, ub, [], 1:nvars, opts);
    t_ga = toc;
    
    % INTLINPROG
    tic
    [x_intlinprog, obj_intlinprog] = intlinprog(c, 1:nvars, A, b, Aeq, beq, lb, ub, x0);
    t_ilp = toc;

    % clean and display results
    % -------------------------
    disp_x_ps = sol_cleanup_ps(m, n, x_ps, A_vec, b, c);
    obj_ps = obj_func(m, n, matrix_to_vector(disp_x_ps), c);
    disp_x_ga = sol_cleanup_ga(m, n, x_ga, A_vec, b, c);
    obj_ga = obj_func(m, n, matrix_to_vector(disp_x_ga), c);
    disp_x_int = vector_to_matrix(m, n, x_intlinprog);

    disp(['PATTERN SEARCH - Objective value (min): ', num2str(obj_ps)]);
    disp(['Time: ', num2str(t_ps), ' sec']);
    disp(['Feasible (after cleanup): ', num2str(sol_check(disp_x_ps, A_vec, b))]);
    if obj_ps ~= 0
        disp(['Solution quality: ', num2str(100*(obj_ps - obj_intlinprog)/obj_ps)]);
        disp(['> Solver efficiency (100/[quality * exec time]): ', num2str(100/(t_ps * 100*(obj_ps - obj_intlinprog)/obj_ps))]);
    end
    if show_solutions
        disp('Solution x =');
        disp(disp_x_ps);
    end

    disp([newline, 'GENETIC ALGORITHM - Objective value (min): ', num2str(obj_ga)]);
    disp(['Time: ', num2str(t_ga), ' sec']);
    disp(['Feasible (after cleanup): ', num2str(sol_check(disp_x_ga, A_vec, b))]);
    if obj_ga ~= 0
        disp(['Solution quality: ', num2str(100*(obj_ga - obj_intlinprog)/obj_ga)]);
        disp(['> Solver efficiency (100/[quality * exec time]): ', num2str(100/(t_ga * 100*(obj_ga - obj_intlinprog)/obj_ga))]);
    end
    if show_solutions
        disp('Solution x =')
        disp(disp_x_ga);
    end

    disp([newline, 'INTLINPROG - Objective value (min): ', num2str(obj_intlinprog)]);
    disp(['Time: ', num2str(t_ilp), ' sec']);
    if show_solutions
        disp('Solution x =');
        disp(disp_x_int);
    end
end
