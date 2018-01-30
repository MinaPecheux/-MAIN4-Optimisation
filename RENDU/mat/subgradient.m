% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% ---------------------------------------------------
% Script utilisant l'algorithme du sous-gradient pour
% résoudre le problème (P).
% ---------------------------------------------------

function subgradient()
    % load data
    f = fopen('../data-1/1/a0507', 'r');
    r_dimensions = textscan(f, '%f', 2);
    r_values = textscan(f, '%f');
    fclose(f);

    m = r_dimensions{1}(1);
    n = r_dimensions{1}(2);
    nvars = m*n;

    tmp = r_values{1};
    c = tmp(1:nvars);
    A_vec = tmp(nvars+1:2*nvars);
    b = tmp(2*nvars+1:end);

    % compose a matrix from a vector
    a = vector_to_matrix(m, n, A_vec);

    % define parameters used in sub_grad algorithm
    epsilon = 0.1;
    ro = 1;
    pi_zero = zeros(m*m,1);
    iterLimit = 100;
    DualNoChangTOL = 2;

    [x_k, theta_pi_k] = sub_grad(epsilon, ro, pi_zero, iterLimit, DualNoChangTOL,a,b,c);

    disp('Objective function value:')
    disp(theta_pi_k)
    disp('x:')
    disp(vector_to_matrix(m, n, x_k(1:m*n)))
    disp('----------------------------------')
    disp('Clean up to get integer solutions:')
    disp('x:')
    disp(sol_cleanup_ps(m, n, x_k(1:m*n), A_vec, c));
end


function [x_k, theta_pi_k] = sub_grad( epsilon, ro, pi_zero, iterLimit, DualNoChangTOL,a,b,c)
    [m,~] = size(a);

    k = 1;
    t = 1;
    beta_k = -inf;
    pi_k = pi_zero;

    condition = true; 

    while(condition)
        [x_k, theta_pi_k] = solve_lagrangien_relaxation(pi_k,a,b,c);
        gamma_k = g_10(x_k, a, b);

        if theta_pi_k > beta_k
            beta_k = theta_pi_k;
        else
           if t < DualNoChangTOL
               t = t + 1;
           else
               ro = ro / 2;
               t = 1;           
           end
        end

        if gamma_k == 0
            break;
        else
            theta_chap = applyMyHeuristic(x_k, a,b,c);
            if (isempty(theta_chap))
                break;
            end

            for j = 1:m*m
                pi_k(j) = max(0, pi_k(j) - (theta_pi_k - theta_chap)/norm(gamma_k)*ro*gamma_k(j) );
            end
            k = k + 1;
            condition = abs(L_10(x_k,pi_k, a, b, c)-theta_pi_k) / theta_pi_k > epsilon && k<iterLimit;
        end

    end
end

function [X, FVAL] = solve_lagrangien_relaxation(pi_k,a,b,c)

    [m,n] = size(a);
    nvars = m*n * (1+m);

    % compose A matrix (for constraints (11) + (13) + (14) + (15))
    A = zeros(m*m + 3*m*m*n, nvars);
                            % SIZE OF THIS MATRIX:
                            % (m*m rows because there are constraints for i,k in {1,...m}
                            % + m*m*n for constraints (13) + m*m*n for constraints(14))
                            % + m*m*n for constraints (15))
    % start A matrix (constraints (11))
    row = 1;
    col_z = m*n + 1;
    k = 0;
    for i = 1:m*m
        A(i, (row-1)*n+1:n*row) = a(row,:);
        A(i, (row-1)*n + 1 + k) = b(row);
        A(i, col_z:col_z+n-1) = -a(row,:);
        A(i, k + col_z) = 0;

        k = mod(k+1, m);
        col_z = col_z + n;
        row = floor(i / m) + 1;
    end
    row = row + m;
    
    % start A matrix (constraints (10))
    %row = 1;
    %col_z = m*n + 1;
    %k = 0;
    %for i = 1:m*m
        %A(i, (row-1)*n + 1 + k) = -b(row) + a(row, k+1);
        %A(i, col_z:col_z+n-1) = a(row,:);

        %k = mod(k+1, m);
        %col_z = col_z + n;
        %row = floor(i / m) + 1;
    %end
    %row = row + m;
    
    % complete A matrix (constraints (13))
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
    % complete A matrix (constraints (14))
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
    B(1:m*m) = repelem(b,m,1);

    % compose Aeq matrix (for constraint (4) on column sum)
    Aeq = zeros(n, nvars);
    Aeq(:,1:m*n) = repmat(eye(n), 1, m);
    beq = ones(n, 1);

    lb = zeros(nvars, 1);
    ub = [ones(m*n, 1);Inf*ones(nvars - m*n, 1)];

    nvars = m*n;
    A0 = zeros(m, nvars);
    for i = 1:m
        for j = 1:n
            A0(i, (i-1)*n+j) = a(i,j);
        end
    end
    X0 = [sol_initial(m,n,c,matrix_to_vector(a),b,A0,Aeq(:,1:nvars),beq);zeros(m*m*n,1)];

    opts = optimoptions('fmincon', 'Display', 'iter');
    [X, FVAL] = fmincon(@(x)L_10(x,pi_k,a,b,c), X0, A,B, Aeq, beq, lb, ub);
end

function theta_chap = applyMyHeuristic(x_0,a,b,c)
    [m,n] = size(a);

    A = zeros(m, m*n);
    j = 1;
    for i=1:m
        A(i, (j-1)*n+1:j*n ) = a(i,:);
        j = j + 1;
    end

    Aeq = repmat(eye(n), 1, m);
    Beq = ones(n,1);

    LB = zeros(m*n, 1);
    UB = ones(m*n, 1);

    f = @(x) obj_func(m, n, x, c);

    opts = optimoptions('patternsearch', 'PollMethod', 'GSSPositiveBasisNp1', 'UseCompletePoll', true, 'ScaleMesh', false, 'MeshTolerance', 0.99);

    [~,theta_chap] = patternsearch(f, x_0(1:m*n), A,b, Aeq, Beq, LB, UB, [], opts);

end

function L = L_10(x, pi, a, b, c)
    % define the Lagrangien relaxation with constraint (10)
    [m,n] = size(a);
    c = [c;zeros(m*n*m,1)];

    g = g_10(x,a,b);
    
    L = sum(c.*x) + sum(pi.*g);
end

function L = L_11(x, pi, a, b, c)
    % define the Lagrangien relaxation with constraint (11)
    [m,n] = size(a);
    c = [c;zeros(m*n*m,1)];

    g = g_11(x,a,b);
    
    L = sum(c.*x) + sum(pi.*g);
end

function g = g_10(x, a, b)
    [m,n] = size(a);
    
    A = zeros(m*m, length(x));
    row = 1;
    col_z = m*n + 1;
    k = 0;
    for i = 1:m*m
        A(i, (row-1)*n + 1 + k) = -b(row) + a(row, k+1);
        A(i, col_z:col_z+n-1) = a(row,:);

        k = mod(k+1, m);
        col_z = col_z + n;
        row = floor(i / m) + 1;
    end
    
    g = A * x;
end

function g = g_11(x, a, b)
    [m,n] = size(a);
    
    A = zeros(m*m, length(x));
    row = 1;
    col_z = m*n + 1;
    k = 0;
    for i = 1:m*m
        A(i, (row-1)*n+1:n*row) = a(row,:);
        A(i, (row-1)*n + 1 + k) = b(row);
        A(i, col_z:col_z+n-1) = -a(row,:);
        A(i, k + col_z) = 0;

        k = mod(k+1, m);
        col_z = col_z + n;
        row = floor(i / m) + 1;
    end
        
    g = A * x;
end
