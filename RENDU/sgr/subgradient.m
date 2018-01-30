% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% -----------------------------------------------------------
% Script utilisant l'algorithme du sous-gradient pour
% résoudre le problème (P).
% -----------------------------------------------------------
function subgradient()
    % load data
    f = fopen('../dat/1/a05100', 'r');
    r_dimensions = textscan(f, '%f', 2);
    r_values = textscan(f, '%f');
    fclose(f);

    m = r_dimensions{1}(1);
    n = r_dimensions{1}(2);
    nvars = m*n;

    tmp = r_values{1};
    c = tmp(1:nvars);
    a_vec = tmp(nvars+1:2*nvars);
    b = tmp(2*nvars+1:end);

    % compose a matrix from a vector
    a = reshape(a_vec, [n,m])';

    % define parameters used in sub_grad algorithm
    epsilon = 0.1;
    ro = 1;
    pi_zero = zeros(m*m,1);
    iterLimit = 100;
    DualNoChangeTOL = 2;

    tic;
    [x_k, theta_pi_k] = sub_grad( epsilon, ro, pi_zero, iterLimit, DualNoChangeTOL,a,b,c);
    tps = toc;
    
    disp('Objective function value:')
    disp(theta_pi_k)
    disp(['Time: ', num2str(tps), ' sec'])
    disp('x:')
    disp(vector_to_matrix(m, n, x_k(1:m*n)))
    disp(['Feasibility: ', num2str(sol_check(vector_to_matrix(m, n, x_k(1:m*n)), a_vec, b))])
    disp('-------------------------------')
    disp('Cleanup to get integer solution')
    x_k_int = sol_cleanup_ps(m, n, x_k(1:m*n), a_vec, b, c);
    disp('Objective function value:')
    disp(sum(c.*matrix_to_vector(x_k_int)))
    disp('x:')
    disp(x_k_int)
    disp(['Feasibility: ', num2str(sol_check(x_k_int, a_vec, b))])
end

% -----------------------------------------------------------
% Fonction appliquant l'algorithme du sous-gradient au
% problème lu avec les critères d'arrêt fournis.
% Input(s):
% - epsilon [réel] : paramètre de précision
% - ro [réel] : 2e paramètre de précision
% - pi_zero [vecteur de réels] : vecteur pi initial
% - iterLimit [entier] : nombre d'itérations maximal
% - DualNoChangeTOL [entier] : nombre d'itérations tolérées
% . sans changement de la valeur de la fonction duale 
% - a [matrice de réels] : matrice de productivité
% - b [vecteur de réels] : vecteur des disponibilités
% - c [vecteur de réels] : vecteur des coûts
% -----------------------------------------------------------
function [x_k, theta_pi_k] = sub_grad(epsilon, ro, pi_zero, iterLimit, DualNoChangeTOL, a, b, c)
    [m,~] = size(a);

    k = 1;
    t = 1;
    beta_k = -inf;
    pi_k = pi_zero;

    condition = true; 

    while( condition )
        [x_k, theta_pi_k] = solve_lagrangien_relaxation(pi_k,a,b,c);
        gamma_k = g_10(x_k, a, b);

        if theta_pi_k > beta_k
            beta_k = theta_pi_k;
        else
           if t < DualNoChangeTOL
               t = t + 1;
           else
               ro = ro / 2;
               t = 1;           
           end
        end

        if gamma_k == 0
            break;
        else
            theta_chap = applymyHeuristic(x_k, a,b,c);
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

% -----------------------------------------------------------
% Fonction résolvant la relaxation lagrangienne de (L) en
% dualisant les contraintes (10).
% Input(s):
% - pi_k [vecteur de réels] : vecteur pi courant
% - a [matrice de réels] : matrice de productivité
% - b [vecteur de réels] : vecteur des disponibilités
% - c [vecteur de réels] : vecteur des coûts
% -----------------------------------------------------------
function [X, FVAL] = solve_lagrangien_relaxation(pi_k, a, b, c)
    [m,n] = size(a);
    x_length = m*n + m * m * n;

    % start A matrix (constraints (10))
    A_10 = zeros(m*m, x_length);
    row = 1;
    col_z = m*n + 1;
    k = 0;
    for i = 1:m*m
        A_10(i, (row-1)*n + 1 + k) = -b(row) + a(row, k+1);
        A_10(i, col_z:col_z+n-1) = a(row,:);

        k = mod(k+1, m);
        col_z = col_z + n;
        row = floor(i / m) + 1;
    end

    % Compose A_11 matrix for constraints (11) 
%     A_11 = zeros(m*m, x_length);
%     row_a = 1;
%     z_variables = m*n + 1;
%     k = 0;
%     for i = 1:m*m
% 
%         A_11(i, (row_a-1)*n+1:n*row_a) = a(row_a,:);
%         A_11(i,(row_a-1)*n+1 + k) = b(row_a);
% 
%         A_11(i, z_variables:z_variables+n-1) = - a(row_a,:);
%         A_11(i, k+z_variables) = 0;
% 
%         k = mod( k+1, m );
%         z_variables = z_variables + n;
%         row_a = floor( i / m ) + 1;
%     end

    % Compose A_13 matrix for constraints (13) 
    A_13 = zeros( m*m*n ,x_length);
    for i = 1:m*n:m*m*n
        j = (i-1)/m;
        A_13( i: i + m*n -1, j+1: j+n ) = repmat( - eye(n), m, 1);
    end
    A_13(:,m*n+1:end) = eye(m*m*n);

    % Compose A_14 matrix for constraints (14) 
    A_14 = zeros( m*m*n ,x_length);
    j = 1;
    for i = 1:n:m*m*n
        A_14( i : i+n-1, j)= -1;
        j = j + 1; 
    end
    A_14(:,m*n+1:end) = eye(m*m*n);

    % Compose A_15 matrix for constraints (15) 
    A_15 = zeros(m*m*n, x_length);
    A_15(:,m*n+1:end) = -eye(m*m*n);

    pattern = zeros(m*n, n); % matrix that will be repeated inside A_15 matrix
    k = 1;
    for i=1:n:m*n
        pattern( i: i + n -1, k) = ones(n,1);
        pattern( i: m*n+1: end) = pattern( i: m*n+1: end) + 1; 
        k = k + 1;
    end

    for i= 1:m*n:m*m*n
        j = (i-1)/m;
        A_15( i: i + m*n -1, j+1: j+n ) = pattern;
    end
    % remove overlapping (doubled) constraints
    A = [A_10; A_13; A_14; A_15];
    A = unique(A, 'stable', 'rows');

    % compose B matrix (for constraints (11) + (13) + (14) + (15))
    B = [zeros(size(A,1) - m*m*n,1);ones(m*m*n,1)];
%    B(1:m*m) = repelem(b,m,1);

    Aeq = [repmat( eye(n), 1, m) zeros(n, x_length - m*n)];
    Beq = ones(n,1);

    LB = zeros(x_length,1);
    UB = ones(x_length,1);

    nvars = m*n;
    A0 = zeros(m, nvars);
    for i = 1:m
        for j = 1:n
            A0(i, (i-1)*n+j) = a(i,j);
        end
    end
    X0 = [sol_initial(m,n,c,matrix_to_vector(a),b,A0,Aeq(:,1:nvars),Beq);zeros(m*m*n,1)];

    opts = optimoptions('fmincon', 'MaxFunctionEvaluations', 5000);
    [ X, FVAL ] = fmincon( @(x)L_11(x,pi_k,a,b,c), X0, A,B, Aeq, Beq, LB, UB, [], opts);
end

% -----------------------------------------------------------
% Fonction déterminant un majorant theta_chap de theta avec
% une heuristique.
% Input(s):
% - x_0 [vecteur de réels] : point inital pour l'heuristique
% - a [matrice de réels] : matrice de productivité
% - b [vecteur de réels] : vecteur des disponibilités
% - c [vecteur de réels] : vecteur des coûts
% -----------------------------------------------------------
function theta_chap = applymyHeuristic(x_0, a, b, c)
    [m,n] = size(a);

    A = zeros(m, m*n);
    j = 1;
    for i=1:m
        A(i, (j-1)*n+1:j*n ) = a(i,:);
        j = j + 1;
    end

    Aeq = repmat( eye(n), 1, m);
    Beq = ones(n,1);

    LB = zeros(m*n, 1);
    UB = ones(m*n, 1);

    f = @(x) sum(c .* x);

    opts = optimoptions('patternsearch', 'PollMethod', 'GSSPositiveBasisNp1', 'UseCompletePoll', true, 'ScaleMesh', false, 'MeshTolerance', 0.99);

    [~,theta_chap] = patternsearch(f, x_0(1:m*n), A,b, Aeq, Beq, LB, UB, [], opts);
end

% -----------------------------------------------------------
% Fonction objectif contenant les contraintes (10) dualisées.
% Input(s):
% - x [vecteur de réels] : point où calculer la valeur
% . objectif
% - pi [vecteur de réels] : multiplicateur de Lagrange
% . courant
% - a [matrice de réels] : matrice de productivité
% - b [vecteur de réels] : vecteur des disponibilités
% - c [vecteur de réels] : vecteur des coûts
% -----------------------------------------------------------
function L = L_10(x, pi, a, b, c)
    % define the Lagrangien relaxation with constraint (10)
    [m,n] = size(a);
    C = [ c; zeros(m*m*n,1) ];

    g = g_10(x,a,b);

    L = sum(C.*x) + sum(pi .* g);
end

% -----------------------------------------------------------
% Fonction renvoyant la valeur des fonctions dualisées avec
% les contraintes (10) au point donné.
% Input(s):
% - x [vecteur de réels] : point où calculer la valeur
% . objectif
% - a [matrice de réels] : matrice de productivité
% - b [vecteur de réels] : vecteur des disponibilités
% -----------------------------------------------------------
function g = g_10(x, a, b)
    % define g vector containing all constraints represented by (10)
    [m,n] = size(a);

    g = zeros(m*m,1);
    z_variables = m*n + 1;

    for i = 1:m
        for k=1:m
            g((i-1)*m+k) = sum( a(i,:)' .* x( z_variables: z_variables + n - 1) ) - (b(i) - a(i,k))*x((i-1) *m + k);
            z_variables = z_variables + n;
        end
    end
end

% -----------------------------------------------------------
% Fonction objectif contenant les contraintes (11) dualisées.
% Input(s):
% - x [vecteur de réels] : point où calculer la valeur
% . objectif
% - pi [vecteur de réels] : multiplicateur de Lagrange
% . courant
% - a [matrice de réels] : matrice de productivité
% - b [vecteur de réels] : vecteur des disponibilités
% - c [vecteur de réels] : vecteur des coûts
% -----------------------------------------------------------
function L = L_11(x, pi, a, b, c)
    % define the Lagrangien relaxation with constraint (11)
    [m,n] = size(a);
    c = [c;zeros(m*n*m,1)];

    g = g_11(x,a,b);
    
    L = sum(c.*x) + sum(pi.*g);
end

% -----------------------------------------------------------
% Fonction renvoyant la valeur des fonctions dualisées avec
% les contraintes (11) au point donné.
% Input(s):
% - x [vecteur de réels] : point où calculer la valeur
% . objectif
% - a [matrice de réels] : matrice de productivité
% - b [vecteur de réels] : vecteur des disponibilités
% -----------------------------------------------------------
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
