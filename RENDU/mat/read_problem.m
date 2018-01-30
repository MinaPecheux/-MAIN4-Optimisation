% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% -----------------------------------------------------------
% Fonction lisant les données d'un fichier d'instance fourni
% pour en extraire les matrices et les vecteurs du problème
% (P). Elle retourne également une solution initiale x0
% réalisable.
% Input(s):
% - file [string] : nom du fichier à charger
% -----------------------------------------------------------
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

    % get feasible initial solution
    x0 = sol_initial(m, n, c, A_vec, b, A, Aeq, beq);
end
