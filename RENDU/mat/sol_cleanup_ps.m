% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% ----------------------------------------------------------
% Fonction nettoyant une solution trouvée par patternsearch
% qui peut contenir des colonnes non entières.
% Input(s):
% - m, n [entiers] : dimensions du problème (resp. nombre
% . d'agents et nombre de tâches)
% - dirty_x [vecteur de réels] : solution trouvée pour (P)
% . par patternsearch (à nettoyer)
% - A_vec [vecteur de réels] : vecteur de productivité
% - c [vecteur de réels] : vecteur des coûts
% -----------------------------------------------------------
function x = sol_cleanup_ps(m, n, dirty_x, A_vec, c)
    % get A and c as matrix
    A_mat = vector_to_matrix(m, n, A_vec);
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
        col = round(col, 4);
        % if column has non-integer values, clean it
        if any(col ~= floor(col) > 0)
            % get matching ratio column
            r = ratios(:, j);
            % get first min index
            ind = find(r == min(r), 1, 'first');
            % apply 1 coeff to x vector at the pos corresponding to r min
            % value
            x(ind, j) = 1;
        % else just copy the column
        else
            x(:,j) = dirty_x_mat(:,j);
        end
    end
end
