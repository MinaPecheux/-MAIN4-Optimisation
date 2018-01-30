% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% -----------------------------------------------------------
% Fonction nettoyant une solution trouvée par ga qui peut
% contenir des colonnes ne respectant pas les contraintes
% d'égalité (4) (un seul 1 par colonne).
% Input(s):
% - m, n [entiers] : dimensions du problème (resp. nombre
% . d'agents et nombre de tâches)
% - dirty_x [vecteur de réels] : solution trouvée pour (P)
% . par patternsearch (à nettoyer)
% - A_vec [vecteur de réels] : vecteur de productivité
% - c [vecteur de réels] : vecteur des coûts
% -----------------------------------------------------------
function x = sol_cleanup_ga(m, n, dirty_x, A_vec, c)
    % get A and c as matrix
    A_mat = vector_to_matrix(m, n, A_vec);
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
        col = round(col, 4);
        % get all non-zero values
        col_nonzeros_ind = find(col);
        % if there is none, clean it
        if isempty(col_nonzeros_ind)
            % get matching ratio columns at non-zero values indexes
            r = ratios(col_nonzeros_ind, j);
            % apply 1 coeff to x vector at the pos corresponding to r max
            % value
            x(r == max(r), j) = 1;
        % if there are more than one, clean it
        elseif length(col_nonzeros_ind) > 1
            % get matching ratio columns at non-zero values indexes
            r = ratios(col_nonzeros_ind, j);
            % apply 1 coeff to x vector at the pos corresponding to r max
            % value
            x(col_nonzeros_ind(r == max(r)), j) = 1;
        % else just copy the column
        else
            x(:,j) = dirty_x_mat(:,j);
        end
    end
end
