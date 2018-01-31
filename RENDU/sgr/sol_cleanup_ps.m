% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents � des t�ches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. P�cheux, C. Voisembert
% ----------------------------------------------------------
% Fonction nettoyant une solution trouv�e par patternsearch
% qui peut contenir des colonnes non enti�res.
% Input(s):
% - m, n [entiers] : dimensions du probl�me (resp. nombre
% . d'agents et nombre de t�ches)
% - dirty_x [vecteur de r�els] : solution trouv�e pour (P)
% . par patternsearch (� nettoyer)
% - A_vec [vecteur de r�els] : vecteur de productivit�
% - b [vecteur de r�els] : vecteur des disponibilit�s
% - c [vecteur de r�els] : vecteur des co�ts
% -----------------------------------------------------------
function x = sol_cleanup_ps(m, n, dirty_x, A_vec, b, c)
    % get A and c as matrix
    A_mat = vector_to_matrix(m, n, A_vec);
    c_mat = vector_to_matrix(m, n, c);
    % compute ratio matrix
    ratios = -c_mat./A_mat;
    checked_ratios = zeros(m, n);
    
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
            checked_ratios(ind, j) = 1;
        % else just copy the column
        else
            x(:,j) = dirty_x_mat(:,j);
        end
    end
    
    tol = 1e-2;
    % check weight constraints (3)
    for i = 1:m
        if sum(A_mat(i,:).*x(i,:)) - tol > b(i)
            % get non-zero values in row
            row_nonzeros_ind = find(x(i,:));
            for k = row_nonzeros_ind
                % get matching ratio column
                r = ratios(:, k);
                % get matching ratio check
                check_r = checked_ratios(:, k);
                % get first non-checked min index
                ind = find(r == min(r(check_r ~= 1)), 1, 'first');
                % check for new solution feasibility on found row
                new_row = x(ind,:);
                new_row(k) = 1;
                if sum(A_mat(ind,:).*new_row) - tol <= b(ind)
                    % apply 1 coeff to x vector at the pos corresponding to r min
                    % value
                    x(i, k) = 0;
                    x(ind, k) = 1;
                    checked_ratios(ind, k) = 1;
                
                    % check for new solution feasibility on checked row
                    if sum(A_mat(i,:).*x(i,:)) - tol <= b(i)
                        break;
                    end
                end
            end
        end
    end
end
