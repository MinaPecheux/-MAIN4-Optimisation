% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents � des t�ches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. P�cheux, C. Voisembert
% -----------------------------------------------------------
% Fonction renvoyant une solution initiale r�alisable pour le
% probl�me (P). Elle essaie d'abord l'algorithme de regret
% puis, si la solution trouv�e n'est pas r�alisable, utilise
% intlinprog en l'interrompant � la premi�re solution trouv�e.
% Input(s):
% - m, n [entiers] : dimensions du probl�me (resp. nombre
% . d'agents et nombre de t�ches)
% - c [vecteur de r�els] : vecteur des co�ts
% - A_vec [vecteur de r�els] : vecteur de productivit�
% - b [vecteur de r�els] : sous-vecteur des disponibilit�s
% . correspondant aux contraintes d'in�galit�
% - A [matrice de r�els] : sous-matrice de productivit�
% . correspondant aux contraintes d'in�galit�
% - Aeq [matrice de r�els] : sous-matrice de productivit�
% . correspondant aux contraintes d'�galit�
% - beq [vecteur de r�els] : sous-vecteur des disponibilit�s
% . correspondant aux contraintes d'�galit�
% -----------------------------------------------------------
function x0 = sol_initial(m, n, c, A_vec, b, A, Aeq, beq)
    [x0_mat, ~] = heuristique_regret(m, n, c, A_vec, b);
    x0 = matrix_to_vector(x0_mat);

    % if x0 is not feasible:
    % find initial feasible solution with interrupted intlinprog
    if sol_check(x0_mat, A_vec, b) == 0
        opts = optimoptions('intlinprog', 'MaxFeasiblePoints', 1);
        nvars = m*n;
        lb = zeros(nvars,1);
        ub = ones(nvars,1);
        [x0, ~] = intlinprog(c, 1:nvars, A, b, Aeq, beq, lb, ub, [], opts);
    end
end
