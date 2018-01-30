% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% -----------------------------------------------------------
% Fonction renvoyant une solution initiale réalisable pour le
% problème (P). Elle essaie d'abord l'algorithme de regret
% puis, si la solution trouvée n'est pas réalisable, utilise
% intlinprog en l'interrompant à la première solution trouvée.
% Input(s):
% - m, n [entiers] : dimensions du problème (resp. nombre
% . d'agents et nombre de tâches)
% - c [vecteur de réels] : vecteur des coûts
% - A_vec [vecteur de réels] : vecteur de productivité
% - b [vecteur de réels] : sous-vecteur des disponibilités
% . correspondant aux contraintes d'inégalité
% - A [matrice de réels] : sous-matrice de productivité
% . correspondant aux contraintes d'inégalité
% - Aeq [matrice de réels] : sous-matrice de productivité
% . correspondant aux contraintes d'égalité
% - beq [vecteur de réels] : sous-vecteur des disponibilités
% . correspondant aux contraintes d'égalité
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
