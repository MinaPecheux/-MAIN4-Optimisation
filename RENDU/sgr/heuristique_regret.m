% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% -----------------------------------------------------------
% Fonction cherchant une solution initiale réalisable pour le
% problème (P) (sans garantie).
% Input(s):
% - m, n [entiers] : dimensions du problème (resp. nombre
% . d'agents et nombre de tâches)
% - c [vecteur de réels] : vecteur des coûts
% - A_vec [vecteur de réels] : vecteur de productivité
% - b [vecteur de réels] : sous-vecteur des disponibilités
% . correspondant aux contraintes d'inégalité
% -----------------------------------------------------------
function [sol, realisable, zeta] = heuristique_regret(m, n, c, A_vec, b)
    % get A and c as matrix
    A_mat = vector_to_matrix(m, n, A_vec);
    c_mat = vector_to_matrix(m, n, c);
    % compute preferences matrix
    f = -c_mat./A_mat;
    
    % initialisation des variables
    realisable = 1;
    beta = b;
    I = 1:m;
    J = 1:n;
    x = zeros(n, 1);
    zeta = 0;
    % find a first best assignment
    while ~isempty(J) && realisable == 1
        regmax = -Inf;
        for j = J
            Rj = I(A_mat(:,j) <= beta);
            if isempty(Rj)
                realisable = 0;
            else
                [~, k] = max(f(Rj,j));
                if isempty(Rj(Rj ~= k))
                    reg = Inf;
                else
                    reg = f(k,j) - max(f(Rj(Rj ~= k),j));
                end
                if reg > regmax
                    regmax = reg;
                    imax = k;
                    jmax = j;
                end
            end
        end
        if realisable == 1
            x(jmax) = imax;
            zeta = zeta - c_mat(imax, jmax);
            beta(imax) = beta(imax) - A_mat(imax, jmax);
            J = J(J ~= jmax);
        end
    end
    
    sol = zeros(m,n);
    if realisable == 1
        for j = 1:n
            % get index of the best agent for task j found yet
            k = x(j);
            % compute indices i matching the two conditions:
            % - i != k
            % - agent i can perform the task j: a_ij <= beta_i
            I_bis = [];
            last_i = 1;
            for i = 1:m
                if i ~= k && A_mat(i,j) <= beta(i)
                    I_bis(last_i) = i;
                    last_i = last_i + 1;
                end
            end
            % get corresponding utilities
            K = -c_mat(I_bis,j);
            % if at least one other agent can perform the task j,
            if ~isempty(K)
                % get the best utility from all the possible replacements
                [minval, qtmp] = min(K);
                % if it is better than the current one, reassign task j
                if minval < -c_mat(k,j)
                    q = I_bis(qtmp);
                    x(j) = q;
                    zeta = zeta + c_mat(k,j) - c_mat(q,j);
                    beta(k) = beta(k) + A_mat(k,j);
                    beta(q) = beta(q) - A_mat(q,j);
                end
            end
        end
    
        % create matrix solution
        for j = 1:n
            sol(x(j),j) = 1;
        end
    end
end
