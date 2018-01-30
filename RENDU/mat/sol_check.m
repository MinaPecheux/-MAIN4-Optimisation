% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents � des t�ches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. P�cheux, C. Voisembert
% -----------------------------------------------------------
% Fonction lisant les donn�es d'un fichier d'instance fourni
% pour en extraire les matrices et les vecteurs du probl�me
% (P). Elle retourne �galement une solution initiale x0
% r�alisable.
% Input(s):
% - x [vecteur de r�els] : solution � v�rifier
% - A_vec [vecteur de r�els] : vecteur de productivit�
% - b [vecteur de r�els] : sous-vecteur des disponibilit�s
% . correspondant aux contraintes d'in�galit�
% -----------------------------------------------------------
function ok = sol_check(x, A_vec, b)
    disp('Checking solution for feasibility.');
    ok = 1;
    [m, n] = size(x);
    A = vector_to_matrix(m, n, A_vec);
    tol = 1e-3;
    % check weight constraints (3)
    for i = 1:m
        if sum(A(i,:).*x(i,:)) - tol > b(i)
            disp(['Solution is infeasible: weight constraints (3) not ok for agent i=',num2str(i),':']);
            disp(['A*x = ',num2str(sum(A(i,:).*x(i,:))),' > b = ',num2str(b(i))])
            ok = 0;
            return;
        end
    end
    % check equality constraints (4)
    for j = 1:n
        if abs(sum(x(:,j)) - 1) > tol
            disp(['Solution is infeasible: equality constraints (4) not ok for task j=',num2str(j),'.']);
            disp(['Sum of column j = ',num2str(j),': ',num2str(sum(x(:,j)))])
            ok = 0;
            return;
        elseif abs(sum(x(:,j)) - 1) > 0
            x(:,j) = round(x(:,j));
        end
    end
end
