% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% -----------------------------------------------------------
% Fonction utile pour transformer un vecteur en matrice de
% m lignes et n colonnes.
% Input(s):
% - m, n [entiers] : dimensions de la matrice retournée
% - v [vecteur de réels] : vecteur à transformer
%
% Remarque : si l'une des dimensions vaut 1, on renvoie
% simplement le vecteur "dans le bon sens" (colonne ou ligne)
% -----------------------------------------------------------
function M = vector_to_matrix(m, n, v)
    if m == 1
        if isrow(v)
            M = v;
            return;
        else
            M = v';
            return;
        end
    elseif n == 1
        if isrow(v)
            M = v';
            return;
        else
            M = v;
            return;
        end
    end

    M = zeros(m,n);
    for i = 1:m
        M(i,:) = v((i-1)*n+1:(i-1)*n+n);
    end
end