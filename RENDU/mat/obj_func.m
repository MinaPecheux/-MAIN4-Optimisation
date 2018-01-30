% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% -----------------------------------------------------------
% Fonction objectif pour le problème (P).
% Input(s):
% - m, n [entiers] : dimensions du problème (resp. nombre
% . d'agents et nombre de tâches)
% - x [vecteur de réels] : solution réalisable où calculer
% . la valeur objectif
% - c [vecteur de réels] : vecteur des coûts
% -----------------------------------------------------------
function obj = obj_func(m, n, x, c)
    obj = 0;
    for i = 1:m
        for j = 1:n
            obj = obj + c((i-1)*n+j)*x((i-1)*n+j);
        end
    end
end
