% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents � des t�ches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. P�cheux, C. Voisembert
% -----------------------------------------------------------
% Fonction objectif pour le probl�me (P).
% Input(s):
% - m, n [entiers] : dimensions du probl�me (resp. nombre
% . d'agents et nombre de t�ches)
% - x [vecteur de r�els] : solution r�alisable o� calculer
% . la valeur objectif
% - c [vecteur de r�els] : vecteur des co�ts
% -----------------------------------------------------------
function obj = obj_func(m, n, x, c)
    obj = 0;
    for i = 1:m
        for j = 1:n
            obj = obj + c((i-1)*n+j)*x((i-1)*n+j);
        end
    end
end
