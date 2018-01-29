% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents à des tâches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. Pêcheux, C. Voisembert
% -----------------------------------------------------------
% Fonction utile pour transformer une matrice en vecteur en
% concaténant ses lignes les unes à la suite des autres.
% Input(s):
% - M [matrice de réels] : matrice à transformer
% -----------------------------------------------------------
function v = matrix_to_vector(M)
     v = reshape(M', [], 1);
end
