% [MAIN4 - Optimisation continue]
% Cours de H. Ouzia
% Projet : Affectation optimale d'agents � des t�ches
% ========
% 
% A. Khizar, R. Kanyamibwa, M. P�cheux, C. Voisembert
% -----------------------------------------------------------
% Fonction utile pour transformer une matrice en vecteur en
% concat�nant ses lignes les unes � la suite des autres.
% Input(s):
% - M [matrice de r�els] : matrice � transformer
% -----------------------------------------------------------
function v = matrix_to_vector(M)
     v = reshape(M', [], 1);
end
