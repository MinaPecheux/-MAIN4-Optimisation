#------------------------------------------
#
# Based on:
# (c) Hacene Ouzia, Polytech'Paris UPMC
#
# Updated by:
# A. Khizar, R. Kanyamibwa, M. Pecheux
# C. Voisembert
#
#------------------------------------------
# Exemple de modele avec gmpl
#------------------------------------------

#--- Nombre de variables
param M, integer, > 0;
param N, integer, > 0;

#--- Indices des colonnes
set I := 1..M;

#--- Indices des lignes
set J := 1..N;

#--- Matrice des contraintes
param a{i in I, j in J};

#--- Second membre
param b{i in I};

#--- Coefficients de la fonc. objectif
param c{i in I, j in J};

#--- Variables de decision
var x{i in I, j in J} >= 0;
var z{i in I, j in J, k in I} >= 0;

#--- Contraintes du probleme
s.t. tasks_to_agents{j in J}: sum{i in I}x[i,j] = 1;
# (P)
s.t. agents_capacities{i in I}: sum{j in J}a[i,j]*x[i,j] <= b[i];
# (M)
s.t. agents_capacities_M{i in I, k in I}: sum{j in J}a[i,j]*x[i,j]*x[i,k] <= (b[i] - a[i,k])*x[i,k];
s.t. cap_M{i in I, k in I}: sum{j in J : j<>k}a[i,j]*x[i,j] - sum{j in J : j<>k}a[i,j]*x[i,j]*x[i,k] <= b[i] - b[i]*x[i,k];
# (L)
s.t. agents_capacities_lin{i in I, k in I}: sum{j in J}a[i,j]*z[i,j,k] <= (b[i] - a[i,k])*x[i,k];
s.t. cap_lin{i in I, k in I}: sum{j in J : j<>k}a[i,j]*x[i,j] - sum{j in J : j<>k}a[i,j]*z[i,j,k] <= b[i] - b[i]*x[i,k];
s.t. lin_const1{i in I, j in J, k in I}: z[i,j,k] <= x[i,j];
s.t. lin_const2{i in I, j in J, k in I}: z[i,j,k] <= x[i,k];
s.t. lin_const3{i in I, j in J, k in I}: x[i,k] + x[i,j] - z[i,j,k] <= 1;

#--- Critere a optimiser
minimize Objective: sum{i in I, j in J} c[i,j]*x[i,j];
