function [k_1, k_2, k_3] = combElasMatr(p)
% Auto-generated constant matrix
k_1 = project.elas.odr_1.K_matr(p);
k_2 = project.elas.odr_2.K_matr(p);
k_3 = project.elas.odr_3.K_matr(p);
