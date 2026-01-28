function [k_1, k_2, k_3] = combGyrMatr(p)
% Auto-generated constant matrix
k_1 = project.inertia.odr_1.gyr_massModel(p);
k_2 = project.inertia.odr_2.gyr_massModel(p);
k_3 = project.inertia.odr_3.gyr_massModel(p);
