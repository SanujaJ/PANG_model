function [k_0, k_1, k_2, k_3] = combGravMats_2(p)
% Auto-generated constant matrix
k_0 = project.grav.odr_0.massModel_dir_2(p);
k_1 = project.grav.odr_1.massModel_dir_2(p);
k_2 = project.grav.odr_2.massModel_dir_2(p);
k_3 = project.grav.odr_3.massModel_dir_2(p);
