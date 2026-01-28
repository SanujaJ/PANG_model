function lifGrad = clGradFcn(alp, U)

fctr = 1;

lifGrad = fctr*(5.975+2.3*exp(-0.4*(U-10)));

end