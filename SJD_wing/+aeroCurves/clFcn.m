function cl = clFcn(alp, U)

fctr = 1;
trg=1;

lifGrad = fctr*(5.975+2.3*exp(-0.4*(U-10)));
kinkPoint = (pi/180)*(0.2139*U-1.756);
clg2 = 2.675;
fac = 50;
da = 0.5*pi/180;


cl = lifGrad.*alp +...
    fctr*trg*0.5*(1+tanh(fac*(alp-kinkPoint-da))).*(...
    clg2*(alp-kinkPoint-da));

end