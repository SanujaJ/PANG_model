function Cm = cmFcn(alp, U)

fctr = 1;

kinkPoint = (pi/180)*(0.2139*U-1.756);
cmg2 = -fctr*0.675;
fac = 50;
da = 0.5*pi/180;

Cm_grad = fctr*(0.2218-0.0038*(U-16.25)-0.25*exp(-0.5*(U-10)));
Cm = Cm_grad.*alp +...
    0.5*(1+tanh(fac*(alp-kinkPoint-da))).*(...
    cmg2*(alp-kinkPoint-da));



end