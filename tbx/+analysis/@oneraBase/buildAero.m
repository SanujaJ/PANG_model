function [M_glob,F_glob] = buildAero(runObj,M_glob,F_glob,q_all,p0,p)

beta = runObj.geom.sweep;

Minf = p0(6)*cos(beta);
mach_bet = sqrt(1-Minf^2);


%%onera stuff....
SL = pi;
lam = runObj.lam;
ML = runObj.ML;
sigM = -0.25*pi;
sigL = 2*pi; 
sM = -0.5890;
KL = 0.5*pi + 1.96*pi*(mach_bet-1);
sigM_bar = -0.25*pi*(1 + 1.4*Minf.^2);


%% runObj...

q = runObj.transF*q_all(runObj.structDisp);
qt = runObj.transF*q_all(runObj.structVel);
qa = q_all(runObj.unsAeroIdx);

if isempty(runObj.unsCircIdx)   
else
    qgam = q_all(runObj.unsCircIdx); %unsteady circulations are explicitly tracked....propegation eq is ON
end

if nargin<5
    p = [];
end

%%

%geometry...semi-chords.,,
bi = project.aero.geom.bi_fcn(p);

%flow components...chordwise...
vy = analysis.aero.flow.Vy(q,qt,p0,p,beta);
vyBlk = diag(vy);

%.... W0/W1 and ders..
w0 = analysis.aero.flow.W0(q,qt,p0,p,beta);
w1 = analysis.aero.flow.W1(q,qt,p0,p,beta);
[w0_dt_LHS, w0_dt_RHS] = analysis.aero.flow.W0_dt(q,qt,p0,p,beta);
[w1_dt_LHS, w1_dt_RHS] = analysis.aero.flow.W1_dt(q,qt,p0,p,beta);

indFac = project.aero.flow.u_ind(p);

% total velocity in z
vz = w0 + 0.5*w1 - mach_bet*project.aero.flow.u_ind(p)*qgam;
vzBlk = diag(vz); %as a diagonal operator..

%total airspeed...
U_tot = sqrt(vz.^2 + vy.^2);
    U_Mat = diag(U_tot);
    invU = inv(U_Mat);
    

%% get all lift components...

[Cl, Cl_alp] = runObj.Cl_fcn(vz/mach_bet,U_tot);
[Cm, ~] = runObj.Cm_fcn(vz/mach_bet,U_tot);
[Cd, ~] = runObj.Cd_fcn(vz/mach_bet,U_tot);

%% get displacement functions...

proj_lift = analysis.aero.disps.liftDisp(q,qt,p0,p,vyBlk,vzBlk);
proj_moment = analysis.aero.disps.momentDisp(q,qt,p0,p,vyBlk,vzBlk);
proj_drag = analysis.aero.disps.dragDisp(q,qt,p0,p,vyBlk,vzBlk);

%% aero forces on structure.. 

totLifStr_R = p0(5)*qa +...
    p0(5)*(invU)*(bi.^2)*(SL*w0_dt_RHS + KL*w1_dt_RHS);

totMoM_R = 2*p0(5)*(bi.^2)*(...
    bi*(sigM_bar*w0_dt_RHS+sM*w1_dt_RHS)+...
    U_Mat*sigM*w1+(U_Mat.^2)*Cm);

totDrg_R = p0(5)*bi*U_Mat*Cd;

totLifStr_L = ...
    -p0(5)*(invU)*(bi.^2)*(SL*w0_dt_LHS + KL*w1_dt_LHS);

totMoM_L = -2*p0(5)*(bi.^2)*(bi*(sigM_bar*w0_dt_LHS+sM*w1_dt_LHS));

%% unsteady aero eqns..

RHS_pot = -lam*qa +...
    lam*bi*(Cl.*U_tot + sigL.*w1) +...
    (invU)*ML*(bi.^2)*(sigL.*w0_dt_RHS +...%sigL this line origanlly is lifgrad - same
    sigL.*w1_dt_RHS);

LHSpot_str = -(invU)*ML*(bi.^2)*(sigL.*w0_dt_LHS +... %sigL this line origanlly is lifgrad - same
    sigL.*w1_dt_LHS);

LHSpot_aer = (invU)*bi;

dwnWsh_pot = -(invU)*ML*(bi.^2)*(-Cl_alp.*indFac);

%% apply structural aerodynamic forces to structure...

F_glob(runObj.structVel) = F_glob(runObj.structVel) + runObj.transF'*(...
    proj_lift*totLifStr_R + proj_moment*totMoM_R + proj_drag*totDrg_R);

F_glob(runObj.unsAeroIdx) = F_glob(runObj.unsAeroIdx) + RHS_pot; %locally 2D unsteady comps

M_glob(runObj.structVel, runObj.structVel) = M_glob(runObj.structVel, runObj.structVel) + runObj.transF'*(...
    proj_lift*totLifStr_L + proj_moment*totMoM_L)*runObj.transF;

M_glob(runObj.unsAeroIdx, runObj.structVel) = M_glob(runObj.unsAeroIdx, runObj.structVel) +...
    LHSpot_str*runObj.transF;

M_glob(runObj.unsAeroIdx, runObj.unsAeroIdx) = M_glob(runObj.unsAeroIdx, runObj.unsAeroIdx) + LHSpot_aer;

if isempty(runObj.unsCircIdx)
else
    
%     M_glob(runObj.structVel, runObj.unsCircIdx) =  M_glob(runObj.structVel, runObj.unsCircIdx) +...
%         runObj.transF'*[proj_lift*p0(5)*inv(U_Mat)*(bi.^2)*(SL*indFac)+proj_mom*2*p0(5)*(bi.^2)*(bi*(sigM_bar*indFac))];

    M_glob(runObj.unsAeroIdx, runObj.unsCircIdx) =  M_glob(runObj.unsAeroIdx, runObj.unsCircIdx) + dwnWsh_pot;

    M_glob(runObj.unsCircIdx, runObj.unsCircIdx) = bi*(invU);

    F_glob(runObj.unsCircIdx) =  F_glob(runObj.unsCircIdx) + qa-qgam;
end

end

















