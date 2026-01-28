clear all; clc;


load('run_ONERA.mat'); run = run_ONERA;
run = run.setPars('alpha0', 0*pi/180, 'alpha', 0, 'g', 9.81, 'mach', 0);

run.dispFlag = false; %this suppresses displayed messages

%get linear modes and matrices....
[shp, evals, Kmat, Cmat, Mmat] = run.getStructModes(run.q0_struct); 

%add damping
Dmat = 0.07*Mmat + (1e-4)*Kmat;
run.dampMatr = Dmat; 

%set to work with modal form...
run = run.setTransform('modal', 6); %...using the first 6 modes

%% equilirbium tracing.. using coco...

%startong points...
alp0 = (2-0.8)*pi/180;
U0 =  11;

%initial aerelastic solution at U0, alp0
q0 = fsolve(@(q_all)(run.aero_structDeriv(q_all, 'U', 11, 'alpha', alp0)),...
    run.q0_aeroStruct);

%ep tool box problem structure...
%see coco po example...
prob_equib = coco_prob();
prob_equib = coco_set(prob_equib,'cont','PtMX', [200,200]);
prob_equib = coco_set(prob_equib,'ode','vectorized', false);


%see coco po example....
prob_equib = ode_isol2ep(prob_equib, '',...
    @(q,p)(...
    run.aero_structDeriv(q, 'U', p(1), 'alpha', p(2), 'alpha0', p(2))),...
    q0,...
       {'U', 'alp'}, [11, alp0]);

%runner call...see coco 
bd_TrvEquib = coco(prob_equib, 'equilibrium_branch', [],...
    1, {'U'}, [10, 31]);

%.read coco solutions....see coco examples...
equib = coco_bd_read('equilibrium_branch');
q_stat = coco_bd_col(equib, 'x');
U_stat = coco_bd_col(equib, 'U');


%plot equilibrium path....
for i=1:length(q_stat(1,:))
    [x,y,z] = run.getDisplField(q_stat(:,i), run.geom.L,...
        'beamModel');
    tipDisp_equib(i) = z(1,1)/run.geom.L;
end

figure; 
xlabel('Airspeed, [U]'); ylabel('w_L/L, [-]'); hold on;
plot(U_stat, tipDisp_equib, 'k-'); hold on;

%% co-dim2 HB tracing...

prob_HB = coco_prob();
prob_HB = coco_set(prob_HB,'cont','PtMX', [-500,500],...
    'ResTOL',1e-1, 'TOL', 5e1, 'h_min', 10,'h_max', 1000);
prob_HB = coco_set(prob_HB,'corr','PtMX', [0,2000], 'ItMX', 250,...
    'ResTOL',1e-1, 'TOL', 5e1,'h_max', 1000);

HB_labs = coco_bd_labs(bd_TrvEquib, 'HB');
prob_HB = ode_HB2HB(prob_HB, '', 'equilibrium_branch', HB_labs(1));

%runner call...Ip
bd_HB = coco(prob_HB, 'Flutter_boundary', [], 1, {'U', 'alp'},...
    {[12.5, 35],[0.8, 5]*pi/180-0.8*pi/180});

%% .collect alpha-U results...

HBEquib = coco_bd_read('Flutter_boundary');

U_flut = coco_bd_col(HBEquib, 'U');
alp_flut = coco_bd_col(HBEquib, 'alp');
q_flut = coco_bd_col(HBEquib, 'x');

%get displacements...
for i=1:length(q_flut(1,:))
    [x,y,z] = run.getDisplField(q_flut(:,i), run.geom.L,...
        'beamModel');
    tipDisp_flut(i) = z(1,1)/run.geom.L;
end

hold on; 
plot(U_flut, tipDisp_flut, 'r-', 'linewidth', 2);
