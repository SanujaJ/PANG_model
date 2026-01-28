%% example for building the model...

%note: 'tbx' fodler must be added to the path

clear all; close all;

buildOb = buildSystem.buildBase; %call bass class for model building...
buildOb.fileLocation = fileparts(mfilename('fullpath')); %set the file location where a folder with a model-specific items could be created...
buildOb.name = 'SJD_Wing'; %a name for the model

%% some example 

L = 0.655; %semi-span

%% planform geometry...

geom = buildSystem.geom;
geom.L = L;
geom.a = @(x)(0*ones(size(x))); %position of the elastic axis... theodorsen definition
geom.b = @(x)(0.075*ones(size(x))); %semi-chord distribution
geom.sweep = 0*pi/180; %sweep (of beam line!)
geom.twist = @(x)(0*ones(size(x))); %wing twist distribution

buildOb.geom = geom; %write to buildBass class...

%% basis properties

basis = buildSystem.basis;
basis.Nw = 4; %size of shape basis..out of plane bending
basis.Nv = 4; %.. chord wise bending
basis.Nthet = 4;%...torsion
basis.xi = linspace(0, L, 21); %aerodynamic grid (spanwise)... collocation points assigned at the mid-points of xi

buildOb.basis = basis; %write to buildBass class...

%% elastic properties

% this is an example of a case where we want to change the Young's modulus
% when calling the buit model as a parameter... note that the rigidities
% EI1, etc depends in this..

elas(1) = buildSystem.structure.elasBase; %call elastic property class

%stiffness properties... matrices with these rigidities will be multiplied
%by a user dfined parameters...
elas(1).EI1 = @(x)(2.268*ones(size(x)));
elas(1).EI2 = @(x)(59.136*ones(size(x)));
elas(1).EI12 = @(x)(0);
elas(1).GJ = @(x)(3.38); %GJ handled seperately as we want this to be independednt of E
elas(1).name = 'K_matr'; %namee for this matrix..
elas(1).fctrId = []; %this property, which must be one of that defined among user parameters, will scale this matrix computed for elas(1)

buildOb.elas = elas; %write to buildBass class...

%% add inertia

inertia = buildSystem.structure.inertiaBase; %inertia property class
inertia.e = @(x)(0); %position of mass axis.. theodorsen 
inertia.m = @(x)(0.1242); %mass per unit length [kg/m]
inertia.mxx = @(x)(2.7400e-06); %twisting inertia about chordwise about the CoM [kgm]
inertia.mzz = @(x)(2.6500e-06); %chord-plane rotation inertia [kgm]
inertia.myy = @(x)(0); %rotation inertia..bending rotation [kgm]

for i=1:10
    elem(i) = buildSystem.structure.descrElem; %descrete inertia for the engine....handled as a point mass
    %all properties for this follow the same definitions as that for
    %inertia...but NOT PER UNIT LENGTH AND ARE NOT FUNCTION HANDLES
    elem(i).m = 0.0613;
    elem(i).mxx = 1.1796e-04;
    elem(i).mzz = 1.3197e-04;
    elem(i).myy = 1.8399e-05;
    elem(i).e = 0.05;
    elem(i).xp = 0.0375 + (i-1)*0.0650; %attachment potision...say 30% of the semi-span
end

inertia(1).elem = elem;
inertia(1).name = 'massModel';
inertia(1).fctrId = []; %no scalling is imposed on the mass matrix.. however can be implemented if needed

buildOb.inertia = inertia; %write to buildBass class...

%% gravity..

grav = buildOb.inertia2grav(inertia); % this method in the buildBass class automatically generates a class with gravitational/weigth properties using the assigned inertia properties
buildOb.grav = grav; %write to buildBass class...

%% create the model...

buildOb.prepFolder; %this method sets up the folder and the necesary sub-content to 
buildOb.writeModel;

%% derive analysis models... the buildBass instance is passed in 

run_ONERA = analysis.oneraBase(buildOb); %analysis module foe the embedded-aerodynamic model
run_ext = analysis.extAeroBase(buildOb); %analysis module foe the external aerodynamic model

run_ONERA.Cl_grad = @(alp, U)(aeroCurves.clGradFcn(alp,U));
run_ONERA.Cl = @(alp, U)(aeroCurves.clFcn(alp,U));
run_ONERA.Cm = @(alp, U)(aeroCurves.cmFcn(alp,U));

%save these analysis modules for calling later...
save('run_ONERA.mat', 'run_ONERA')
save('run_ext.mat', 'run_ext')

%%