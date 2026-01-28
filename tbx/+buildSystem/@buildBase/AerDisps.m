function  obj = AerDisps(obj)

disp([])
disp('Running mass matrix computations......................')
disp('  Extracting saved parameters and basis fcns...')

fPath = ['+project\+aero\+disps'];

%extract attachment positions from EoMs...
a = obj.geom.a; b = obj.geom.b;
WT = obj.geom.WTWidth;
beta = obj.geom.sweep;

%open saved basis info..
N_tot = obj.basis.Ntot;

Aer = obj.aer;
xGrid = Aer.xGrid;
NA = Aer.NA;
N_gam = length(xGrid);
xTrail = Aer.xTrail;
bx = Aer.bAer;
bx = diag(bx);
dx = Aer.dx;
dx = diag(dx);

varAlloc = {@(x)(-project.basis.fv(x)-project.basis.fw(x)) , @(x)(project.basis.V(x)), @(x)(project.basis.W(x)), @(x)(project.basis.V1(x)), @(x)(project.basis.W1(x)), @(x)(project.basis.T(x)),...
    @(x)(project.basis.W(x)), @(x)(project.basis.V(x)), @(x)(project.basis.V1(x)), @(x)(project.basis.W1(x)), @(x)(project.basis.T(x)), @(x)(project.basis.T1(x))};

%extract saved symbolic expressions....
exprs = load('AerDispl_2.mat');
disps = exprs.disps;
varDef = exprs.varDef; parDef = exprs.parDef;
varDef = [disps, varDef];

%create a combined symbolic function for RHS..
Q = sym('q',[N_tot,1]); assume(Q, 'real');
Qt = sym('qt',[N_tot,1]); assume(Qt, 'real');

U = sym('U'); alp = sym('alp'); assume([U,alp], 'real');

disp(['generating aerodynamic force projections...\n'])

%% D_W0 terms...

trms = exprs.AerTrm_Mom;

%assign null matrices...
F{1} = zeros(N_gam,N_tot);
F{2} = zeros(1,N_tot,N_gam,N_tot);
F{3} = zeros(N_tot,N_tot,N_gam,N_tot);

for odr=1:3

    if odr<=min([3,length(trms)])
        for numTrms = 1:length(trms{odr})

            %ident scaling fcators and contributing terms...
            [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
            scalF = @(x)(scalFcn(b(x),a(x)));

            %ident the comtributing shape functionals....
            allocPos = 1;
            for varNum=1:length(varDef)
                pow = polynomialDegree(Var_sym,varDef(varNum));
                if pow>0
                    for powCount=1:pow
                        var{allocPos} = varAlloc{varNum};

                        %record the dimensionality of the assigned vector...
                        if varNum==1, dim{allocPos} = 2; else, dim{allocPos} = 1; end
                        allocPos = allocPos+1;
                    end
                end
            end

            argin.odr = odr; argin.scalF = scalF;
            spaces = zeros(1,odr); % indicator to record free vectors argument spaces..
            for j=1:length(var)
                avail = find(spaces~=1); %identify free spaces..
                if dim{j}==1
                    argin.vctr{avail(1)} = @(x,J)(var{j}(x));
                    spaces(avail(1))=1;
                end
                if dim{j}==2
                    if spaces(1)==0
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(3)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,[avail(1),avail(3)]) = 1; %lock assigmend space
                    else
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(2)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,avail(1:2)) = 1; % lock assigned space
                    end
                end
            end
            Matr0 = matCalc_Aer(argin,N_tot,xGrid); %call integrator...
            clear argin; clear var; clear dim;
            F{odr} = F{odr}  + Matr0;
        end
    end

    writeFunction(F{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], 'Mom_Proj', [])
end


% Fmom_sym = evalMatr_K(F,Q)'*sym(dx);
% Fmom_sym = simplify(Fmom_sym);
% Fmom = matlabFunction(Fmom_sym, 'File','AeroEqns\Mom_Proj','vars', {Q,U,alp});

%%
trms = exprs.AerTrm{1};

%assign null matrices...
F{1} = zeros(N_gam,N_tot);
F{2} = zeros(1,N_tot,N_gam,N_tot);
F{3} = zeros(N_tot,N_tot,N_gam,N_tot);

for odr=1:3%min([3,length(trms)])
    if odr<=min([3,length(trms)])
        for numTrms = 1:length(trms{odr})

            %ident scaling fcators and contributing terms...
            [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
            scalF = @(x)(scalFcn(b(x),a(x)));

            %ident the comtributing shape functionals....
            allocPos = 1;
            for varNum=1:length(varDef)
                pow = polynomialDegree(Var_sym,varDef(varNum));
                if pow>0
                    for powCount=1:pow
                        var{allocPos} = varAlloc{varNum};

                        %record the dimensionality of the assigned vector...
                        if varNum==1, dim{allocPos} = 2; else, dim{allocPos} = 1; end
                        allocPos = allocPos+1;
                    end
                end
            end

            argin.odr = odr; argin.scalF = scalF;
            spaces = zeros(1,odr); % indicator to record free vectors argument spaces..
            for j=1:length(var)
                avail = find(spaces~=1); %identify free spaces..
                if dim{j}==1
                    argin.vctr{avail(1)} = @(x,J)(var{j}(x));
                    spaces(avail(1))=1;
                end
                if dim{j}==2
                    if spaces(1)==0
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(3)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,[avail(1),avail(3)]) = 1; %lock assigmend space
                    else
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(2)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,avail(1:2)) = 1; % lock assigned space
                    end
                end
            end
            Matr0 = matCalc_Aer(argin,N_tot,xGrid); %call integrator...
            clear argin; clear var; clear dim;
            F{odr} = F{odr}  + Matr0;
        end
    end

    writeFunction(F{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], 'Lif_proj1', [])
end


% F_sym = evalMatr_K(F,Q)'*sym(dx);
% F_sym = simplify(F_sym);
%     F = matlabFunction(F_sym, 'File','AeroEqns\Lif_proj1','vars', {Q,U,alp});

%horizontal component...
clear F
trms = exprs.AerTrm{2};

%assign null matrices...
F{1} = zeros(N_gam,N_tot);
F{2} = zeros(1,N_tot,N_gam,N_tot);
F{3} = zeros(N_tot,N_tot,N_gam,N_tot);

for odr=1:3%min([3,length(trms)])
    if odr<=min([3,length(trms)])
        for numTrms = 1:length(trms{odr})

            %ident scaling fcators and contributing terms...
            [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
            scalF = @(x)(scalFcn(b(x),a(x)));

            %ident the comtributing shape functionals....
            allocPos = 1;
            for varNum=1:length(varDef)
                pow = polynomialDegree(Var_sym,varDef(varNum));
                if pow>0
                    for powCount=1:pow
                        var{allocPos} = varAlloc{varNum};

                        %record the dimensionality of the assigned vector...
                        if varNum==1, dim{allocPos} = 2; else, dim{allocPos} = 1; end
                        allocPos = allocPos+1;
                    end
                end
            end

            argin.odr = odr; argin.scalF = scalF;
            spaces = zeros(1,odr); % indicator to record free vectors argument spaces..
            for j=1:length(var)
                avail = find(spaces~=1); %identify free spaces..
                if dim{j}==1
                    argin.vctr{avail(1)} = @(x,J)(var{j}(x));
                    spaces(avail(1))=1;
                end
                if dim{j}==2
                    if spaces(1)==0
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(3)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,[avail(1),avail(3)]) = 1; %lock assigmend space
                    else
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(2)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,avail(1:2)) = 1; % lock assigned space
                    end
                end
            end
            Matr0 = matCalc_Aer(argin,N_tot,xGrid); %call integrator...
            clear argin; clear var; clear dim;
            F{odr} = F{odr}  + Matr0;
        end
    end
    writeFunction(F{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], 'Lif_proj2', [])
end

%
% F_sym = evalMatr_K(F,Q)'*sym(dx);
% F_sym = simplify(F_sym);
% F = matlabFunction(F_sym, 'File','AeroEqns\Lif_proj2','vars', {Q,U,alp});

%% Drag wise terms..

trms = exprs.DrgTrm{1};

clear F
%assign null matrices...
F{1} = zeros(N_gam,N_tot);
F{2} = zeros(1,N_tot,N_gam,N_tot);
F{3} = zeros(N_tot,N_tot,N_gam,N_tot);

for odr=1:3%min([3,length(trms)])
    if odr<=min([3,length(trms)])
        for numTrms = 1:length(trms{odr})

            %ident scaling fcators and contributing terms...
            [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
            scalF = @(x)(scalFcn(b(x),a(x)));

            %ident the comtributing shape functionals....
            allocPos = 1;
            for varNum=1:length(varDef)
                pow = polynomialDegree(Var_sym,varDef(varNum));
                if pow>0
                    for powCount=1:pow
                        var{allocPos} = varAlloc{varNum};

                        %record the dimensionality of the assigned vector...
                        if varNum==1, dim{allocPos} = 2; else, dim{allocPos} = 1; end
                        allocPos = allocPos+1;
                    end
                end
            end

            argin.odr = odr; argin.scalF = scalF;
            spaces = zeros(1,odr); % indicator to record free vectors argument spaces..
            for j=1:length(var)
                avail = find(spaces~=1); %identify free spaces..
                if dim{j}==1
                    argin.vctr{avail(1)} = @(x,J)(var{j}(x));
                    spaces(avail(1))=1;
                end
                if dim{j}==2
                    if spaces(1)==0
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(3)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,[avail(1),avail(3)]) = 1; %lock assigmend space
                    else
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(2)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,avail(1:2)) = 1; % lock assigned space
                    end
                end
            end
            Matr0 = matCalc_Aer(argin,N_tot,xGrid); %call integrator...
            clear argin; clear var; clear dim;
            F{odr} = F{odr}  + Matr0;
        end
    end
    writeFunction(F{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], 'Drg_proj1', [])
end


% F_sym = evalMatr_K(F,Q)'*sym(dx);
% F_sym = simplify(F_sym);
% F = matlabFunction(F_sym, 'File','AeroEqns\Drg_proj1','vars', {Q});

%horizontal component...

trms = exprs.DrgTrm{2};

clear F
%assign null matrices...
F{1} = zeros(N_gam,N_tot);
F{2} = zeros(1,N_tot,N_gam,N_tot);
F{3} = zeros(N_tot,N_tot,N_gam,N_tot);

for odr=1:3%min([3,length(trms)])
    if odr<=min([3,length(trms)])
        for numTrms = 1:length(trms{odr})

            %ident scaling fcators and contributing terms...
            [scalFac_sym, Var_sym] = coeffs(trms{odr}(numTrms), varDef);
            scalFcn = matlabFunction(scalFac_sym, 'vars', parDef);
            scalF = @(x)(scalFcn(b(x),a(x)));

            %ident the comtributing shape functionals....
            allocPos = 1;
            for varNum=1:length(varDef)
                pow = polynomialDegree(Var_sym,varDef(varNum));
                if pow>0
                    for powCount=1:pow
                        var{allocPos} = varAlloc{varNum};

                        %record the dimensionality of the assigned vector...
                        if varNum==1, dim{allocPos} = 2; else, dim{allocPos} = 1; end
                        allocPos = allocPos+1;
                    end
                end
            end

            argin.odr = odr; argin.scalF = scalF;
            spaces = zeros(1,odr); % indicator to record free vectors argument spaces..
            for j=1:length(var)
                avail = find(spaces~=1); %identify free spaces..
                if dim{j}==1
                    argin.vctr{avail(1)} = @(x,J)(var{j}(x));
                    spaces(avail(1))=1;
                end
                if dim{j}==2
                    if spaces(1)==0
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(3)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,[avail(1),avail(3)]) = 1; %lock assigmend space
                    else
                        argin.vctr{avail(1)} = @(x,J)(ones(N_tot,1));
                        argin.vctr{avail(2)} = @(x,J)(dimDelim(var{j}(x),J));
                        spaces(1,avail(1:2)) = 1; % lock assigned space
                    end
                end
            end
            Matr0 = matCalc_Aer(argin,N_tot,xGrid); %call integrator...
            clear argin; clear var; clear dim;
            F{odr} = F{odr}  + Matr0;
        end
    end
    writeFunction(F{odr}, odr, [fPath,'\+odr_',num2str(odr),'\'], 'Drg_proj2', [])
end


% F_sym = evalMatr_K(F,Q)'*sym(dx);
% F_sym = simplify(F_sym);
% F = matlabFunction(F_sym, 'File','AeroEqns\Drg_proj2','vars', {Q});


%% downwash equation...

Q_lam = sym('q_lam',[NA,1]); assume(Q_lam, 'real');
f_del = eye(length(xGrid));
for j=1:length(xGrid)-1
    f_del(j,j+1) = -1;
end

% for i=1:N_gam
%     for j=1:N_gam
%         f_ind(i,j) = (0.25/pi)*(1./(xTrail(j)-xGrid(i)) + 1./(xTrail(j)+xGrid(i)));
%     end
% end

%sweep version
if beta~=0
    for i=1:N_gam
        for j=1:N_gam

            swpFct1 = (1-sin(beta))./(cos(beta));

            %sweep - left wing
            dx1 = (xTrail(j)+xGrid(i))*cos(beta);
            dy = (xTrail(j)-xGrid(i))*sin(beta);
            dr = sqrt(dx1^2+dy^2);
            sB = (dy./dr);
            cB = dx1./dr;
            swpFct2 = (1-sB)./cos(beta);

            %2nd effect
            cs1=cB*cos(beta)+sB*sin(beta);

            f_ind(i,j) = (0.25/pi)*(swpFct1./(xTrail(j)-xGrid(i)) + swpFct2./(xTrail(j)+xGrid(i)));

            %...
            %sweep - left wing
            if j>1
                dx2 = (xTrail(j-1)+xGrid(i))*cos(beta);
                dy = (xTrail(j-1)-xGrid(i))*sin(beta);
            else
                dx2 = (0+xGrid(i))*cos(beta);
                dy = (0-xGrid(i))*sin(beta);
            end
            dr = sqrt(dx2^2+dy^2);
            sB = (dy./dr);
            cB = dx2./dr;

            %2nd effect
            cs2=-cB*cos(beta)-sB*sin(beta);
            dp = sin(2*beta)*xGrid(i);

            f_ind2(i,j) = (0.25/pi)*(cs1+cs2)./dp;
        end
    end
else
    f_ind2 = 0;
    for i=1:N_gam
        for j=1:N_gam
            f_ind(i,j) = (0.25/pi)*(1./(xTrail(j)-xGrid(i)) + 1./(xTrail(j)+xGrid(i)));
        end
    end
end

u_ind_matr = f_ind*f_del + f_ind2;

writeFunction(u_ind_matr, 1, ['+project\+aero\+flow\'], 'u_ind', [])
writeFunction(bx, 1, ['+project\+aero\+geom\'], 'bi_fcn', [])
writeFunction(dx, 1, ['+project\+aero\+geom\'], 'dx_fcn', [])

%% function assign coordinates to function...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function outVec = evalMatr_K(matr, q)
        sz = size(matr{1}); N = sz(2);
        for J=1:length(xGrid)
            for I=1:N
                subElem{1}(J,I) = matr{2}(1,:,J,I)*q +...
                    q'*matr{3}(:,:,J,I)*q;
            end
        end
        outVec = (subElem{1} + matr{1});
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% function to output selected columns from pre-integrated matrices...
    function out = dimDelim(mat, dim)
        if length(dim)==1
            out = mat(1:N_tot,dim);
        end
        if length(dim)==3
            out = mat(1:N_tot,dim);
        end
    end
end