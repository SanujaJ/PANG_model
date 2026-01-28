classdef analysisBase

    properties
        transF
        par
        par0
        p = [];
        p0 = [20, 0, 0, 9.81, 1.225, 0];
        dampMatr

        dispFlag = true;
    end

    properties (SetAccess = {?buildSystem.buildBase})
        basis
        geom
        file
    end

    properties (Dependent)
        Nstr
        structDisp
        structVel
        q0_struct
    end

    methods

        function obj = analysisBase(obj_in)
            arguments
                obj_in buildSystem.buildBase
            end
            obj.basis = obj_in.basis;
            obj.geom = obj_in.geom;
            obj.transF = eye(obj_in.basis.Ntot,obj_in.basis.Ntot);
            obj.par = obj_in.par;
            obj.par0 = obj_in.par0;
            obj.dampMatr = zeros(obj_in.basis.Ntot);
            obj.file = obj_in.file;
        end

        %%
%                 function obj = set.dampMatr(obj,val)
%                     sz = size(val);
%                     if sz(1)~=obj.basis.Ntot
%                         error(['Incorrect size for damping matrix, must be: ', num2str(obj.basis.Ntot), ' by ', num2str(obj.basis.Ntot)]);
%                     else
%                         if sz(2)~=obj.basis.Ntot
%                             error(['Incorrect size for damping matrix, must be: ', num2str(obj.basis.Ntot), ' by ', num2str(obj.basis.Ntot)]);
%                         end
%                     end
%                     obj.dampMatr = val;
%                 end

        function q0_struct = get.q0_struct(obj)
            q0_struct = zeros(2*obj.Nstr,1);
        end

        function Nstr = get.Nstr(obj)
            Nstr = length(obj.transF(1,:));
        end

        function structDisp = get.structDisp(obj)
            Nstr = length(obj.transF(1,:));
            structDisp = [1:Nstr];
        end
        function structVel = get.structVel(obj)
            Nstr = length(obj.transF(1,:));
            structVel = [Nstr+1:2*Nstr];
        end

        %%
        function obj = setTransform(obj, type, DoF)

            if nargin<3 %if not specified, use the complete set of DoFs
                DoF = obj.basis.Ntot;
            end

            if DoF>obj.basis.Ntot
                error('Requested degrees of freedom exceeds model size');
            end

            if strcmp(type, 'default')
                if DoF>obj.basis.Ntot
                    DoF = obj.basis.Ntot;
                    warning(['default transformation requires the entire model size...re setting to ',num2str(obj.basis.Ntot)])
                end
            end

            %first, re-set to default..
            obj.transF = eye(obj.basis.Ntot);

            %now change..
            switch type
                case 'default'
                    obj.transF = [eye(DoF); zeros(obj.basis.Ntot - DoF, DoF)];
                case 'modal'
                    obj_temp = obj;
                    obj_temp.dampMatr = zeros(obj.basis.Ntot);
                    [shp, evals] = getStructModes(obj_temp,zeros(2*obj.basis.Ntot,1));
                    for j=1:DoF
                        mc_r=real(shp(:,j));
                        mc_i=imag(shp(:,j));
                        [uu,ss,vv]=svd([mc_r,mc_i]');
                        dq=uu(:,1).'*[mc_r.';mc_i.'];   % real mode
                        shpMat(:,j) = dq';
                    end
                    obj.transF = shpMat;
            end
        end

        %% function for eigenvalue analysis - structural..
        function [shp, evals, K, C, M] = getStructModes(obj,q0,varargin)

            [Mglob_0,Fglob_0] = obj.structModel(q0, varargin{:});
            epsilon = 1e-6;

            K = zeros(obj.Nstr, obj.Nstr);
            M = Mglob_0(obj.Nstr+1:end, obj.Nstr+1:end);
            C = obj.transF'*obj.dampMatr*obj.transF;

            for j=1:obj.Nstr
                q=q0; q(j) = q(j)+epsilon;
                [~,Fglob] = obj.structModel(q,varargin{:});
                K(:,j) = -(Fglob(obj.Nstr+1:end) - Fglob_0(obj.Nstr+1:end))/epsilon;
            end

            [V,D] = polyeig(K, C, M);

            [~,idx] = sort(abs(D));
            for vecNum = 1:length(V(1,:))
                %normalise
                normFac = sqrt(V(:,vecNum)'*M*V(:,vecNum));
                V(:,vecNum) = V(:,vecNum)./normFac; %mass normalised vectors
            end
            evals = D(idx); shp = V(:,idx);

            evals = evals(2:2:end); shp = shp(:,2:2:end);
        end

        function [shp, evals, K, C, M] = plotStructModes(obj,q0,varargin)

            %call the eigenvalue solver...
            [shp, evals, K, C, M] = obj.getStructModes(q0,varargin{:});

            %get 'real' modes.. 
            for j=1:length(evals)
                mc_r=real(shp(:,j));
                mc_i=imag(shp(:,j));
                [uu,ss,vv]=svd([mc_r,mc_i]');
                dq=uu(:,1).'*[mc_r.';mc_i.'];   % real mode
                shpMat(:,j) = dq';
            end

            %loop and plot..
            q0 = q0(obj.structDisp);
            figure;
            for j=1:min([12,length(evals)])
                vctr = shpMat(:,j);

                %get max displacement for scaling..
                [x0, y0, z0] = obj.getMesh(q0, 'aircraft');
                [x, y, z] = obj.getMesh(q0+vctr, 'aircraft');

                dx = max(max(abs(x-x0))); dz = max(max(abs(z-z0)));
                vctr = vctr*0.05*obj.geom.L/max([dz,dx]);
                [x, y, z] = obj.getMesh(q0+vctr, 'aircraft');

                subplot(3,4,j);
                plot3([x0(1,:), flip(x0(2,:))], [y0(1,:), flip(y0(2,:))], [z0(1,:), flip(z0(2,:))], 'k-', 'linewidth', 1); hold on;
                plot3([x(1,:), flip(x(2,:))], [y(1,:), flip(y(2,:))], [z(1,:), flip(z(2,:))], 'r-', 'linewidth', 2); hold on;
                set(gca , 'dataaspectratio', [1,1,1]);
                title(['\omega = ', num2str(abs(evals(j)/(2*pi)), 3), 'Hz, \zeta = ',...
                    num2str(-real(evals(j))./abs(evals(j)),3)]);
            end
        end



        %%
        function obj = setPars(obj,varargin)

            %reference parametert values
            p = obj.p;
            p0 =  obj.p0;

            %update with user provided values..
            if isempty(varargin)
            else
                for p_inpt=1:length(varargin)/2
                    idx0 = find(strcmp(varargin{2*p_inpt-1}, obj.par0));
                    if isempty(idx0)
                        idx = find(strcmp(varargin{2*p_inpt-1}, obj.par));
                        if isempty(idx)
                            error(['Undefined variable ',varargin{2*p_inpt-1}]);
                        else
                            p(idx) = varargin{2*p_inpt};
                        end
                    else
                        p0(idx0) = varargin{2*p_inpt};
                    end
                end
            end

            obj.p = p;
            obj.p0 = p0;

        end

        %%

        function qdt = structDeriv(obj,q_all,varargin)
            [M_glob,F_glob] = obj.structModel(q_all,varargin{:});
            qdt = M_glob\F_glob;
        end

        function [M_glob,F_glob] = structModel(obj,q_all,varargin)

            %reference parametert values
            p = obj.p;
            p0 = obj.p0;

            %update with user provided values..
            if isempty(varargin)
            else
                for p_inpt=1:length(varargin)/2
                    idx0 = find(strcmp(varargin{2*p_inpt-1}, obj.par0));
                    if isempty(idx0)
                        idx = find(strcmp(varargin{2*p_inpt-1}, obj.par));
                        if isempty(idx)
                            error(['Undefined variable ',varargin{2*p_inpt-1}]);
                        else
                            p(idx) = varargin{2*p_inpt};
                        end
                    else
                        p0(idx0) = varargin{2*p_inpt};
                    end
                end
            end

            %assign memory
            M_glob = zeros(2*obj.Nstr, 2*obj.Nstr);
            F_glob = zeros(2*obj.Nstr, 1);
            [M_glob,F_glob] = obj.buildStruct(M_glob,F_glob,q_all,p0,p);

        end

        %% function to get displacements..
        function [x, y, z] = getDisplField(obj,q_all,xStat,frame)
            
            [x, y, z] = analysis.plotFcns.wingDefl(obj, q_all, xStat);

            availTypes = {'aircraft', 'beamModel'};
            if isempty(find(strcmp(frame, availTypes)))
                error(['unexpected frame type ',frame,', Expected: aircraft or beamModel']);
            else
                switch frame
                    case 'beamModel'
                    case 'aircraft'
                        R(1,:) = [y(1,:), y(2,:)];
                        R(2,:) = [x(1,:), x(2,:)];
                        R(3,:) = [-z(1,:), -z(2,:)];

                        beta = obj.geom.sweep;
                        R = [cos(beta), -sin(beta), 0;...
                            sin(beta), cos(beta), 0;...
                            0, 0, 1]*R;

                        x = [R(1,1:end/2); R(1,end/2+1:end)];
                        y = [R(2,1:end/2); R(2,end/2+1:end)];
                        z = [R(3,1:end/2); R(3,end/2+1:end)];
                end
            end
        end
    end

    methods (Static)
        function obj = loadobj(obj)
            addpath([obj.file], '-begin');
            clear persistent
            if obj.dispFlag
                clc;
                fprintf(['...Accessing matrices found in %s\n'], obj.file);
            end
        end
    end
end

