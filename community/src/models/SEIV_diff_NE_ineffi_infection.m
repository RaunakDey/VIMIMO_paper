classdef SEIV_diff_NE_ineffi_infection < ode_funs
    
    properties
        name;
        statevars = ["S","E","I","V"];
        NH;
        NV;
        NE;
        id;
    end
    
    methods
        
        function obj = SEIV_diff_NE_ineffi_infection(NH,NV,NE)
            obj.name = sprintf('SEIV%d',NE);
            obj.NH = NH;
            obj.NV = NV;
            obj.NE = NE;
            id.S = 1:NH;
            id.E = (1:NH*NV*NE)+id.S(end);
            id.I = (1:NH*NV)+id.E(end);
            id.V_infec = (1:NV)+id.I(end);
            
            id.V_noninfec = (1:NV) + id.V_infec(end);

            id.D = 1+id.V_noninfec(end); % debris
            id.Emat = reshape(id.E,[NH NV NE]);
            id.Imat = reshape(id.I,[NH NV]);
            obj.id = id;
        end

        function vec = zeros(obj)
            %vec = zeros(obj.id.V(end),1);
            vec = zeros(obj.id.D,1);
        end
        
        function S = sum_hosts(obj,y) % all host cells including infected
            S = y(:,obj.id.S);
            for i = 1:obj.NH % could use bsxfun or arrayfun here?
               S(:,i) = S(:,i) + sum(y(:,obj.id.Emat(i,:)),2) + sum(y(:,obj.id.Imat(i,:)),2); 
            end
        end
        
        function V = sum_viruses(obj,y) % free phage only
           V = y(:,obj.id.V_infec) + y(:,obj.id.V_noninfec);
        end
        
        function dydt = ode(obj,t,y,pars)
            
            OH = ones(obj.NH,1);
            OV = ones(obj.NV,1);
            S = y(obj.id.S);
            Emat = y(obj.id.Emat);
            Imat = y(obj.id.Imat);
            V_infec = y(obj.id.V_infec);
            V_noninfec = y(obj.id.V_noninfec);
            D = y(obj.id.D);
            N = S+sum(Emat,3)*OV+Imat*OV;
            etaeff = pars.eta.*(pars.NE+1);

            exposed_transition_fun = obj.exposed_transition_fun;
            host_growth_fun = obj.host_growth_fun;
            viral_decay_fun = obj.viral_decay_fun;
            viral_adsorb_fun = obj.viral_adsorb_fun; % only use in dV
            lysis_reset_fun = obj.lysis_reset_fun;
            debris_inhib_fun = obj.debris_inhib_fun;
                  
            Sdeb = S*debris_inhib_fun(pars,D); % debris
            
            dS = host_growth_fun(pars,S,N) - Sdeb.*((pars.M.*pars.phi)*V_infec);
            dEmat = (pars.M.*pars.phi).*(Sdeb*V_infec') - etaeff.*Emat(:,:,1);
            dEmat2 = exposed_transition_fun(etaeff,Emat);
            for i = 1:obj.NH
                for j = 1:obj.NV
                    if (pars.NE(i,j) ~= obj.NE  && pars.NE(i,j) ~=0) 
                        Emat(i,j,pars.NE(i,j)+1:end) = Emat(i,j,pars.NE(i,j));
                    end
                end
            end
            
            dImat = etaeff.*Emat(:,:,end) - etaeff.*Imat ;
       
            dV_infec = (pars.prob_infec.*pars.beta.*etaeff.*Imat)'*OH - V_infec.*(viral_adsorb_fun(pars)'*N);
            dV_noninfec = ((1-pars.prob_infec).*pars.beta.*etaeff.*Imat)'*OH - V_noninfec.*(viral_adsorb_fun(pars)'*N);

            dD = sum(etaeff(:).*Imat(:)); % sum across all pairs for net lysis rate
            
            dydt = [dS; dEmat(:); dEmat2(:); dImat(:); dV_infec; dV_noninfec; dD];
            
        end
    end
end

