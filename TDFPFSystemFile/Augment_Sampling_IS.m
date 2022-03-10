function [Xfail,fzhz,PF,Xall,Sall] = Augment_Sampling_IS(Probs,N,MethodC)
% IS sampling in augmented space
% update 2020-12-23
% update 2021-08-07, fitting to Problem_FPF definition

    methodinfo = 'IS method for Augmented problem';
    % NE =123;   Probs = OPTModel2020(NE); N = 1000;  lb = Probs.lb; ub = Probs.ub; XP = Probs.XP;  Methods = 1;
    Para = Probs.Para;
    Case = 'CaseB';
    
    if ~isfield(Probs,'XP')
        Probs.XP = PF_AFORM_Norm_Main(Probs);
    end
    XP = Probs.XP;
    
    if     strcmp(Case,'CaseA') == 1 
        Probs.S0 = Probs.lb;    % \theta^(0), initial design
    elseif strcmp(Case,'CaseB') == 1 
        Probs.S0 = (Probs.lb+Probs.ub)/2 ; 
    elseif strcmp(Case,'CaseC') == 1 
        Probs.S0 = Probs.ub; 
    end
    
    %% generate samples from s-U[lb,ub], design parameters£¬Sita
    if strcmp(MethodC,'WIS') == 1   % WIS
        Si = repmat(Probs.S0,N,1);
        for i = Probs.dvLocal    
             x(:,i) = normrnd(XP(i),Probs.sd(i),[N,1]);  % usually norm ISD is used!!
        end
    else
        % ASI-IS
        for i=1: Probs.Nd
            Si(:,i) = unifrnd(Probs.lb(:,i),Probs.ub(:,i),N,1);
        end
        % generate x_r by f(x|s)f(s)
        % prepare the para for sampling 
        for i = 1:Probs.nxr 
            ixr = Probs.dvLocal(i);              % index of x_r in x
            i_sinPara = Probs.si_Para{i};        % index of s_i in Para, [1, 2]
            i_sinSita = Probs.si_Sita{i};        % index of s_i in Sita, [3]
            ParaN{ixr} = repmat(Para{ixr},N,1);         
            ParaN{ixr}(:,i_sinPara) = Si(:,i_sinSita);  % put s_i into ParaN = {[N,nx], [N,nx]}       
        end
        % x_r: generate x_r samples from f(x|s)f(s), conditional distribution, related variables
        for i = Probs.dvLocal
            [~,ns] = size(Para{i});
            if ns==1            x(:,i) = random(Probs.Dist{i}, ParaN{i}(:,1)); end % Note: just suitable for the first parameter!!
            if ns==2            x(:,i) = random(Probs.Dist{i}, ParaN{i}(:,1), Para{i}(:,2)); end
        end
    end

    %% x_u: generate samples from H(z|XP), those not related with s 
    for i = Probs.xuLocal    
         x(:,i) = normrnd(XP(i),Probs.sd(i),[N,1]);  % usually norm ISD is used!!
    end
 
    % evalutation
    if Probs.NE ==2020
        Sfail = Si;
        [Xfail, fzhz, PF] =  BenchProb2_AuIS(Probs,N);
    else
        gx = Probs.Fung(x);
        % Count the number of failure points
        XFI = find(gx<=0);
        Nf = length(XFI);
        Xfail = x(XFI,:);   
        Sfail = Si(XFI,:);
        % The augment failure Probsability
        PF = double(Nf)/(N);
        fzhz = 1;
   end
    
    Xall = x;
    Sall = Si;
    fprintf('\nThe augment failre Probsability(IS) is : \t PF=%f\n', PF);
    
    