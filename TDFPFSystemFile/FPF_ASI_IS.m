function [Pfsita, CovPfsita] = FPF_ASI_IS(Sita,N,Xfail,fzhz,Probs,MethodC)
% FPF by ASI-MCS, importance sampling in augmented space
% updated: 2020-12-23
% updated: 2021-08-07

    dvLocal = Probs.dvLocal; Para = Probs.Para;  
    XP = Probs.XP;  lb = Probs.lb;  ub = Probs.ub; 
    
    [Nf, nxf] = size(Xfail);
    if nxf ~= Probs.Nd    
        Xfail_dv = Xfail(:,dvLocal); 
    else
        Xfail_dv = Xfail; 
    end

    % Function that in denominator
    if strcmp(MethodC,'WIS') == 1
        deta_X = prod(ub-lb)*prod(normpdf(Xfail_dv,repmat(XP(:,dvLocal),Nf,1),  repmat(Probs.sd(:,dvLocal),Nf,1)),2);% WA is f_x_S0
    else
        deta_X = Deta_X(Xfail,Probs);   % [Nf,1], not related with Sita
    end

    [Ns, ~] = size(Sita);
    for p = 1:Ns  % for each different sita[1*Nd]
        % f(x_r| s_i)
        f_xr_sita = fxr_sita(Xfail,Sita(p,:),Probs);

        % H(x_u|x_u*)
        if ~isempty(Probs.xuLocal) 
            for i = Probs.xuLocal    
                 f_xuj(:,i) = pdf(Probs.Dist{i}, Xfail(:,i), Para{i}(:,1), Para{i}(:,2));
                 H_xuj(:,i) = normpdf( Xfail(:,i),XP(i),Probs.sd(i));  % usually norm ISD is used!!
            end
            f_xu    = prod(f_xuj(:,Probs.xuLocal),2);        
            H_xu_xp = prod(H_xuj(:,Probs.xuLocal),2);      
            fxuhxu  = f_xu./H_xu_xp;
        else
            fxuhxu = 1;
        end
        % FPF estimate
        Ratio = prod(ub-lb) .* ( f_xr_sita./deta_X )  .*  fxuhxu .* fzhz;  % 2021 updated
        Pfsita(p,:) = 1/N*sum(Ratio);

        % cov
        VarPfsita = 1/(N-1)* (1/N* sum((Ratio).^2)  - Pfsita(p,:).^2  );
        CovPfsita(p,:) = sqrt(VarPfsita)/Pfsita(p,:);
    end
    
    
    
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%    
function f_xr_sita = fxr_sita(Xfail,Sita_p,Probs)
    % compute the PDF of x corresponding to sita    
        % prepare the para for x_r
        dvLocal = Probs.dvLocal;  Para = Probs.Para;
        
        [Nf, ~] = size(Xfail);
        for i = 1:Probs.nxr 
            ixr = dvLocal(i);              % index of x_r in x
            i_sinPara = Probs.si_Para{i};        % index of s_i in Para 1 or 2
            i_sinSita = Probs.si_Sita{i};       
            ParaN{ixr} = repmat(Para{ixr},Nf,1);         
            ParaN{ixr}(:,i_sinPara) = repmat(Sita_p(:,i_sinSita),Nf,1);         
        end
        
        % f(x_r| s_i)
        fxrj_sita = ones(size(Xfail));
        for i = dvLocal  % Note: not (1:Probs.nxr), as it is not from 1 
            fxj_sita(:,i) = pdf(Probs.Dist{i}, Xfail(:,i), ParaN{i}(:,1), ParaN{i}(:,2));  % f(x^(j) | sita)
        end
        f_xr_sita = prod(fxj_sita(:,dvLocal),2);
    
    
function deta_X = Deta_X(Xfail,Probs)
% calculate Deta(X) for each X sample, based on certain Sita 
% updated: 2020-11-28, estended to any design parameter, any location
%  Xfail : failure samples
%  Sita  : those not intergral at certain value, Sita0, [1, ns]
%  lb,ub : The design domain of sita
% 
   
    [Nf, ~] = size(Xfail);
    for j = 1:Nf
        for i = 1:Probs.nxr  % number of xr
            ixr = Probs.dvLocal(i);      % index of x_r in x
            isP = Probs.si_Para{i};      % index of s_i in Para 1 or 2
            Domain_is = Probs.Domain_x_s{i};
            [~,n_sofx] = size(isP);
            if n_sofx == 1  % one parameter for single x
                D_s1 = Domain_is{1};
                % only for single sita, note there may be two design parameters (sita) for single x, so sita should be all inputed 
                % deta_xi(i) = quadgk(@(s) fxi_sita(s,Sita,Xfail(j,ixr),i),D_s1(1),D_s1(2));     
                 deta_xi(i) = quadgk(@(s) fxi_sita_1d(s,Xfail(j,ixr),i,Probs),D_s1(1),D_s1(2));     
                
                % deta_xi(i) = integral(fxtheta,sitaL,sitaU)
            elseif n_sofx == 2 % two parameter for single x changing at the same time, it needs two-dimension integral
                D_s1 = Domain_is{1};D
                D_s2 = Domain_is{2};
                deta_xi(i) = quad2d(@(s1,s2) fxi_sita_2d(s1,s2,Xfail(j,ixr),i),D_s1(1),D_s1(2),D_s2(1),D_s2(2),Probs);     

            else  % more than two (three) parameters for single x
                error('too many domain for sita');
            end
            
            % produce all the xr_i together
            deta_X(j,:) = prod(deta_xi,2);

        end
    end
    
function fxi = fxi_sita_2d(s1,s2,X_i,ith_xr,Probs)
% For integral to calculte Deta_x, case: two design s in one x_i
% s_i : the one of s for integral,  note that the imput of s_i is a line vector for intergral calculation.
% Sita: all the sita [1, N_s], as there maybe two s in single x 
% X_i : Xr(i) related with s 
% index_s: index of s
% updated: 2020-12-06

        Para = Probs.Para;  
        
        [n1_s1,n2_s1] = size(s1); % s_i is a line vector for intergral calculation
        [n1_s2,n2_s2] = size(s2);
        if n1_s1~=n1_s2  % checking 
            fprintf('n1_s1~=n1_s2'); error('s1 and s2 not matched ')
        elseif n2_s1~=n2_s2
            fprintf('n2_s1~=n2_s2'); error('s1 and s2 not matched ')
        end

        index_xr = Probs.dvLocal(ith_xr);
        index_s  = Probs.si_Sita{ith_xr};       
        for i = 1:n1_s1 
            for j = 1:n2_s1 
                Para{index_xr} = [s1(i,j),s2(i,j)];
                fxi(i,j) = pdf(Probs.Dist{index_xr},X_i,Para{index_xr}(:,1),Para{index_xr}(:,2));    % f(x|sita) for integral   
            end
        end
        
function fxi = fxi_sita_1d(s1,X_i,ith_xr, Probs)
% For integral to calculte Deta_x, case: one design s in one x_i
% s_i : the one of s for integral,  note that the imput of s_i is a line vector for intergral calculation.
% Sita: all the sita [1, N_s], as there maybe two s in single x 
% X_i : Xr(i) related with s 
% index_s: index of s
% updated: 2020-12-06

        Para = Probs.Para;   

        [n1_s1,n2_s1] = size(s1);
        if n1_s1~=1  % checking 
            fprintf('n1_s1~=1'); error('s1 is not a line vector ')
        end

        index_xr = Probs.dvLocal(ith_xr);
        index_s  = Probs.si_Sita{ith_xr};       
        is  = Probs.si_Para{ith_xr};      % index of s_i
        for i = 1:n2_s1 
            Para{index_xr}(is(1)) = s1(:,i);
            fxi(:,i) = pdf(Probs.Dist{index_xr},X_i,Para{index_xr}(:,1),Para{index_xr}(:,2));    % f(x|sita) for integral   
        end

