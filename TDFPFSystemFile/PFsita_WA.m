function [Pfinequ, Pfequ, gradc, gradceq]= PFsita_WA(Sita,Xfail,fxhx_no,Eu,Betaint,Probs,N,Method,flag)
% 加权抽样的失效概率函数估计-即RBDO概率约束,WMS WIS and WLS
% Smu  : input design varibles value
% Xfail: fail samples, all x or x_dv is both OK, [Nfail, Nx or Nv]
% fxhx_no : the part that not relevent with the desgin varibales, should contain all of them
% Probs: the setting associated with the instrumental PDF
% flag : 'Pf&Cov' for output the cov, otherwise for the fmincon fuction
% updated : 2021-08-12, put WIS and WLS together
% updated : 2021-08-13, Changed from Cons_PFsita_WA,fix it to Problem_FPF.m
% updated : 2021-09-18, Changed WLS part according to Marcos

    global Prob indt
    Idv = Prob.dvLocal;
    % if nargin < 6,  flag = 'Pfequ'; end
    
    Hmu = Probs.mu;  Hsd = Probs.sd; % instrumental sampling PDF which is embedded in Probs. 
    
    [m,ndv] = size(Sita);
    if ndv ~= Prob.Nd     Sita = Sita(:,Idv);   end
    
    [Nfail, nxf] = size(Xfail);
    if nxf ~= Prob.Nd            
        Xfail_dv = Xfail(:,Idv); 
    else
        Xfail_dv = Xfail;
    end

    %compute the IS weight (f_x_r/h_x_r) of X_r, variable corresponding to design paremeter
    Probk = rmfield(Probs,'mu');  % Probk is changed as each Sita value 2021-09-21
    for k = 1:m  % for each different smu = sita
        % 2021-08-13 put Sita into Probs 
        fv=zeros(Nfail,Probk.nxr);
        for i = 1:Probk.nxr 
            si_Para = Probk.si_Para{i};        % index of s_i in Para, [1] or [1, 2]
            si_Sita = Probk.si_Sita{i};        % index of s_i in Sita, [3] or [3, 4]
            ixr = Probk.dvLocal(i);              % index of x_r in x, i.e., [1 3 5]
            Probk.Para{ixr}(:,si_Para) = Sita(k,si_Sita);  % put s_i into ParaN = {[N,nx], [N,nx]}   
            
            ParaN{ixr} = Probk.Para{ixr};      
            ParaN{ixr}(:,si_Para) = Sita(k,si_Sita); % put s_i into ParaN = {[N,nx], [N,nx]}       
            fv(:,i) = pdf(Probk.Dist{ixr},Xfail_dv(:,i),ParaN{ixr}(:,1),ParaN{ixr}(:,2));  %2021-08-12
        end
        Probk = PrepareProb(Probk); % 2021-09-21
        
        if    strcmp(Method,'WMCS') || strcmp(Method,'WMCS_BOC_Ave') || strcmp(Method,'WMCS_BOC_Var') || strcmp(Method,'WMCS_BOC_Cov') 
            for i = 1:Probk.nxr   % changed from Probs.Nd  2021-8-13
                 hv(:,i) = pdf(Probk.Dist{i},Xfail_dv(:,i),repmat(Hmu(:,i),Nfail,1),repmat(Hsd(i),Nfail,1));
            end
            fxhx_no  =  1;
        elseif strcmp(Method,'WIS') || strcmp(Method,'WIS_BOC_Ave')  || strcmp(Method,'WIS_BOC_Var')  || strcmp(Method,'WIS_BOC_Cov')
            if ~isfield(Prob,'XPt')
            Hmu = Probk.XP;  % instrumental sampling PDF
            else
             Hmu = Probk.XPt(indt,:);
            end
            hv=zeros(Nfail,Probk.nxr);
           for i = 1:Probk.nxr   % Note that: only for two parameters case
                 hv(:,i) = normpdf(Xfail_dv(:,i),repmat(Hmu(i),Nfail,1),repmat(Hsd(i),Nfail,1));
           end
        elseif strcmp(Method,'WLS') || strcmp(Method,'WLS_BOC_Ave')  || strcmp(Method,'WLS_BOC_Var')  || strcmp(Method,'WLS_BOC_Cov')   % LS
            % Ref. Song .etl.. SS, 84(2020) 101936
            for i = 1:Probk.nxr % Note that: only for two parameters case
                 hv(:,i) = fv(:,i); % pdf(Probs.Dist{i},Xfail_dv(:,i),repmat(Hmu(:,i),Nfail,1),repmat(Hsd(i),Nfail,1));
            end
            % 计算积分项 2021-09-15
            fmu = Probk.mu;  fsd = Probk.sd;  % original pdf
            % uSita_int  =  (Xint - repmat(Hmu,N,1))./repmat(Hsd,N,1);  % u value in the H(x|sita) space cooresponding to Xint
            ks = Hsd(:,Idv)./fsd(:,Idv);                    
            epsino = prod(ks,2);                                % Eq. A2    
            zeta  = sum((1- ks.^2).*Eu(:,Idv).^2/2 , 'all');    % Eq. A4
           for j = 1:N
                uj = (Xfail_dv(j,:) - Hmu(:,Idv))./Hsd(:,Idv); % Standard normal samples u, 
                up = uj - Eu(:,Idv)*uj'.*Eu(:,Idv) ;           % Vector perpendicular u^p to importance direction Eu 
                
                % ku = ( Sita(k,:)-Hmu(:,Idv))./Hsd(:,Idv);      % ki: offset 偏差 if only mu is changed, Marcos
                % namda(j) = sum( ku.*up - ku.^2/2, 'all');      % namda_j 指数项 if only mu is changed, Marcos
                % K(j)     = sum(ku.*Eu(:,Idv),'all');           % 
                
                ku = ( Hmu(:,Idv) + Hsd(:,Idv).*up(:,Idv) - fmu(:,Idv))./fsd(:,Idv);      % ki: offset 偏差 
                namda(j) = sum(up.^2/2 - ku.^2/2, 'all');       % namda_j 指数项 Eq. A4
                K(j)     = sum(up.*Eu(:,Idv) - ku.*ks.*Eu(:,Idv),'all');  % k

                betaj(j)  = Betaint(j,:); 
                %Inte(j,:) = exp( namda(j)+ K(j).^2/2  ) * normcdf( K(j)-betaj(j) ) ;  % if only mu is changed, Marcos
                Inte(j,:) = epsino/sqrt(1-2*zeta) * exp( namda(j) + K(j).^2/(2-4*zeta)  ) * normcdf( (K(j)-(1-2*zeta)*betaj(j))/sqrt(1-2*zeta)  ) ;  % Eq. A9
            end
           fxhx_no  = Inte;%  normcdf(-BetaSita);   %Pfj =  normcdf(-BetaSita);
        else
            error('No such method');
        end
        
        fvhv = prod(fv./hv,2);
        
        % 失效概率函数 及 置信带
        Pf(k)    = 1/N *sum(fxhx_no.*fvhv); 
        VarPf(k) = 1/(N-1)*(1/N*sum((fxhx_no.*fvhv).^2,'all')-Pf(k).^2);   % it is same with the following 2021-05-01
        % VarPf(k) = 1/(N*(N-1))*sum((Pfj.*fvhv-Pf(k)).^2,'all');
        CovPf(k) = sqrt(VarPf(k))/Pf(k);
    end
        
    % return constraint value
    Pfinequ = Pf - Probk.Pftol;    
    if isfield(Probk,'dcons')  % if there are determinstic constraints, included here    
        Pfinequ =[Pfinequ; Probk.dcons1(Sita)]; 
    end
    Pfequ = [];
    
    % return the real Pf and Pfcov, for some cases used
    if strcmp (flag, 'Pf&Cov')
        Pfinequ = Pf;
        Pfequ   = CovPf;
    end
    
    % output the grad if necessary
%     if nargout > 2
%         gradc = Pfdmu';  % in 'fmincon', vertical row are accepted.
%         gradceq = [];
%         fprintf('\nCurrent Sita: ');fprintf('%6.4f\t',Smu); fprintf('\n');
%         fprintf('Grad_Pf(sita): [');fprintf('%-12.4e\t ',gradc);fprintf(']\n');
%     end
