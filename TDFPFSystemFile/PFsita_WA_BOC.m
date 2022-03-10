function [Pfinequ, Pfequ] = PFsita_WA_BOC(Sita,Xfail,fxhx_no,Eu,Betaint,Probk,N,Method,flag)
% Reused algorithm, through the Best optimal Combination (according to minimization of Cov) based on Weighted approach
% updated: 2021-04-29
% updated: 2021-05-29, change to Problem_FPF problem which suit for all dv
% Xfail: cell
% hx   : cell
% updated: 2021-08-13: put FPF_ASI_WMCS PFsita_WIS and FPF_ASI_WLS together
% updated: 2021-09-13: Betaint

    global Prob 
    [~,nks] = size(Xfail);
    
    for ks = 1:nks  % 2021-08-13 all in one file WA
        [PFs(:,ks), CovPFs(:,ks)] = PFsita_WA(Sita,Xfail{ks},fxhx_no{ks},Eu{ks},Betaint{ks},Probk{ks},N,Method,'Pf&Cov');
    end
    CovPFs(isnan(CovPFs))=1;
    % Best optimite combining: using the weights that minimizes the Cov(PF(sita)) 
    ncomb = 1; % Prob.ncomb;% 1;  %  
    if nks> ncomb
        if     strcmp(Method,'WLS_BOC_Ave') == 1 || strcmp(Method,'WIS_BOC_Ave') == 1  || strcmp(Method,'WMCS_BOC_Ave') == 1    
             wk =  1/(nks-ncomb+1);
            % EcvPFs = PFs.*CovPFs.^2;
            % wk =  EcvPFs(:,ncomb:nks).^(-1)./sum(EcvPFs(:,ncomb:nks).^(-1),2);
        elseif strcmp(Method,'WLS_BOC_Var') == 1 || strcmp(Method,'WIS_BOC_Var') == 1 || strcmp(Method,'WMCS_BOC_Var') == 1
            VarPFs = (PFs.*CovPFs).^2;
            wk =  VarPFs(:,ncomb:nks).^(-1)./sum(VarPFs(:,ncomb:nks).^(-1),2);
        elseif strcmp(Method,'WLS_BOC_Cov') == 1 ||  strcmp(Method,'WIS_BOC_Cov') == 1 || strcmp(Method,'WMCS_BOC_Cov') == 1
            wk =  CovPFs(:,ncomb:nks).^(-2)./sum(CovPFs(:,ncomb:nks).^(-2),2);
        else
            error('Method selection is wrong! ');
        end
        
        PFs_BOC = sum(wk.*PFs(:,ncomb:nks),2);
        VarPFs_BOC = sum(wk.^2.*(CovPFs(:,ncomb:nks).*PFs(:,ncomb:nks)).^2,2);
        CovPFs_BOC = sqrt(VarPFs_BOC)./PFs_BOC;
    else
        PFs_BOC = PFs(:,1);
        CovPFs_BOC = CovPFs(:,1);
    end
    
    % for outloop of optimization
    Pfinequ = PFs_BOC - Prob.Pftol;
    if isfield(Prob,'dcons')    
        Pfinequ =[Pfinequ; Prob.dcons(Sita)']; 
    end
    Pfequ = 0;
    
    % return the real Pf and Pfcov
    if strcmp (flag, 'Pf&Cov')
        Pfinequ = PFs_BOC;
        Pfequ   = CovPFs_BOC;
    end        