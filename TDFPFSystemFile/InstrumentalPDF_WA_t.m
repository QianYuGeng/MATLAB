function [Xfail_dv,fxhx_no,Eu,Betaint,Ifail] = InstrumentalPDF_WA_t(Probs,N,Method,ind_t,Ifail_old)
% WA: weighted approach with intrumental PDF H(x) which is embeded in Probs, WA including WMCS, WIS and WLS
% Xfail_dv : failure samples for x_r 
% fxhx_no  : rato of PDFs for x_n: 
% updated: 2021-08-12 putting WMCS, WIS and WLS in one m-file
% updated: 2021-09-15 WLS:change Xint to Betaint

    % addpath('..\AProblemDefinition');
    ts=0:Probs.dt:Probs.tmax;
    
    if Probs.NE == 2018
        [Xfail_dv, fxhx_no] = BenchProb2_AuIS(Probs,N);
        Eu = 1; Betaint = 1 ;
        return
    end
    
    nsample = 1;  % updated : 2021-08-12
    while nsample <= N   
        if strcmp(Method,'WMCS') == 1  || strcmp(Method,'WMCS_BOC_Ave') == 1  || strcmp(Method,'WMCS_BOC_Var')== 1 || strcmp(Method,'WMCS_BOC_Cov') == 1  % MCS 
            for i = 1:Probs.Nx
                Xpoint(nsample,i) = random(Probs.Dist{i}, Probs.Para{i}(1),Probs.Para{i}(2)); 
            end
        elseif strcmp(Method,'WIS') == 1 || strcmp(Method,'WIS_BOC_Ave') == 1 || strcmp(Method,'WIS_BOC_Var') == 1 || strcmp(Method,'WIS_BOC_Cov') == 1  % IS
            if isfield(Probs,'XPt')
            Hmu = Probs.XPt(ind_t,:);
            else
                Hmu = Probs.XP;
            end
            Hsd = Probs.sd;
            Xpoint(nsample,:) = normrnd(Hmu,Hsd); 
        elseif strcmp(Method,'WLS') == 1 || strcmp(Method,'WLS_BOC_Ave') == 1 || strcmp(Method,'WLS_BOC_Var') == 1 || strcmp(Method,'WLS_BOC_Cov') == 1  % LS
            Hmu = Probs.mu;    Hsd = Probs.sd;
            Xpoint(nsample,:) = normrnd(Hmu,Hsd); 
        end

        if  Probs.PositiveCheck == 1    % Check the minus value
            Ifsum = sum(Xpoint(nsample,:)<0,2);
            if Ifsum==0   % no minus values, receive it
                nsample = nsample+1;
            end
        else
            nsample = nsample+1;
        end
    end

    % fxhx_no: rato of PDFs for x_no: variables for un-related with design parameters 
    if strcmp(Method,'WMCS') == 1  || strcmp(Method,'WMCS_BOC_Ave') == 1  || strcmp(Method,'WMCS_BOC_Ave') == 1 || strcmp(Method,'WMCS_BOC_Cov') == 1  % MCS
        gx = Fung(Xpoint);
        Ifail = gx<0;
        Xfail = Xpoint(Ifail,:);
        fxhx_no = 1;
        Eu = 1; Betaint = 1 ;
    elseif strcmp(Method,'WIS') == 1 || strcmp(Method,'WIS_BOC_Ave') == 1 || strcmp(Method,'WIS_BOC_Var') == 1 || strcmp(Method,'WIS_BOC_Cov') == 1  % IS
        gx = Fung(Xpoint,ts(ind_t),Probs);
        Ifail = gx<0;
        Ifail=or(Ifail,Ifail_old);
        Xfail = Xpoint(Ifail,:);
       [Nfail, ~]= size(Xfail);
        fxhxi = ones(Nfail,Probs.Nx);
        for i = Probs.xuLocal   % Probs.Nd+1 : Probs.Nx  updated: 2021-05-30
            fxi = pdf(Probs.Dist{i},Xfail(:,i), Probs.Para{i}(1),Probs.Para{i}(2));
%             hxi = normpdf(Xfail(:,i), Probs.XP(i),Probs.Hsd(i));
            hxi = normpdf(Xfail(:,i), Hmu(i),Hsd(i));
            fxhxi(:,i) = fxi./hxi;
        end
        fxhx_no = prod(fxhxi,2);
        Eu = 1; Betaint = 1 ;
        
    elseif strcmp(Method,'WLS') == 1 || strcmp(Method,'WLS_BOC_Ave') == 1 || strcmp(Method,'WLS_BOC_Var') == 1 || strcmp(Method,'WLS_BOC_Cov') == 1  % LS   
    end
            
    % return
    Xfail_dv = Xfail(:,Probs.dvLocal);
    Betaint ;
    fxhx_no;
end    
 function gx=Fung(x,tk,Probs) %function gx = Fungxt(x,tk,Probk)
%     output=Probk.Fungt(x,t,F);
    if ~isfield(Probs,'Fungt')
    gx=Fung(x);
    else
    
    ts  = 0: Probs.dt : Probs.tmax; 
    IF = find(ts==tk);  % 指定时间处的LSF
    
    X =  x(:,1:Probs.Nr); Z = [] ; Ftk = []; gx = [];
    [N,~] =  size(x); 
    for j = 1:N
        if isfield(Probs,'NF') && Probs.NF>0     
            Z =  x(j,Probs.Nr+1:end); 
            Ft{j} = Probs.ZpToFt(Z,Probs);       % Ft: [NF,Nt]
            Ftk = Ft{j}(:,IF);%时间在列Ftk [:,Nt]
        end
        gx(j,:) = Probs.Fungt(X(j,:),tk,Ftk) ;

    end
    end
 end
        

    