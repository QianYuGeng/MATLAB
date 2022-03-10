function [Xfail_dv,fxhx_no,Eu,Betaint] = InstrumentalPDF_WA(Probs,N,Method)
% WA: weighted approach with intrumental PDF H(x) which is embeded in Probs, WA including WMCS, WIS and WLS
% Xfail_dv : failure samples for x_r 
% fxhx_no  : rato of PDFs for x_n: 
% updated: 2021-08-12 putting WMCS, WIS and WLS in one m-file
% updated: 2021-09-15 WLS:change Xint to Betaint

    % addpath('..\AProblemDefinition');
    
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
            Hmu = Probs.XP;    Hsd = Probs.sd;
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
        gx = Probs.Fung(Xpoint);
        Ifail = gx<0;
        Xfail = Xpoint(Ifail,:);
        fxhx_no = 1;
        Eu = 1; Betaint = 1 ;
    elseif strcmp(Method,'WIS') == 1 || strcmp(Method,'WIS_BOC_Ave') == 1 || strcmp(Method,'WIS_BOC_Var') == 1 || strcmp(Method,'WIS_BOC_Cov') == 1  % IS
        gx = Probs.Fung(Xpoint);
        Ifail = gx<0;
        Xfail = Xpoint(Ifail,:);
       [Nfail, ~]= size(Xfail);
        fxhxi = ones(Nfail,Probs.Nx);
        for i = Probs.xuLocal   % Probs.Nd+1 : Probs.Nx  updated: 2021-05-30
            fxi = pdf(Probs.Dist{i},Xfail(:,i), Probs.Para{i}(1),Probs.Para{i}(2));
%             hxi = normpdf(Xfail(:,i), Probs.XP(i),Probs.Hsd(i));
            hxi = normpdf(Xfail(:,i), Probs.XP(i),Hsd(i));
            fxhxi(:,i) = fxi./hxi;
        end
        fxhx_no = prod(fxhxi,2);
        Eu = 1; Betaint = 1 ;
        
    elseif strcmp(Method,'WLS') == 1 || strcmp(Method,'WLS_BOC_Ave') == 1 || strcmp(Method,'WLS_BOC_Var') == 1 || strcmp(Method,'WLS_BOC_Cov') == 1  % LS
        uP = (Probs.XP-Hmu)./Hsd;    % 在h(x)空间中求解重要方向    
        Beta = sqrt(sum(uP.^2));    % 可靠度指标
        Eu = uP/Beta ;              % 单位重要方向

        c1 = Beta-0.5;    c2 = Beta;  c3 = Beta+1.0;
        for j=1:N  
            ui = (Xpoint(j,:)- Hmu)./Hsd;
            Vuj = ui-(ui*Eu')*Eu ;       % 求垂直于单位重要方向的向量

            DOT1(j,:)=(c1*Eu+Vuj).*Hsd + Hmu;
            DOT2(j,:)=(c2*Eu+Vuj).*Hsd + Hmu;
            DOT3(j,:)=(c3*Eu+Vuj).*Hsd + Hmu;
            g1 = Probs.Fung(DOT1(j,:));
            g2 = Probs.Fung(DOT2(j,:));
            g3 = Probs.Fung(DOT3(j,:));

            % 三点插值
            betaj(j,:)=c1.*g2.*g3./((g1-g3).*(g1-g2))+c2.*g1.*g3...
                ./((g2-g3).*(g2-g1))+c3.*g2.*g1./((g3-g1).*(g3-g2));

            % save
            Xfail(j,:) = Xpoint(j,:);   % 抽取的样本
            % Uint (j,:) = ui + (betaj(j,:)-(ui*Eu'))*Eu; % updated: 2021-04-26, see MSSP Yuan Eq.23
            % Xint (j,:) = Uint (j,:).*Hsd + Hmu; 
            
            fxhx_no = 1;
        end 
        Betaint  = betaj;
    end
            
    % return
    Xfail_dv = Xfail(:,Probs.dvLocal);
    Betaint ;
    fxhx_no;
        
        

    