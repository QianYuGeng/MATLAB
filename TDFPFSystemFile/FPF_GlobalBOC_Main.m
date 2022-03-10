% function [niter, NMPP] = FPF_GlobalBOC_Main(NE,Case,N,Method,Seed0)    
% Global FPF by BOC and adaptive algorithm 
% updeated : 2021-05-29
% updated: 2021-09-15 put InstrumentalPDF_IS LS in one file
 

     clear all;
%      load('03081629NE41101.mat','Pfsita');
load('03092045NE60301.mat','Pfsita');
     %figure(19);subplot(2,1,1);delete(get(gca,'title'));title(' ');legend('BOC-AWIS');
     %close all;
%      NE=5501;N=500;
%      NE=41101;N=500;
%      NE=5201;N=100;
     NE=060301;N=500;
%      NE=5501;N=50;
%      NE=5201;N=300;
%      NE=5201;N=300;
%      NE=41101;N=100;
%      NE=101;N=100;
%      NE=40101;N=100;
%         NE=401;N=100; 
%      NE=1401;N=300; 
% NE=11101;N=300; 
%      NE = 111;  N = 100;        % Ex. 1
    % NE = 1005; N = 100;         % Ex. 2
    % NE = 123; N = 200;         % Ex. 2
    % NE = 2020; N = 100;        % Ex. 3
    
    % NE = 111; N = 300;         % WLS, Ex. 1
    % NE = 124; N = 300;         % WLS, Ex. 1
    
%     figure(20);subplot(2,1,1);ylim([1e-3,1e4]);
%figure(15);view(90,0);
    
    Method  = 'WIS_BOC_Cov';    % Methods: 'WIS_BOC_Var', 'WIS_BOC_Ave'
    Case    = 'CaseB';          % 'CaseA', 'CaseC'
    Seed0=unidrnd(100);
%     Seed0   =55;               % Seed for random nuber generation
    XPflag=2;
% for XPflag=1:2
%     XPflag=3;
    PlotChoice = 1;     % 1: plot the FPF for each iteration, 0: not plot
    Hfactor    = 1;     % the enlarge factor of H(x) 
    Cov_tol    = 0.2;   % the stop creterior
    Kmax       = 40;     % the max number of iterations
    
    LineWid = 1; FontSZ     = 12;    % the size in figures
    
    addpath('..\AProblemDefinition');  % include the file OPTModel.m and MCS-Step.m   % addpath('./qcm')
    global Prob Nmcs
    Nmcs=1e7;
    Prob = Problem_FPF(NE); 
    Idv = Prob.dvLocal ; 
    lb = Prob.lb;       ub = Prob.ub;       
    
    if     strcmp(Case,'CaseA') == 1 
        Sopt(1,:) = Prob.lb;    % \theta^(0), initial design
    elseif strcmp(Case,'CaseB') == 1 
        Sopt(1,:) = (Prob.lb+Prob.ub)/2 ; 
    elseif strcmp(Case,'CaseC') == 1 
        Sopt(1,:) = Prob.ub; 
    end
   %% 
    rng(Seed0); 
     % generate the samples of sita, candidate points for minimizing the CovPF 
    if Prob.Nd == 1
        Ns = 200;
        Sseeds(:,1) = linspace(lb(:,1),ub(:,1),Ns);
        % Sseeds(:,1) = unifrnd(lb(:,1),ub(:,1),Ns,1);
    elseif Prob.Nd == 2
        Sseeds = []; Ns = 50;%Ns = 50;
        Stemp1 = linspace(lb(:,1),ub(:,1),Ns); Stemp2 = linspace(lb(:,2),ub(:,2),Ns);
        [StempX,StempY] = meshgrid(Stemp1,Stemp2);
        for i = 1:Ns
            for j = 1:Ns
                % Sseeds = [Sseeds; StempX(i,j),StempY(i,j)];  !!!
                Sseeds = [Sseeds; StempX(j,i),StempY(j,i)];
            end
        end
    else  % Nd>3, using the unif distribution
        Ns = 1000;
        for i=1: Prob.Nd
            Sseeds(:,i) = unifrnd(lb(:,i),ub(:,i),Ns,1);
        end
    end
    
    k = 1; 
    flag = 0;       Covmax(1) = 1;  tidmax(1)=1; t_all=0:Prob.dt:Prob.tmax;
    while k<= Kmax  && (flag~=1|| k<=2/3*Kmax ) 
        Sita = Sopt(k,:);
        fprintf('\nThe %d-th sampling center : ',k);fprintf('%6.4f\t',Sopt(k,:)); fprintf('\n');
        
        % instrumental PDF setting, updated into Probs 
        Probs = rmfield(Prob,{'mu','sd'});
        for i = 1:Prob.nxr 
            si_Para = Prob.si_Para{i};        % index of s_i in Para, [1, 2]
            si_Sita = Prob.si_Sita{i};        % index of s_i in Sita, [3]
            ixr = Idv(i);            % index of x_r in x
            Probs.Para{ixr}(:,si_Para) = Sita(:,si_Sita);  % put s_i into ParaN = {[N,nx], [N,nx]}       
        end
        Probs = PrepareProb(Probs);
        if Hfactor ~= 1 % 2021-09-21 Hfactor effect into Probs
            Probs.sd(Idv) = Hfactor*Probs.sd(Idv);  
            Probs = rmfield(Probs,'Para');
            Probs = PrepareProb(Probs);
        end
        
        % updated the XP at each iteration
        if NE == 123
            [Probs.XP, PfAFORM,BetaAFORM]= PF_RF_Main(Probs); 
        elseif NE == 2020
            Probs.XP = Probs.mu;
        else
            
%             [Probs.XP,~,~,Ncall(k)]= PF_AFORM_Norm_Main(Probs); 
             if  XPflag==1
%                 [Probs.XP, PfAFORM(k),BetaAFORM,Ncall(k)]= PF_AFORM_Norm_Main(Probs);
%                 XP= TPF_AFORM_Norm_Main(Probs,0);
 [Probs.XP, PfAFORM(k),BetaAFORM,Ncall(k)]= TPF_AFORM_Norm_Main(Probs,0);
%                 [Probs.XP, PfAFORM(k),BetaAFORM,Ncall(k)]= TPF_AFORM_Norm_Main(Probs,t_all(ceil(Probs.Nt/2)));
%                 if XP~=Probs.XP
%                     pause;
%                 end
             elseif  XPflag==2
                [Probs.XP, PfAFORM,BetaAFORM,Ncall(k)]= TPF_AFORM_Norm_Main(Probs,t_all(tidmax(k))); 
%                 [Probs.XP, PfAFORM,BetaAFORM,Ncall(k)]= TPF_AFORM_Norm_Main(Probs,0);
             else%XPflag==3
                 [Probs.XP,~,~,Ncall(k)]= PF_AFORM_Norm_Main(Probs); 
                 for  iii=1:Prob.Nt
                    [Probs.XPt(iii,:), PfAFORM(iii),BetaAFORM(iii),Ncall(k)]= TPF_AFORM_Norm_Main(Probs,t_all(iii));
                 end
             end
        end
        Probk{k} = Probs;     % save all setting into Probk
        
        % Sampling based on Probs
        [Xfail(:,k),fxhx_no(:,k),Eu{k},Betaint{k},Ifail{k}] = TDInstPDF_WA(Probs,N,Method); % 2021-08-13
%         [Xfail{k},fxhx_no{k},Eu{k},Betaint{k}] = InstrumentalPDF_WA(Probs,N,Method); % 2021-08-13
        fprintf('The number of failure samples Nf is : \t Nf=%d\n\n', size(Xfail{k},1));

        % find the next sampling center which minimizing the Cov
        [PFseeds,   CovPFseeds{k}]    = TDPFsita_WA_BOC(Sseeds,Xfail,fxhx_no,Eu,Betaint,Probk,N,Method,'Pf&Cov');  
%         [PFseeds,   CovPFseeds{k}]    = PFsita_WA_BOC(Sseeds,Xfail,fxhx_no,Eu,Betaint,Probk,N,Method,'Pf&Cov');
%         [Covmax(k+1), Idmax(k+1) ]    = max(CovPFseeds{k}(:,1));
%         [Covmax(k+1), Idmax(k+1) ]    = max(CovPFseeds{k});
            [maxres,idres ]    = max(CovPFseeds{k});
         [maxmaxres,idmaxres]=max(maxres);
         Covmax(k+1)=maxres(idmaxres); Idmax(k+1)=idres(idmaxres);
         tidmax(k+1)=idmaxres;
            Sopt(k+1,:) = Sseeds(Idmax(k+1),:);
        
%         if Covmax(k) < Cov_tol % stop 
        if Covmax(k+1) < Cov_tol % stop 
                flag = 1;
        end
        k = k+1;
    end
    
    % return
    niter = k-1;
    if NE == 123 ||  NE == 2020,    NMPP  = 100;    
    else,                           NMPP  = sum(Ncall,'all'); % number of calls of LSF in XP solving
    end
    % Sopt = Sopt(1:niter,:);

    % save( ['NE' num2str(NE) 'N' num2str(N) 'C' num2str(Cov_tol*10) 'Kmax' num2str(Kmax)  'Hfact' num2str(Hfactor*10) 'Rng' num2str(Seed0) '.mat']) 
        %%
        global Pfsita_all nowtime
%        Pfsita_all=Pfsita;
        %%
        nowtime=datestr(now,'mmddHHMM'); 
    %plot optimal results
   
    if PlotChoice == 1
        % load('NE5N2000C3Kmax10Hfact10Rng36.mat');
        % load('NE2020N1000C2Kmax10Hfact10Rng38.mat');
        global flagPlot
        flagPlot=3;
        [Pfsita,PFs_BOC,CovPFs_BOC]=Plot_FPF_iterations(Xfail,fxhx_no,Eu,Betaint,Probk,N,Method,Sopt,Covmax);
        
    end   
    %%
    save([nowtime 'NE' num2str(Probk{1}.NE) '.mat']);
    %% plot the covPF Vs iteration   //Sita-Sita
    if Prob.Nd == 2
        figure('Name',nowtime);
        if niter > 7,      ks_plot =  1:5:niter;%[1,4,7,10] 
        else,              ks_plot = 1:1:nks; 
        end
        n_plot = length(ks_plot);    
        for j = 1:n_plot % 1: niter
            i = ks_plot(j); 
            %t loop
            for it=1:Prob.Nt
            subplot(n_plot,3,Prob.Nt*(j-1)+it); 
            ZCovPFseeds = reshape(CovPFseeds{i}(:,it),[Ns,Ns]);
            FigCov = imagesc(Stemp1,Stemp2, ZCovPFseeds);
            colorbar;  % ('location','SouthOutside')
            hold on
            % imagesc(Sopt(k-1,1),Sopt(k-1,2), Covmax(k-1));
            scatter(Sopt(i+1,1),Sopt(i+1,2), 58,'r','filled');  % 
            xlabel('\theta_1', 'FontSize', FontSZ);    ylabel('\theta_2','FontSize', FontSZ);
            title([ num2str(i) '-th iteration'],'FontSize', FontSZ);
            end
        end
        saveas(gcf, [datestr(now,'mmddHHMM') 'NE' num2str(Probk{1}.NE) 'CovPF_Vs_niter.fig']);
    end
    
     %% plot the covPF Vs iteration   //Sita-t
     if Prob.Nd == 1
        figure('Name',nowtime);
        if niter > 16  ,      ks_plot =  1:5:niter;%[1,4,7,10] 
        elseif niter > 7 && niter < 16  ,      ks_plot =  1:3:niter;%[1,4,7,10] 
        else,              ks_plot = 1:1:niter; 
        end
        n_plot = length(ks_plot);    
        for j = 1:n_plot % 1: niter
            i = ks_plot(j); 
            subplot(ceil(n_plot/2),2,j); 
            ZCovPFseeds = CovPFseeds{i};
            FigCov = imagesc(Sseeds,0:Prob.dt:Prob.tmax, ZCovPFseeds');
            colorbar;  % ('location','SouthOutside')
            hold on
            scatter(Sopt(i+1),t_all(tidmax(i+1)),58,'r','filled');  % 
            xlabel('\theta_1', 'FontSize', FontSZ);    ylabel('t','FontSize', FontSZ);
            title([ num2str(i) '-th iteration'],'FontSize', FontSZ);            
        end
        saveas(gcf, [nowtime 'NE' num2str(Probk{1}.NE) 'CovPF-t_Vs_niter.fig']);
     end
%%  
if Prob.Nd==1  
    [sfX,sfY]=meshgrid(0:Prob.dt:Prob.tmax,Sseeds);
         %PFs_BOC_all{1};%[22,2]
         
%          subplot(1,2,1);
%         surf(sfX,sfY,PFs_BOC_all{1})
%         xlabel('t','FontSize', FontSZ); ylabel(['\theta_i','i=',num2str(1)],'FontSize', FontSZ);   zlabel('P_F(\theta_i,t)','FontSize', FontSZ);
%         set(gca,'ZScale','log')
%         subplot(1,2,2);
           
%          for i=1:Prob.Nt
             figure('Name',nowtime);
                 surf(sfX,sfY,CovPFseeds{niter},'FaceAlpha',0.7,'EdgeColor','none');hold on;
            plot3(t_all(tidmax(end)),Sopt(end),Covmax(end),'r.','markersize',10);
            fprintf(' %12.4e %12.4e %12.4e \n  ',Sopt(end),t_all(tidmax(end)),Covmax(end));
%          end
        xlabel('t','FontSize', FontSZ); ylabel(['\theta_1' ],'FontSize', FontSZ);   zlabel('Cov[P_F(\theta_1,t)]','FontSize', FontSZ);
end
        saveas(gcf, [nowtime 'NE' num2str(Probk{1}.NE) 'CovPF-tSurf.fig']);
% end
