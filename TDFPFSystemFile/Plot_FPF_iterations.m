function [PfsitaRes,PFs_BOC,CovPFs_BOC]=Plot_FPF_iterations(Xfail,fxhx_no,Eu,Betaint,Probk,N,Method,Sopt,Covmax) 
% Plot the FPF in each iteration
% update : 2021-04-25
% updated: 2021-05-29

    PfsitaRes=[];

    addpath('..\AProblemDefinition');  % include the MCS/IS/LS methods for comparesion
   global Prob nowtime Pfsita_all Nmcs flagPlot% PFs_BOC CovPFs_BOC
%     flagPlot = 3; % 1: plot single PF(s), 2: plot combined PFC(s); 3: plot methods comparing WIS or ASI_IS; 
    flagMCS  = 1; % 0: plot without MCS exact point-wise values;  1: plot with exact values
    
    flagMat  = 0; % 0: no Exact result beforehand, need MCS; 1: yes, we have, no MCS needed. 
    flagMatIS= 0; % 0: no IS and ASI result beforehand ; 1: yes, we have, no IS needed. 
    
    WISflag=1;
    
    if Probk{1}.NE == 2020 || Probk{1}.NE == 1005
        MethodsExact = 'IS';
    else
        MethodsExact = 'MCS';
    end
  %%  
    if flagPlot  == 3  % ASI, WIS
        N_WA = N*length(Probk)/2;    %  using the same total number of samples to compare
        if flagMatIS == 0
            % WIS
            % [Xfail_WIS,   fzhz_WIS] = Augment_Sampling_IS(Probk{1},N_IS,'WIS');
%             [Xfail_WIS, fzhz_WIS,~,~] = InstrumentalPDF_IS(Probk{1},N_WA,'WIS');
%             [Xfail_WIS, fzhz_WIS,  Eu_WLS, Uint_WLS] = InstrumentalPDF_WA(Probk{1},N_WA,'WIS');
        if WISflag==1
                [Xfail_WIS, fzhz_WIS,  Eu_WLS, Uint_WLS] = InstrumentalPDF_WA_t(Probk{1},N_WA,'WIS',ceil(Prob.Nt/2),false(N_WA,1));
        end
                % ASI
% %             [Xfail_ASI_IS,fzhz_ASI] = Augment_Sampling_IS(Probk{1},N_WA,'ASI');
            % save(['WIS_ASI_NE' num2str(Probk{1}.NE) '.mat' ],'Xfail_WIS','Xfail_ASI_IS','fzhz_WIS','fzhz_ASI');
        else
            % load(['WIS_ASI_NE' num2str(Probk{1}.NE) '.mat' ],'Xfail_WIS','Xfail_ASI_IS','fzhz_WIS','fzhz_ASI');
        end
    end
    
    %% plot FPF one-dimension plot
    legendtext = {'k=1','k=2','k=3','k=4','k=5','k=6','k=7','k=8','k=9','k=10','k=11','k=12','k=13','k=14','k=15','k=16','k=17','k=18','k=19','k=20','k=21','k=22','k=23','k=24','k=25','k=26','k=27''k=28''k=29''k=30'};
    figclo = {'k--+','k:v','k-.x','b--s','b:d','b-.^','r-->','r:<','r-.p','c--h','c:*','c-.','g--','g:','g-.','y--','y:','y-.','m--','m:','m-.',}; % mark: o + * . d x s d ^ v < > p h
    LineWid = 1;                FontSZ = 10;
    Nkey = 22;                % for exact results 
    MarkerIndex = 1:2:Nkey ;
    [~,nks] = size(Xfail);   % nks = number of iterations 
    for PlotID = 1:Prob.Nd
        dv = linspace(Prob.dvDomain{PlotID}(1),Prob.dvDomain{PlotID}(2),Nkey);
        Sita = repmat(Prob.S0,size(dv,2),1);   % update: 2021-05-29
        Sita(:,PlotID) = dv';
        
        if flagMCS == 1  % exact values of FPF
            dvmcs = dv([1:3:Nkey]);
            Sitamcs =  Sita([1:3:Nkey],:);
            %Nmcs = 1e6;
            Nss = 1000; Nis = 1000; Nls = 1000;
            if flagMat  == 1  % load results
                load( ['Pf_ExactNE' num2str(Probk{1}.NE) 'PlotID' num2str(PlotID) '.mat'],'Pfsita','CovPfsita');
            elseif ~isempty(Pfsita_all)
                Pfsita=Pfsita_all((PlotID-1)*length(Pfsita_all)/Prob.Nd+1:PlotID*length(Pfsita_all)/Prob.Nd,:);
            else % no exact result file
                for k = 1:size(Sitamcs,1)
                    % Probs = rmfield(Prob,'mu');   Probs.mu = Smu(k,:); Probs = ParaStat(Probs,'StoP');
                    Probmcs = rmfield(Prob,'mu');
                    for i = 1:Prob.nxr 
                        si_Para = Prob.si_Para{i};        % index of s_i in Para, [1, 2]
                        si_Sita = Prob.si_Sita{i};        % index of s_i in Sita, [3]
                        ixr = Prob.dvLocal(i);              % index of x_r in x
                        Probmcs.Para{ixr}(:,si_Para) = Sitamcs(k,si_Sita);  % put s_i into ParaN = {[N,nx], [N,nx]}       
                    end
                    Probmcs = PrepareProb(Probmcs);

                    if    strcmp(MethodsExact,'MCS') == 1      
                        [ Pfsita(k,:), CovPfsita(k,:)] = TDPF_MCS_Main(Probmcs,Nmcs,[0,Probmcs.dt*floor(Probmcs.Nt/2),Probmcs.tmax]);
                        %换成TDPF
                    elseif strcmp(MethodsExact,'SS') == 1    
                        [ Pfsita(k), CovPfsita(k)] = PF_SS_Main(Probmcs,Nss);
                    elseif strcmp(MethodsExact,'IS') == 1    
                        if Probmcs.NE == 2020
                           [ ~, ~, Pfsita(k), CovPfsita(k)] = BenchProb2_AuIS(Probmcs,Nis);
                        else
                           [ Pfsita(k), CovPfsita(k)] = PF_IS_Main(Probmcs,Nis);
                        end
                    elseif strcmp(MethodsExact,'LS') == 1    
                        [ Pfsita(k), CovPfsita(k)] = PF_LS_Main(Probmcs,Nls);
                    else
                        error('Method exact selection is wrong! ');
                    end
                end
                Pfsita(Pfsita==0)=1;
                
                % output to file
                fid = fopen([ 'ExactPointwise_resultsNE' num2str(Probmcs.NE) '.txt'],'at+');
                fprintf(fid,'\nMethodsExact = %s\nPlotID = %d \t sita_i =  %12.4e \n',MethodsExact,PlotID,Sitamcs(k,si_Sita));
                fprintf(fid,'Nmcs = %12d\tNss = %12d\tNis = %12d\tNls = %12d\t\n',Nmcs, Nss, Nis,Nls);
                fprintf(fid,'Pf = %12.4e(CovPf = %12.4f)\n',Pfsita(k), CovPfsita(k));
                fclose(fid);
%                 save( [nowtime 'Pf_ExactNE' num2str(Probmcs.NE) 'PlotID' num2str(PlotID) '.mat'],'Pfsita','CovPfsita');
            end
        end

        %% plot the figures
        figure('Name',nowtime)
        if nks > 7,      ks_plot =1:5:nks;% [1,3,6,8,10];%  
        else,            ks_plot = 1:1:nks; 
        end
        if flagPlot == 1  % 1: single iteration PF(s) not BOC
            subplot(2,1,1);
            % FPF plot 
            for i = 1:length(ks_plot)
                ks = ks_plot(i);
                % [PFs(:,ks), CovPFs(:,ks)] = PFsita_WIS(Sita,Xfail{ks},fxhx_no{ks},Probk{ks},N,'Pf&Cov');  % single iteration PF(s)
                [PFs(:,ks), CovPFs(:,ks)] = PFsita_WA(Sita,Xfail{ks},fxhx_no{ks},Eu{ks},Betaint{ks},Probk{ks},N,Method,'Pf&Cov');  % single iteration PF(s) % 2021-09-17
                semilogy(dv,PFs(:,ks),figclo{i},'LineWidth',LineWid,'MarkerIndices',MarkerIndex); hold on;
            end
            [PFs_BOC,   CovPFs_BOC]     = PFsita_WA_BOC(Sita,Xfail,fxhx_no,Eu,Betaint,Probk,N,Method,'Pf&Cov'); 
            semilogy(dv,PFs_BOC,figclo{length(ks_plot)+1},'LineWidth',LineWid,'MarkerIndices',MarkerIndex);
            if flagMCS == 1  
                semilogy(dvmcs,Pfsita,'o','MarkerFaceColor','r'); % 
                legend(legendtext{ks_plot},'BOC',MethodsExact ,'FontSize', FontSZ);
            else
                legend(legendtext{ks_plot},'BOC','FontSize', FontSZ);
            end
            xlabel(['\theta_' int2str(PlotID)],'FontSize', FontSZ);    ylabel(['P_F(\theta_' int2str(PlotID) ')'],'FontSize', FontSZ);
            
            % covPF plot
            subplot(2,1,2);
            for i = 1:length(ks_plot)
                ks = ks_plot(i);
                plot(dv,CovPFs(:,ks),figclo{i},'LineWidth',LineWid,'MarkerIndices',MarkerIndex); hold on;
                hold on
            end
            plot(dv,CovPFs_BOC,figclo{length(ks_plot)+1},'LineWidth',LineWid,'MarkerIndices',MarkerIndex);
            xlabel(['\theta_' int2str(PlotID)],'FontSize', FontSZ);    ylabel(['Cov[P_F(\theta_' int2str(PlotID) ')]'],'FontSize', FontSZ);
            legend(legendtext{ks_plot},'BOC','FontSize', FontSZ);
        
        elseif  flagPlot == 2  % 2: plot BOC results Vs niter (not single WIS)
            subplot(2,1,1);
            % FPF plot % BOC PFC(s)
            for i = 1:length(ks_plot)
                ks = ks_plot(i);
                [PFs_BOC(:,ks),   CovPFs_BOC(:,ks)]  = PFsita_WA_BOC(Sita,Xfail(1,1:ks),fxhx_no(1,1:ks),Eu(1:ks),Betaint(1:ks),Probk(1:ks),N,Method,'Pf&Cov'); 
               %  [PFs_BOC(:,ks),   CovPFs_BOC(:,ks)]  = TDPFsita_WA_BOC(Sita,Xfail(1:ks),fxhx_no(1:ks),Eu(1:ks),Betaint(1:ks),Probk(1:ks),N,Method,'Pf&Cov'); 
               semilogy(dv,PFs_BOC(:,ks),figclo{i},'LineWidth',LineWid,'MarkerIndices',MarkerIndex);
                hold on
            end
            if flagMCS == 1  
%                 semilogy(dvmcs,Pfsita,'o','MarkerFaceColor','r'); % ,'MarkerFaceColor','b'
                scatter(dvmcs,Pfsita(:,1),'o','MarkerFaceColor','r','DisplayName',MethodsExact);
%                 scatter(dvmcs,Pfsita(:,ceil(Prob.Nt/2)),'o','MarkerFaceColor','r','handlevisibility','off');
%                 scatter(dvmcs,Pfsita(:,end),'o','MarkerFaceColor','r','handlevisibility','off');
                legend(legendtext{ks_plot},MethodsExact,'FontSize', FontSZ);
            else
                legend(legendtext{ks_plot},'FontSize', FontSZ);
            end
            xlabel(['\theta_' int2str(PlotID)],'FontSize', FontSZ);    ylabel(['P_F(\theta_' int2str(PlotID) ')'],'FontSize', FontSZ);
            % cov plot
            subplot(2,1,2);
            for i = 1:length(ks_plot)
                ks = ks_plot(i);
                plot(dv,CovPFs_BOC(:,ks),figclo{i},'LineWidth',LineWid,'MarkerIndices',MarkerIndex)
                hold on
            end
            if Probk{1}.Nd == 1
                scatter(Sopt(ks_plot+1,PlotID),Covmax(ks_plot+1)','bd','filled','SizeData',55);    legend(legendtext{ks_plot},'\theta_H^{(k+1)}','FontSize', FontSZ);
            else
                legend(legendtext{ks_plot},'FontSize', FontSZ);
            end
            xlabel(['\theta_' int2str(PlotID)],'FontSize', FontSZ);    ylabel(['Cov[P_F(\theta_' int2str(PlotID) ')]'],'FontSize', FontSZ);
            saveas(gcf, [nowtime 'NE' num2str(Probk{1}.NE) 'iter.fig']);
            saveas(gcf, [nowtime 'NE' num2str(Probk{1}.NE) 'iter.png']);
        
        else % 3: Plot Comparsion of methods
            % FPF plot
            if WISflag==1
            subplot(2,1,1);hold on;
              [Pfsita2, CovPfsita2] = PFsita_WA(Sita,Xfail_WIS,fzhz_WIS,Eu_WLS,Uint_WLS,Probk{1},N_WA,Method,'Pf&Cov');
            semilogy(dv,Pfsita2,figclo{1},'LineWidth',LineWid,'MarkerIndices',MarkerIndex,'DisplayName','WIS'); 
              subplot(2,1,2);hold on;
            plot(dv,CovPfsita2,figclo{1},'LineWidth',LineWid,'MarkerIndices',MarkerIndex); % ,'MarkerIndices',1:10:ns
            end
            subplot(2,1,1);hold on;title([nowtime 'NE=' num2str(Probk{1}.NE)]);
            set(gca,'yscale','log');
            % WIS and ASI
%             [Pfsita2, CovPfsita2 ] = FPF_ASI_IS (Sita,N,Xfail_WIS,fzhz_WIS,Probk{1},'WIS');
%             [Pfsita2, CovPfsita2] = PFsita_WIS(Sita,Xfail_WIS,fzhz_WIS,Probk{1},N_IS,'Pf&Cov');
%         [Pfsita2, CovPfsita2] = PFsita_WA(Sita,Xfail_WIS,fzhz_WIS,Eu_WLS,Uint_WLS,Probk{1},N_WA,Method,'Pf&Cov');  % single iteration PF(s) % 2021-09-17
        
% %         [Pfsita5, CovPfsita5] = FPF_ASI_IS (Sita,N_WA,Xfail_ASI_IS,fzhz_ASI,Probk{1},'ASI');
%         semilogy(dv,Pfsita2,figclo{1},'LineWidth',LineWid,'MarkerIndices',MarkerIndex,'DisplayName','WIS'); % ,'MarkerIndices',1:10:ns
            hold on
% %         semilogy(dv,Pfsita5 ,figclo{2},'LineWidth',LineWid,'MarkerIndices',MarkerIndex); 
            xlabel(['\theta_' num2str(PlotID)]);    ylabel(['P_F(\theta_',num2str(PlotID), ')']);
%             ylim([1e-6 1e4])
           
            % BOC PFC(s)
%             [PFs_BOC,   CovPFs_BOC]  = PFsita_WA_BOC(Sita,Xfail(1:end),fxhx_no(1:end),Eu,Betaint,Probk,N,Method,'Pf&Cov'); 
            [PFs_BOC,   CovPFs_BOC]  = TDPFsita_WA_BOC(Sita,Xfail,fxhx_no,Eu,Betaint,Probk,N,Method,'Pf&Cov');
            
            semilogy(dv,PFs_BOC(:,1),figclo{3},'LineWidth',LineWid,'MarkerIndices',MarkerIndex,'DisplayName','BOC-AWIS t=t_0');
            semilogy(dv,PFs_BOC(:,ceil(Prob.Nt/2)),figclo{3},'LineWidth',LineWid,'MarkerIndices',MarkerIndex,'DisplayName','BOC-AWIS t=t_{mid}');%'handlevisibility','off');
            semilogy(dv,PFs_BOC(:,end),figclo{3},'LineWidth',LineWid,'MarkerIndices',MarkerIndex,'DisplayName','BOC-AWIS t=t=t_{max}');%'handlevisibility','off');
            if flagMCS == 1  
%                 semilogy(dvmcs,Pfsita,'o','MarkerFaceColor','r','DisplayName',MethodsExact); % ,'MarkerFaceColor','b'
%                scatter(dvmcs,Pfsita,'o','MarkerFaceColor','r','DisplayName',MethodsExact);
                scatter(dvmcs,Pfsita(:,1),'o','MarkerFaceColor','r','DisplayName',MethodsExact);
                scatter(dvmcs,Pfsita(:,ceil(size(Pfsita,2)/2)),'o','MarkerFaceColor','r','handlevisibility','off');
                scatter(dvmcs,Pfsita(:,end),'o','MarkerFaceColor','r','handlevisibility','off');
%                 legend('WIS','BOC-AWIS','FontSize', FontSZ);
                legend('FontSize', FontSZ);
            else
                legend('WIS','ASI-IS','AWIS','FontSize', FontSZ);
            end
%             PfsitaRes=[PfsitaRes;Pfsita];
            
            % cov plot
            subplot(2,1,2);hold on;
%             plot(dv,CovPfsita5,figclo{2},'LineWidth',LineWid,'MarkerIndices',MarkerIndex); 
%             plot(dv,CovPFs_BOC,figclo{3},'LineWidth',LineWid,'MarkerIndices',MarkerIndex);
            plot(dv,CovPFs_BOC(:,1),figclo{3},'LineWidth',LineWid,'MarkerIndices',MarkerIndex);
            plot(dv,CovPFs_BOC(:,ceil(Prob.Nt/2)),figclo{3},'LineWidth',LineWid,'MarkerIndices',MarkerIndex);%'handlevisibility','off');
            plot(dv,CovPFs_BOC(:,end),figclo{3},'LineWidth',LineWid,'MarkerIndices',MarkerIndex);%'handlevisibility','off');
            
            xlabel(['\theta_' num2str(PlotID)]);    ylabel(['Cov[P_F(\theta_',num2str(PlotID), ')]']);
        end
        PfsitaRes=[PfsitaRes;Pfsita];
        
        saveas(gcf, [nowtime 'NE' num2str(Probk{1}.NE) 'PFtheta' num2str(PlotID) '.fig']);
        saveas(gcf, [nowtime 'NE' num2str(Probk{1}.NE) 'PFtheta' num2str(PlotID) '.png']);
    end
    
    
    %% two dimensional plot
    if Prob.Nd == 2
        % two dimension plot PF(s1,s2) using surf 
        [uX,uY] = meshgrid(linspace(Prob.lb(1),Prob.ub(1),Nkey),linspace(Prob.lb(2),Prob.ub(2),Nkey)); % Draw a two-dimensional picture
        Pfus = [];
        for i=1:size(uX,1)
            for j= 1:size(uY,2)
                Sita_2d = [uX(i,i) uY(j,j)];
                [Pfsita_2d(i,j), CovPfsita_2d(i,j) ] = PFsita_WA_BOC(Sita_2d,Xfail(1,:),fxhx_no(1,:),Eu,Betaint,Probk,N,Method,'Pf&Cov'); 
            end
        end
        % [Pfsitamcs_2d(i,j), CovPfsitamcs_2d(i,j)] = MCS_PointFPF_Sita(Sita_2d,Nmcs); 
                
        figure
        subplot(1,2,1);
        surf(uX,uY,Pfsita_2d)
        xlabel('\theta_1','FontSize', FontSZ); ylabel('\theta_2','FontSize', FontSZ);   zlabel('P_F(\theta_1,\theta_2)','FontSize', FontSZ);
        set(gca,'ZScale','log')
        subplot(1,2,2);
        surf(uX,uY,CovPfsita_2d)
        xlabel('\theta_1','FontSize', FontSZ); ylabel('\theta_2','FontSize', FontSZ);   zlabel('Cov[P_F(\theta_1,\theta_2)]','FontSize', FontSZ);
        
        saveas(gcf, [nowtime 'NE' num2str(Probk{1}.NE) 'PFtheta12.fig']);
        % saveas(gcf, ['NE' num2str(Probk{1}.NE) 'PFtheta12.pdf']);
    end
%%
    
    if 0
     legend('k=1 AWIS','k=3 AWIS','k=6 AWIS','k=8 AWIS','Exact','FontSize', 10);

        
        legend('k=1 AWIS','k=3 AWIS','k=6 AWIS','k=8 AWIS','Exact','FontSize', 10);
     legend('k=1 AWIS','k=3 AWIS','k=6 AWIS','k=8 AWIS','\theta_H^{(k+1)}','FontSize', 10);

     legend('k=1 WIS','k=3 WIS','k=6 WIS','k=8 WIS','k=8 AWIS','Exact','FontSize', 10);
     legend('k=1 WIS','k=3 WIS','k=6 WIS','k=8 WIS','k=8 AWIS','Exact','FontSize', 10);

     legend('WIS','ASI-IS','AWIS','Exact','FontSize', 10);
    end
    
    