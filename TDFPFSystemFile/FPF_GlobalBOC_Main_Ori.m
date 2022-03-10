% function [niter, NMPP] = FPF_GlobalBOC_Main(NE,Case,N,Method,Seed0)    
% Global FPF by BOC and adaptive algorithm 
% updeated : 2021-05-29
% updated: 2021-09-15 put InstrumentalPDF_IS LS in one file
 

     clear all;close all;
     NE = 111;  N = 100;        % Ex. 1
    % NE = 1005; N = 100;         % Ex. 2
    % NE = 123; N = 200;         % Ex. 2
    % NE = 2020; N = 100;        % Ex. 3
    
    % NE = 111; N = 300;         % WLS, Ex. 1
    % NE = 124; N = 300;         % WLS, Ex. 1
    
    Method  = 'WIS_BOC_Cov';    % Methods: 'WIS_BOC_Var', 'WIS_BOC_Ave'
    Case    = 'CaseB';          % 'CaseA', 'CaseC'
    Seed0   = 65;               % Seed for random nuber generation

    PlotChoice = 1;     % 1: plot the FPF for each iteration, 0: not plot
    Hfactor    = 1;     % the enlarge factor of H(x) 
    Cov_tol    = 0.2;   % the stop creterior
    Kmax       = 10;     % the max number of iterations
    
    LineWid = 1; FontSZ     = 12;    % the size in figures
    
    addpath('..\AProblemDefinition');  % include the file OPTModel.m and MCS-Step.m   % addpath('./qcm')
    global Prob  
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
    
    rng(Seed0); 
     % generate the samples of sita, candidate points for minimizing the CovPF 
    if Prob.Nd == 1
        Ns = 200;
        Sseeds(:,1) = linspace(lb(:,1),ub(:,1),Ns);
        % Sseeds(:,1) = unifrnd(lb(:,1),ub(:,1),Ns,1);
    elseif Prob.Nd == 2
        Sseeds = []; Ns = 50;
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
    flag = 0;       Covmax(1) = 1; 
    while k<= Kmax  && flag~=1 
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
            [Probs.XP, PfAFORM,BetaAFORM,Ncall(k)]= PF_AFORM_Norm_Main(Probs); 
        end
        Probk{k} = Probs;     % save all setting into Probk
        
        % Sampling based on Probs
        [Xfail{k},fxhx_no{k},Eu{k},Betaint{k}] = InstrumentalPDF_WA(Probs,N,Method); % 2021-08-13
        fprintf('The number of failure samples Nf is : \t Nf=%d\n\n', size(Xfail{k},1));

        % find the next sampling center which minimizing the Cov
        [PFseeds,   CovPFseeds{k}]    = PFsita_WA_BOC(Sseeds,Xfail,fxhx_no,Eu,Betaint,Probk,N,Method,'Pf&Cov');      
        [Covmax(k+1), Idmax(k+1) ]    = max(CovPFseeds{k});
        Sopt(k+1,:) = Sseeds(Idmax(k+1),:);
        
        if Covmax(k) < Cov_tol % stop 
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
       
    %% plot optimal results
    if PlotChoice == 1
        % load('NE5N2000C3Kmax10Hfact10Rng36.mat');
        % load('NE2020N1000C2Kmax10Hfact10Rng38.mat');
        Plot_FPF_iterations(Xfail,fxhx_no,Eu,Betaint,Probk,N,Method,Sopt,Covmax);
    end    
    
    % plot the covPF Vs iteration 
    if Prob.Nd == 20
        figure
        if niter > 7,      ks_plot = [1,4,7,10];% 1:3:nks; 
        else,              ks_plot = 1:1:nks; 
        end
        n_plot = length(ks_plot);
        for j = 1:n_plot % 1: niter
            i = ks_plot(j); 
            subplot(ceil(n_plot/2),2,j);
            ZCovPFseeds = reshape(CovPFseeds{i},[Ns,Ns]);
            % FigCov = pcolor(StempX,StempY, ZCovPFseeds);
            % FigCov.FaceColor = 'interp';
            % FigCov.LineStyle = 'none';
            FigCov = imagesc(Stemp1,Stemp2, ZCovPFseeds);
            colorbar;  % ('location','SouthOutside')
            hold on
            % imagesc(Sopt(k-1,1),Sopt(k-1,2), Covmax(k-1));
            scatter(Sopt(i+1,1),Sopt(i+1,2), 58,'r','filled');  % 
            xlabel('\theta_1', 'FontSize', FontSZ);    ylabel('\theta_2','FontSize', FontSZ);
            title([ num2str(i) '-th iteration'],'FontSize', FontSZ);
        end
        saveas(gcf, ['NE' num2str(Probk{1}.NE) 'CovPF_Vs_niter.fig']);
    end
    
