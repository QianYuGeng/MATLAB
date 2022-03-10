% function BOC_Compare()
% verify the out-performance
% 2021-05-03

    clear
    addpath('..\AProblemDefinition');  % include the file OPTModel.m and MCS-Step.m

     NE = 4;    N = [20:20:100]; % wall
    % NE = 111;    N = 500:500:1500; ncomb = 3*ones(size(N));  % roof 
    % NE = 123;    N = 1000; ncomb = 2*ones(size(N));  % axle
    % NE = 10;    N = 400:300:1600; ncomb = 6*ones(size(N)); % tenbar
    % Case = 'CaseC'; 
    % ncomb = 1:5;  N = 100*ones(size(ncomb));
    % Seed0 = 56;
    Seed0 = 37;
    Method = 'WIS_BOC_Cov';
    LineMark = {'k:o','b-v','r-.d'};
    LineWid = 1; FontSZ     = 10;    % the size in figures
 
    for k = 1:length(N)
        Case = 'CaseA'; Seed0 = 26;
        [niter1(k),Ncalls1(k)]= FPF_GlobalBOC_Main(NE,Case,N(k),Method,Seed0);
        Case = 'CaseB'; Seed0 = 26;
        [niter2(k),Ncalls2(k)]= FPF_GlobalBOC_Main(NE,Case,N(k),Method,Seed0);
        Case = 'CaseC'; Seed0 = 26;
        [niter3(k),Ncalls3(k)]= FPF_GlobalBOC_Main(NE,Case,N(k),Method,Seed0);
    end
    
    figure
    subplot(3,1,1);
    plot(N,niter1,LineMark{1},'LineWidth',LineWid);
    hold on
    plot(N,niter2,LineMark{2},'LineWidth',LineWid);
    plot(N,niter3,LineMark{3},'LineWidth',LineWid);
    xlabel('Number of samples in each iteration N' ,'FontSize', FontSZ);   
    ylabel('Number of iterations','FontSize', FontSZ);
    legend('Initial design 1','Initial design 2','Initial design 3','FontSize', FontSZ);
    ylim([1, 13]);
    
    subplot(3,1,2);
    plot(N,N.*niter1,LineMark{1},'LineWidth',LineWid);
    hold on
    plot(N,N.*niter2,LineMark{2},'LineWidth',LineWid);
    plot(N,N.*niter3,LineMark{3},'LineWidth',LineWid);
    xlabel('Number of samples in each iteration N' ,'FontSize', FontSZ);   
    ylabel('IS samples N_T','FontSize', FontSZ);
    
    subplot(3,1,3);
    plot(N,N.*niter1+Ncalls1,LineMark{1},'LineWidth',LineWid);
    hold on
    plot(N,N.*niter2+Ncalls2,LineMark{2},'LineWidth',LineWid);
    plot(N,N.*niter3+Ncalls3,LineMark{3},'LineWidth',LineWid);
    xlabel('Number of samples in each iteration N' ,'FontSize', FontSZ);   
    ylabel('Total number of calls N_{all}','FontSize', FontSZ);

    save('NE4_niter_Vs_N20-100Rng26.mat')

