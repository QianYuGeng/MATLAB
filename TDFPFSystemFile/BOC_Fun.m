% function  [PFs_BOC, CovPFs_BOC]= BOC_Fun(PFs,CovPFs)
% Investigating the performance of different BOC algorithms: BOC-AVE, BOC-Var, BOC-Cov
% 2021-08-16


        PFs = [3.5e-5,  1.9e-7];
        CovPFs = [0.16,  0.91];

%         ncomb = 1;
%         [~,nks] = size(PFs);   % PFs = [] nsita * nks
% 
%         wk =  CovPFs(:,ncomb:nks).^(-2)./sum(CovPFs(:,ncomb:nks).^(-2),2);
%         PFs_BOC = sum(wk.*PFs(:,ncomb:nks),2);
%         VarPFs_BOC = sum(wk.^2.*(CovPFs(:,ncomb:nks).*PFs(:,ncomb:nks)).^2,2);
%         CovPFs_BOC = sqrt(VarPFs_BOC)./PFs_BOC;
        
        legendtext = {'k=1','k=2','k=3','k=4','k=5','k=6','k=7','k=8','k=9','k=10','k=11','k=12','k=13','k=14','k=15','k=16','k=17','k=18','k=19','k=20','k=21'};
        figclo = {'k--+','k:v','k-.x','b--s','b:d','b-.^','r-->','r:<','r-.p','c--h','c:*','c-.','g--','g:','g-.','y--','y:','y-.','m--','m:','m-.',}; % mark: o + * . d x s d ^ v < > p h
        LineWid = 1;                FontSZ = 10;
        
        Cave = @(r,c1,c2) sqrt(c1^2*r.^2+c2^2)./(r+1);
        Cvar = @(r,c1,c2) c1.*c2.*sqrt(c1^2.*r.^2+c2^2)./(c1^2.*r+c2^2);
        Ccov = @(r,c1,c2) c1.*c2.*sqrt(c2^2.*r.^2+c1^2)./(c2^2.*r+c1^2);
        
        r = linspace(0,1,100);
        C1= [0.1, 0.1 0.3]; C2 =[0.3, 0.1  0.1];
        ncase = length(C1);
        figure
        for i = 1:ncase
            subplot(ncase,1,i)
            c1=C1(i); c2=C2(i);
            y1 = Cave(r,c1,c2);
            y2 = Cvar(r,c1,c2);
            y3 = Ccov(r,c1,c2);
            plot(r,y1,figclo{1},r,y2,figclo{2},r,y3,figclo{3},'LineWidth',LineWid,'MarkerIndices',1:10:100);
            xlabel('The ratio of failure probability r');    ylabel('CovP_F(\theta)');
            legend('BOC-Ave','BOC-Var','BOC-Cov','FontSize', FontSZ);
            axis([0  1 0.0 0.3]) 
            if i ==1;     title('c_1=0.1,c_2=0.3');
            elseif i ==2; title('c_1=0.1,c_2=0.1');
            else          title('c_1=0.3,c_2=0.1'); 
            end
        end
        
        

