function Prob = Problem_FPF(NE)
% Problem definition
% update: 2020-08-10, change the definition of the RBDO problem, it can be any paramters % 1. Prob.dvLocal,Prob.xuLocal, Prob.si_Para, Prob.Domain_x_s are newly added.
% updated: 2021-05-29, Global FPF by BOC: chang nx to Nx, ns to Nd, dsdomain to dvDomain
    Prob.NE = NE; % NE 2021-08-10
    if NE == 5301       % # Example 1 % Two bar frame  
        Prob.Nr = 4;   % number of time independent random variables： D1  l1  l2
        Prob.NF = 1;   % Number of random loads 
        Prob.tmax = 10;           % Total years
        Prob.dt = 0.2;           % time interval
        nKL=fix(Prob.tmax/Prob.dt);
        Prob.nKL =nKL; 
        Prob.L = 1;               % L big，Cov big， Why？
        Prob.Nt = Prob.tmax/Prob.dt+1;   % total steps of time instance, YY

        Prob.Nx = Prob.Nr+Prob.nKL*Prob.NF;
        Prob.Dist = {'norm','norm','norm','logn'};
        for i = Prob.Nr+1:Prob.Nx
            Prob.Dist{i}= 'norm';
        end
        Prob.mu = [0.2      0.4     0.3   2.5e8  zeros(1,Prob.nKL*Prob.NF) ];  % mean : D1  l1  l2
        Prob.sd = [0.002, 0.004, 0.003    2.5e6  ones(1,Prob.nKL*Prob.NF) ];   % standard deviation 
        Prob.muF = [2.2e6];    % Mean value of dynamic load
        Prob.sdF = [2.2e5];

        Prob.Covfun = @(tao,L) {Prob.sdF.^2.*exp(-tao.^2./L^2)};  % 
        [Prob.B, Prob.SigmaB] = KLexpansion(Prob);
        Prob.ZpToFt = @(Zp,Prob) ZpToFt(Zp,Prob);
        Prob.YFunt = @(t) exp(-0.05.*t);% Degenerate formula of mean value of resistance
        
        Prob.Fungt  = @(x,t,F) pi.*x(:,3).*x(:,1).^2.*x(:,4).*Prob.YFunt(t)-4.*F(:,1).*sqrt(x(:,2).^2+x(:,3).^2);  % LSF, response function
       
         Prob.dvLocal =  [ 1 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1]};    
        Prob.si_Sita = {[1]};    
        Prob.Domain_x_s = {{[0.15 0.25]}};  
        
        Prob.PositiveCheck = 0 ;           % prior distribution of sita
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
    end
     
    if NE == 060301%非线性 
     %from zhang2021
        Prob.Nr  = 3;   % number of random variables锛? [h  w  E]
        Prob.NF  = 1;   % Number of random loads  
        
        Prob.tmax = 10;           % Total years
        Prob.dt = 1;           % time interval
        nKL=fix(Prob.tmax/Prob.dt);
        Prob.nKL =nKL; 
        Prob.L  = 1;
        Prob.Nt = Prob.tmax/Prob.dt+1;   % total steps of time instance, YY        
        Prob.Nx = Prob.Nr+Prob.nKL*Prob.NF;
        
        Prob.Dist = {'norm','norm','norm'};
        for i = 1:Prob.Nx     Prob.Dist{i}= 'norm';      end
        Prob.mu = [0.2 0.04 2.4*10^8 zeros(1,Prob.nKL*Prob.NF)];  % mean
        Prob.sd = [0.02 0.004 2.4*10^7 ones(1,Prob.nKL*Prob.NF)];   % standard deviation  
        Prob.muF = [3500  ];    % Mean value of dynamic load
        Prob.sdF = [700   ];
        
%           Prob.dvLocal =  [ 1 2];  % good
%         Prob.si_Para = {[1],[1]};    
%         Prob.si_Sita = {1,2};    
%         Prob.Domain_x_s = {{[0.15,0.25]},{[0.02, 0.04]}};  

        Prob.dvLocal =  [ 1  ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1]};    
        Prob.si_Sita = {[1] };    
        Prob.Domain_x_s = {{[0.1 0.3]}};  
        
        covfun1 = @(tao,L) Prob.sdF(1).^2.*exp(-tao.^2./L^2);
        Prob.Covfun = @(tao,L) {covfun1(tao,L)};%,covfun2(tao,L)};
        C0 = Prob.Covfun(0: Prob.dt : Prob.tmax,Prob.L);
        Prob.C=cell(1,Prob.NF);
        for i=1:Prob.NF
            Prob.C{i}=toeplitz(C0{i});
        end
        
        [Prob.B, Prob.SigmaB] = KLexpansion(Prob);
        Prob.ZpToFt = @(Zp,Prob) ZpToFt(Zp,Prob);
        Prob.Fungt = @(x,t,F) -(F(1,:)*5/4+7.85*10^4*x(:,1).*x(:,2).*25/8)+( ( x(:,1)-2*5*10^-4*t).*(x(:,2)-2*5*10^-4*t).^2.*x(:,3))/4;
        %Prob.Fu%对数正态h0 w0 sigy 即x(:,1) x(:,2) x(:,3);
        
%         Prob.dvLocal =  [ 1  ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1]};    
%         Prob.si_Sita = {[1] };    
%         Prob.Domain_x_s = {{[0.01 0.4]}};  
      
        
        Prob.PositiveCheck = 0 ;           % prior distribution of sita
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
    end
 
    if NE==5501
         Prob.Nr     = 14;   % number of time independent random variables： c        
        Prob.NF     = 6;   % Number of random loads  
        Prob.tmax   = 2;           % Total years
        Prob.dt     = 1;           % time interval
        nKL=fix(Prob.tmax/Prob.dt);
        Prob.nKL =nKL; 
        Prob.L      = 1;
        Prob.Nt     = Prob.tmax/Prob.dt+1;   % total steps of time instance, YY
        Prob.Nx = Prob.Nr+Prob.nKL*Prob.NF;
        for i = 1:Prob.Nx     Prob.Dist{i}= 'norm';      end
        Prob.Dist{Prob.Nr } = 'norm';
%         Prob.mu = [100 200 80  20  200 400 600                800 1000 1200 1400 70  8.75      0.037  zeros(1,Prob.nKL*Prob.NF) ];  % mean : A B C D L1 L2 L3 L4 L5 L6 L Ea Ew
%         Prob.sd = [1   2   0.8 0.2 2   4   6   8   10   12   14   0.7 0.0875    0.00037 ones(1,Prob.nKL*Prob.NF) ];   % standard deviation       
        Prob.mu = [100 200 80  20  200 400 600 800 1000 1200 1400 70  8.75      0.025  zeros(1,Prob.nKL*Prob.NF) ];  % mean : A B C D L1 L2 L3 L4 L5 L6 L Ea Ew
        Prob.sd = [1   2   0.8 0.2 2   4   6   8   10   12   14   0.7 0.0875    0.0025 ones(1,Prob.nKL*Prob.NF) ];   % standard deviation
        Prob.muF = [15 15 15 15 15 15];    % Mean value of dynamic load
        Prob.sdF = [3  3  3  3  3  3];

        covfun1 = @(tao,L) Prob.sdF(1).^2.*exp(-tao.^2./L^2);
        %covfun1 = @(tao,L) Prob.sdF(1).^2.*cos(1*pi*tao/L^2);
        Prob.Covfun = @(tao,L) {covfun1(tao,L),covfun1(tao,L),covfun1(tao,L),covfun1(tao,L),covfun1(tao,L),covfun1(tao,L)};
        [Prob.B, Prob.SigmaB] = KLexpansion(Prob);
        Prob.ZpToFt = @(Zp,Prob) ZpToFt(Zp,Prob);
        Prob.YFunt = @(t) exp(-0.02.*t);% Degenerate formula of mean value of resistance
        
        PL   = @(x,t,F) F(1,:).*(x(:,11)-x(:,5))+ F(2,:).*(x(:,11)-x(:,6))+ F(3,:).*(x(:,11)-x(:,7))+ F(4,:).*(x(:,11)-x(:,8))+ F(5,:).*(x(:,11)-x(:,9))+ F(6,:).*(x(:,11)-x(:,10));
        Ymax = @(x,t,F) ( 0.5* x(:,1).* x(:,2).^2+ x(:,12)./ x(:,13).* x(:,4).* x(:,3).*( x(:,2)+ x(:,4)))./( x(:,1).* x(:,2)+x(:,12)./ x(:,13).* x(:,4).* x(:,3));
        Iw   = @(x,t,F) 1/12*x(:,1).*x(:,2).^3 + x(:,1).*x(:,2).*(Ymax(x,t,F)-0.5*x(:,2)).^2+ 1/12*x(:,12)./ x(:,13).* x(:,3).* x(:,4).^3 + x(:,12)./ x(:,13).* x(:,3).* x(:,4).*(0.5*x(:,4)+x(:,2)-Ymax(x,t,F)).^2 ;
        Sigma3 =@(x,t,F) (PL(x,t,F).*x(:,7)./x(:,11)-F(1,:).*(x(:,7)-x(:,5))- F(2,:).*(x(:,7)-x(:,6))).*Ymax(x,t,F)./ Iw(x,t,F) ;
        Prob.Fungt  = @(x,t,F) x(:,14).*Prob.YFunt(t)  - Sigma3(x,t,F);  % LSF, response function
      
        Prob.dvLocal =  [ 1 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1]};    
        Prob.si_Sita = {[1]};    
        Prob.Domain_x_s = {{[90 110]}};  
        
%         Prob.dvLocal =  [1 2 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1],[1]};    
%         Prob.si_Sita = {[1],[2]};    
%         Prob.Domain_x_s = {{[8e-4 12e-4]},{[ 0.03  0.05]}}; 
       
        Prob.PositiveCheck = 0 ;           % prior distribution of sita
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
    end
    
    if NE==5201
        Prob.Nr  = 3;   % number of random variables： [h  w  E]
        Prob.NF  = 2;   % Number of random loads       
       
        Prob.tmax = 2;           % Total years
        Prob.dt = 1;           % time interval
        nKL=fix(Prob.tmax/Prob.dt);
        Prob.nKL =nKL; 

        Prob.L  = 1;
        Prob.Nt = Prob.tmax/Prob.dt+1;   % total steps of time instance, YY     
        Prob.Nx = Prob.Nr+Prob.nKL*Prob.NF;
        
        Prob.Dist = {'norm','norm','norm'};
        for i = Prob.Nr+1:Prob.Nx
            Prob.Dist{i}= 'norm';
        end
        Prob.mu = [1.5    2    300000 zeros(1,Prob.nKL*Prob.NF) ];  % mean :   w h  E
        Prob.sd = [0.015 0.02  30000  ones(1,Prob.nKL*Prob.NF) ];   % standard deviation  
        Prob.muF = [1000  500];    % Mean value of dynamic load
        Prob.sdF = [200   100];
        
        covfun1 = @(tao,L) Prob.sdF(1).^2.*exp(-tao.^2./L^2);
        covfun2 = @(tao,L) Prob.sdF(2).^2.*exp(-tao.^2./L^2);
        % covfun2 = @(tao,L) Prob.sdF(2).^2.*cos(pi*tao.*L^1);
        Prob.Covfun = @(tao,L) {covfun1(tao,L),covfun2(tao,L)};
        [Prob.B, Prob.SigmaB] = KLexpansion(Prob);
        Prob.ZpToFt = @(Zp,Prob) ZpToFt(Zp,Prob);
%         Prob.YFunt = @(t) exp(-0.1.*t);   % Degenerate formula of resistance
        Prob.Fungt  = @(x,t,F) x(:,3).*exp(-0.1.*t)-(600.*F(1,:)./(x(:,1).*x(:,2).^2)+600.*F(2,:)./(x(:,1).^2.*x(:,2)));  % LSF, response function
         
        % Case 1
        Prob.dvLocal =  [ 2 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1]};    
        Prob.si_Sita = {[1]};    
        Prob.Domain_x_s = {{[1 2.1]}};  
       
        Prob.PositiveCheck = 0 ;           % prior distribution of sita
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
    end
    
    if NE == 41101        
        Prob.Nr = 2;  % number of variables, x
        Prob.NF  = 1;   % Number of random loads       
        Prob.nKL = 0;  % numbers of KL expension
        Prob.tmax = 2;           % Total years
        Prob.dt = 0.08;           % time interval
        nKL=fix(Prob.tmax/Prob.dt);
        Prob.nKL =nKL; 
        Prob.Dist = {'norm','norm'}; % distribution of x
        Prob.Nx = Prob.Nr+Prob.nKL*Prob.NF;
        for i = Prob.Nr+1:Prob.Nx    Prob.Dist{i}= 'norm';     end
        Prob.mu = [10    5    zeros(1,Prob.nKL*Prob.NF) ];  % mean of [x z]
        Prob.sd = [1    0.5   ones(1,Prob.nKL*Prob.NF) ];   % standard deviation  
%         Prob.Para = {[10, 1], [0, 0.5]};
        Prob.muF = [2 ];    % Mean value of stochastic load
        Prob.sdF = [0.2  ];
        
        Prob.Nx = Prob.Nr+Prob.nKL*Prob.NF;
        Prob.Nt = Prob.tmax/Prob.dt+1;   % total steps of time instance, YY
    
       covfun1 = @(tao,L) Prob.sdF(1).^2.*exp(-0.05*tao.^1./L^2);  % stochastic load 
        Prob.Covfun = @(tao,L) {covfun1(tao,L)};
        [Prob.B, Prob.SigmaB] = KLexpansion(Prob);
        Prob.ZpToFt = @(Zp,Prob) ZpToFt(Zp,Prob);
        Prob.Fung  = @(x,t,F) 17-x(:,1)+2*x(:,2).*exp(-0.1.*t)-5*F(1,:);
        Prob.Fungt  = @(x,t,F) 17-x(:,1)+2*x(:,2).*exp(-0.1.*t)-5*F(1,:);  % LSF, x: [N,Nx], t:[1,Nt], F: [NF,Nt]

        % Case 1
        Prob.dvLocal =  [ 1 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1]};    
        Prob.si_Sita = {[1]};    
        Prob.Domain_x_s = {{[8  14]}};  
       
        Prob.PositiveCheck = 0 ;           % prior distribution of sita
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
    end
    
    if NE == 101        
        Prob.Nr = 2;  % number of variables, x
        Prob.NF  = 0;   % Number of random loads       
        Prob.nKL = 0;  % numbers of KL expension
        Prob.tmax = 10;           % Total years
        Prob.dt = 5;           % time interval
        Prob.Dist = {'norm','norm'}; % distribution of x
        Prob.Para = {[0, 1], [0, 1]};
         
        Prob.Nx = Prob.Nr+Prob.nKL*Prob.NF;
        Prob.Nt = Prob.tmax/Prob.dt+1;   % total steps of time instance, YY
        
%         Prob.at0 = @(t) 5;
        Prob.at1 = @(t) -1*exp(0.2.*t);   % Degenerate formula 
%         Prob.at2 = @(t) -1;% 1-0.01*log(t+1);
        Prob.Fungt  = @(x,t,F)  5+ -1*exp(0.1.*t).*x(:,1) + -1.*x(:,2);  % LSF, response function        
        % Case 1
        Prob.dvLocal =  [ 1 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1]};    
        Prob.si_Sita = {[1]};    
        Prob.Domain_x_s = {{[-1  6]}};  
       
        Prob.PositiveCheck = 0 ;           % prior distribution of sita
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
    end
    
    if NE == 40101   % for testing of WLS-BOC , illustration  % 2021-05-13
       Prob.dt=   0.5;
         Prob.tmax=1;
        Prob.Nt=Prob.tmax/Prob.dt+1;
          Prob.Nx = 2;  % number of variables, x
          Prob.Nr = 2;  % number of variables, x
        Prob.Dist = {'norm','norm'}; % distribution of x
        Prob.Para = {[0, 1], [0, 1]};
%         Prob.Dist = {'norm'}; % distribution of x
%          Prob.Para = {[0, 1]};

        Prob.Fung  = @(x) 3+ exp(-0.3*x(:,1))-x(:,2);  % LSF, response function
        Prob.Fungt  = @(x,t,F) 3*exp(-0.5*t)+ exp(-0.3*x(:,1))-x(:,2);  % LSF, response function


        % Case 1
        Prob.dvLocal =  [ 1 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1]};    
        Prob.si_Sita = {[1]};    
        Prob.Domain_x_s = {{[-6  6]}};  
% %         Prob.dvLocal =  [ 1 ];  % location index of the variable related to sita, index_xr 
% %         Prob.si_Para = {[1]};    
% %         Prob.si_Sita = {[1]};    
% %         Prob.Domain_x_s = {{[-5  5]}};  

        
        Prob.PositiveCheck = 0 ;           % prior distribution of sita
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
    end
    
    if NE == 201  % for testing of WLS-BOC , illustration  % 2021-05-13
         Prob.dt=   0.5;
         Prob.tmax=1;
        Prob.Nt=Prob.tmax/Prob.dt+1;
        Prob.Nx = 2;  % number of variables, x
        Prob.Nr=2;
        Prob.Dist = {'norm','norm'}; % distribution of x
        Prob.Para = {[0, 0.5], [0, 0.5]};
        Prob.Funt=@(t) exp(-0.01*t);
        Prob.fun1=@(x) exp(0.04* x(:,1)+7);
        Prob.fun2=@(x) x(:,2).*exp(0.3* x(:,1).^2+5);
        Prob.td=1;
        Prob.Fung  = @(x) Prob.fun1(x).^1-Prob.fun2(x) ;  % LSF, response function
        Prob.Fungt  = @(x,t,F) exp(0.04* x(:,1)+7).^exp(-0.1*t) -1.*exp(0.3* x(:,1).^2+5) 
%         Prob.Fungt  = @(x,t,F) exp(0.04* x(:,1)+7).^exp(-0.5*t-0.5) - x(:,2).*exp(0.3* x(:,1).^2+5) ;

        % opitimization
        % Case 1: tow dimension sita, two mean
        Prob.dvLocal =  [1 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1]};    
        Prob.si_Sita = {[1]};    
        Prob.Domain_x_s = {{[-2  2]}}; 
        Prob.PositiveCheck = 0;           % prior distribution of sita
        
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);

   end
     
    if NE == 11101  %屋架算例,changed, take As and Ac as the design parameters, 2021-09-18
        Prob.dt=   0.25;
        Prob.tmax=0.5;
        Prob.Nt=Prob.tmax/Prob.dt+1;
        Prob.Nx = 6;
        Prob.Nr=6;
        Prob.Dist = {'norm','norm','norm','norm','norm','norm'};
        Prob.mu = [ 9.82e-4  0.04 1e11 2e10   20000 12  ];  % [As Ac Es Ec q l ]
        Prob.sd = [ 9.82e-5 0.004 1e10 2e9    2000  1.2 ];
        % Prob.Fung = @(x)  0.05- x(:,6).*x(:,1).^2./2.*(  3.81./(x(:,3).*x(:,5))+1.13./(x(:,2).*x(:,4))   );
        Prob.Fung = @(x)  0.05*exp(0.5)- x(:,5).*x(:,6).^2./2.*(  3.81./(x(:,2).*x(:,4))+1.13./(x(:,1).*x(:,3))   );
%         Prob.Fungt = @(x,t,F)  0.05*exp(-1*t.^0.25-0.5)- x(:,5).*x(:,6).^2./2.*(  3.81./(x(:,2).*x(:,4))+1.13./(x(:,1).*x(:,3))   );
        Prob.Fungt = @(x,t,F)  0.05*exp(-1*t+0.5)- x(:,5).*x(:,6).^2./2.*(  3.81./(x(:,2).*x(:,4))+1.13./(x(:,1).*x(:,3))   );
        
% Case 1: passed
%         Prob.dvLocal =  [1];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1]};    
%         Prob.si_Sita = {[1]};    
%         Prob.Domain_x_s = {{[8e-4 12e-4]}}; 
        
        % Case 2: passed
        Prob.dvLocal =  [1 2 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1],[1]};    
        Prob.si_Sita = {[1],[2]};    
        Prob.Domain_x_s = {{[8e-4 12e-4]},{[ 0.03  0.05]}}; 

        Prob.fminobj =@(s) 3.1112.*s(:,1)+1.1859.*s(:,2);  
        Prob.Pftol = 0.0001;
        Prob.Htol = 0.1;    % 相对误差
        Prob.Xtol = 0.02;    % 相对误差
        % Prob.S0 = Prob.ub;   % initial point
        
        Prob.Pfmcs = 0.000670;
        Prob.PositiveCheck = 1;           
        Prob = PrepareProb(Prob);
        Prob.output = 'NE2-truss.output';
    end 
    
    if NE == 501  % for testing of WLS-BOC , illustration  % 2021-05-13
         Prob.dt=   0.5;
         Prob.tmax=1;
        Prob.Nt=Prob.tmax/Prob.dt+1;
        Prob.Nx = 2;  % number of variables, x
        Prob.Nr=2;
        Prob.Dist = {'norm','norm'}; % distribution of x
        Prob.Para = {[0, 0.5], [0, 0.5]};
        Prob.Funt=@(t) exp(-0.01*t);
        Prob.fun1=@(x) exp(0.04* x(:,1)+7);
        Prob.fun2=@(x) x(:,2).*exp(0.3* x(:,1).^2+5);
        Prob.td=1;
%         Prob.Fung  = @(x) Prob.fun1(x).^Prob.td-Prob.fun2(x) ;  % LSF, response function
        Prob.Fung  = @(x) Prob.fun1(x).^1-Prob.fun2(x) ;  % LSF, response function
%        Prob.Fungt  = @(x,t,F) exp(0.04* x(:,1)+7).*Prob.td - x(:,2).*exp(0.3* x(:,1).^2+5) ;
%         Prob.Fungt  = @(x,t,F) exp(0.04* x(:,1)+7).^exp(-0.1*t) - x(:,2).*exp(0.3* x(:,1).^2+5) ;
        Prob.Fungt  = @(x,t,F) exp(0.04* x(:,1)+7).^exp(-0.1*t) -1.*exp(0.3* x(:,1).^2+5) 
%         Prob.Fungt  = @(x,t,F) exp(0.04* x(:,1)+7).^exp(-0.5*t-0.5) - x(:,2).*exp(0.3* x(:,1).^2+5) ;

        % opitimization
        % Case 1: tow dimension sita, two mean
        Prob.dvLocal =  [1 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1]};    
        Prob.si_Sita = {[1]};    
        Prob.Domain_x_s = {{[-2  2]}}; 
        Prob.PositiveCheck = 0;           % prior distribution of sita
        
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);

    end
   
    if NE == 401   % for testing of WLS-BOC , illustration  % 2021-05-13
       Prob.dt=   0.5;
         Prob.tmax=1;
        Prob.Nt=Prob.tmax/Prob.dt+1;
          Prob.Nx = 2;  % number of variables, x
          Prob.Nr = 2;  % number of variables, x
        Prob.Dist = {'norm','norm'}; % distribution of x
        Prob.Para = {[0, 1], [0, 1]};
%         Prob.Dist = {'norm'}; % distribution of x
%          Prob.Para = {[0, 1]};

        Prob.Fung  = @(x) 3+ exp(-0.3*x(:,1))-x(:,2);  % LSF, response function
        Prob.Fungt  = @(x,t,F) 3*exp(-1*t)+ exp(-0.3*x(:,1))-x(:,2);  % LSF, response function


        % Case 1
%         Prob.dvLocal =  [ 1 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1]};    
%         Prob.si_Sita = {[1]};    
%         Prob.Domain_x_s = {{[-6  6]}};  
% %         Prob.dvLocal =  [ 1 ];  % location index of the variable related to sita, index_xr 
% %         Prob.si_Para = {[1]};    
% %         Prob.si_Sita = {[1]};    
% %         Prob.Domain_x_s = {{[-5  5]}};  

        % Case 2
%         Prob.dvLocal =  [ 1 ];  % to be continue
%         Prob.si_Para = {[1,2]};    
%         Prob.si_Sita = {[1,2]};    
%         Prob.Domain_x_s = {{[-6  8],[0.1, 1.5]}};  

%         Case 3
        Prob.dvLocal =  [ 1 2];  % good
        Prob.si_Para = {[1],[1]};    
        Prob.si_Sita = {1,2};    
        Prob.Domain_x_s = {{[-5  5]},{[-3, 3]}};  
        
        Prob.PositiveCheck = 0 ;           % prior distribution of sita
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
    end

    if NE == 12401         %汽车前轴结构 log-normal 2021-09-18, used for WLS
        % front axles with lognorm x(5) and x(6): X=lnY~N(mu=P1,simga=P2),
        % mean(Y) = 3.5; CovY = 0.35;
        Prob.dt=   0.5;
         Prob.tmax=1;
        Prob.Nt=Prob.tmax/Prob.dt+1;
        Prob.Nx = 6; % Pf= 0.1369(MCS,1e5)
        Prob.Nr = 6;
        Prob.Dist = {'norm','norm','norm','norm','norm','norm'};
%         Prob.mu = [12  14 65 85 1.2515 1.1302]; % [a t b h M T]=[3.5,3.1] Cov = 0.05 for all
%         Prob.cov = [0.05*ones(1,4)  0.049968  0.049968 ];    %
        Prob.mu = [12  14   65    85   1.2478  1.1264]; % [a t b h M T] Cov = 0.1 for all
        % Prob.cov = [0.1*ones(1,4)  0.09975  0.09975 ];    % 
        Prob.sd = [1.2  1.4  6.5  8.5   0.09975  0.09975] ;
        Wx   = @(x) x(:,1).*(x(:,4)-2.*x(:,2)).^3./(6.*x(:,4))+x(:,3).*(x(:,4).^3-(x(:,4)-2.*x(:,2))).^3./(6.*x(:,4));
        Wp =   @(x) 0.8.*x(:,3).*x(:,2).^2+0.4.*(x(:,1).^3.*(x(:,4)-2.*x(:,2))./x(:,2));
        Sigma = @(x) exp(x(:,5))./Wx(x);    % Note: exp transfromation is used for x(5)
        tao = @(x) exp(x(:,6))./Wp(x);      % Note: exp transfromation is used for x(6)
        Prob.Fung  = @(x) 0.68/1000 - sqrt(Sigma(x).^2+3.*tao(x).^2);
        Prob.Fungt  = @(x,t,F) 0.68/1000*exp(-0.1*t) - sqrt(Sigma(x).^2+3.*tao(x).^2);
        
        % Case 1: general cases
%         Prob.dvLocal    = [ 1   ];         % index of x_r in x; which x is related to sita
%         Prob.si_Para    = {[1]  };         % index of distribution parameters (1, 2 or 3) of x
%         Prob.si_Sita    = {[1]  };         % define the index of sita_i in the vector of Sita
%         Prob.Domain_x_s = {{[10 14]}};     % design region for each sita_i, cell{} for x_r, [] for component
        
        % Case 2: general cases
        Prob.dvLocal    = [ 1,  2 ,  3 ,4  ];         % index of x_r in x; which x is related to sita
        Prob.si_Para    = {[1],     [1],    [1],    [1]  };         % index of distribution parameters (1, 2 or 3) of x
        Prob.si_Sita    = {[1],     [2],    [3],    [4]  };         % define the index of sita_i in the vector of Sita
        Prob.Domain_x_s = {{[10 14]},{[12 16]},{[60, 70]},{[80 90 ]}};     % design region for each sita_i, cell{} for x_r, [] for component
        
        % opitimization
        Prob.output = 'NE12.output';
        Prob.fminobj =@(s) s(:,1).*s(:,4)+2.*s(:,3).*(s(:,2)-s(:,1));
        % Prob.Cons = @(s) (s(:,1)*s(:,4)+2*(s(:,3)-s(:,1)).*s(:,2))-2350;  % V = ah+2(b-a)t in [2154 2880]；
        % Prob.GradCons = @(s) [s(:,4)-2*s(:,2);  2*(s(:,3)-s(:,1)); 2*s(:,2); s(:,1)];
        Prob.Pftol = 0.001;
        Prob.Htol = 50;
        % Prob.XP = [12.9998   14.9998   61.5944   80.1664 1.2478 1.1264];    % D-L Case B 1/2 ： Pf=0.001336   278*10^7  (default) OK
        Prob.Pfmcs =  0.1372;
        % for augmented space (x,sita)
        Prob.Sstyle = 'Unif'; % 'Beta';
        Prob.Spar1 = [ 1 ];
        Prob.Spar2 = [ 1 ];
        
        Prob.PositiveCheck = 1;           
        Prob = PrepareProb(Prob);
    end
    
    if NE == 1401  % 算例 5 A composite beam 更新修改过的
        % A composite beam (X Du: probabilistic uncertainty analysis by  mean-value FOS)
         Prob.dt=   0.5;
         Prob.tmax=0;
        Prob.Nt=Prob.tmax/Prob.dt+1;
        Prob.Nr = 19;
        Prob.Nx = 19;
        %Prob.Dist = {'norm','norm','norm','norm','norm','norm','norm','norm','norm','norm',     'norm','norm','norm','norm','norm','norm','norm','norm','norm','norm'};
        Prob.Dist = {'norm','norm','norm','norm','norm','norm','norm','norm','norm','norm',     'norm','ev','ev','ev','ev','ev','ev','norm','norm'};
        Prob.mu = [100 200 80 20   200 400 600 800 1000 1200  1400   15 15 15 15 15 15     70 8.75];
        Prob.sd = [5 10 4  1       2   4    6  8   10   12    14     3 3 3 3 3 3          0.7  0.0875 ];
        % Prob.cov = [0.01*ones(1,11),0.2*ones(1,6),0.01*ones(1,2)];
        PL   = @(x) x(:,12).*(x(:,11)-x(:,5))+ x(:,13).*(x(:,11)-x(:,6))+ x(:,14).*(x(:,11)-x(:,7))+ x(:,15).*(x(:,11)-x(:,8))+ x(:,16).*(x(:,11)-x(:,9))+ x(:,17).*(x(:,11)-x(:,10));
        Ymax = @(x) ( 0.5* x(:,1).* x(:,2).^2+ x(:,18)./ x(:,19).* x(:,4).* x(:,3).*( x(:,2)+ x(:,4)))./( x(:,1).* x(:,2)+x(:,18)./ x(:,19).* x(:,4).* x(:,3));
        Iw   = @(x) 1/12*x(:,1).*x(:,2).^3 + x(:,1).*x(:,2).*(Ymax(x)-0.5*x(:,2)).^2+ 1/12*x(:,18)./ x(:,19).* x(:,3).* x(:,4).^3 + x(:,18)./ x(:,19).* x(:,3).* x(:,4).*(0.5*x(:,4)+x(:,2)-Ymax(x)).^2 ;
        Sigma2 =@(x) (PL(x).*x(:,6)./x(:,11)-x(:,12).*(x(:,6)-x(:,5))).*Ymax(x)./ Iw(x) ;
        Sigma3 =@(x) (PL(x).*x(:,7)./x(:,11)-x(:,12).*(x(:,7)-x(:,5))- x(:,13).*(x(:,7)-x(:,6))).*Ymax(x)./ Iw(x) ;
        Sigma4 =@(x) (PL(x).*x(:,8)./x(:,11)-x(:,12).*(x(:,8)-x(:,5))- x(:,13).*(x(:,8)-x(:,6))- x(:,14).*(x(:,8)-x(:,7))).*Ymax(x)./ Iw(x) ;
        Sigma5 =@(x) (PL(x).*x(:,9)./x(:,11)-x(:,12).*(x(:,9)-x(:,5))- x(:,13).*(x(:,9)-x(:,6))- x(:,14).*(x(:,9)-x(:,7))- x(:,15).*(x(:,9)-x(:,8))).*Ymax(x)./ Iw(x) ;
          Prob.Fung  = @(x) 19.8/1000 - max( [Sigma2(x), Sigma3(x), Sigma4(x), Sigma5(x)],[],2 );
          Prob.Fungt  = @(x,t,F) 19.8/1000 * exp(-0.1*t) - max( [Sigma2(x), Sigma3(x), Sigma4(x), Sigma5(x)],[],2 );
        %Prob.Fung=@(x) 19.5/1000-max( [Sigma4(x)],[],2 );
        
        % opitimization
        % Case 1: tow dimension sita, two mean
%         Prob.dvLocal =  [1 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1]};    
%         Prob.si_Sita = {[1]};    
%         Prob.Domain_x_s = {{[90  100]}}; 

        % Case 1: passed
        Prob.dvLocal =  [1 2 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1],[1]};    
        Prob.si_Sita = {[1],[2]};    
        Prob.Domain_x_s = {{[80  120]},{[170  250]}}; 
        
        % Case : not good for sigma 
%         Prob.dvLocal =  [1 2 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1],[1 2]};    
%         Prob.si_Sita = {[1],[2 3]};    
%         Prob.Domain_x_s = {{[80  120]},{[180  220],[2 20]}}; 
        
        % Case : not good for sigma 
%         Prob.dvLocal =  [1 12 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1],[1]};    
%         Prob.si_Sita = {[1],[2]};    
%         Prob.Domain_x_s = {{[80  110]},{[12 18]}}; 

        % Case : not good for sigma 
%         Prob.dvLocal =  [14 15 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1],[1]};    
%         Prob.si_Sita = {[1],[2]};    
%         Prob.Domain_x_s = {{[10  25]},{[10 25]}}; 

        % Prob.dvDomain ={[95 105],[190 210], [75 85],[18 22],[],[],[],[],[],[],    [],[],[],[],[],[],[],[],[]} ;
        Prob.fminobj =@(s) s(:,1)*s(:,2);
        Prob.Pftol = 0.0001;
        Prob.Htol = 100;
        Prob.X0 = Prob.mu;
        % Prob.S0 = Prob.lb;
         Prob.XP = [105   210  85 22 200 400 600 800 1000 1200  1400   15 15 15 15 15 15     70 8.75];
        Prob = ParaStat(Prob,'StoP');
        Prob.Pfmcs = 0.010160 ;
        
        Prob.Sstyle = 'Unif'; % 'Beta';
        % Prob.Sstyle={ 'Unif'};
        Prob.Ssd = [ 5/3  10/3];
        Prob.PositiveCheck = 1;           
        Prob = PrepareProb(Prob);
    end


%--------------------------------------------------------------------------------------
%----------------------------------------时变 静态分割线--------------------------
%--------------------------------------------------------------------------------------
   
if NE == 4   % for testing of WLS-BOC , illustration  % 2021-05-13
        Prob.Nx = 2;  % number of variables, x
        Prob.Dist = {'norm','norm'}; % distribution of x
        Prob.Para = {[0, 1], [0, 1]};
        Prob.Fung  = @(x) 3+ exp(-0.3*x(:,1))-x(:,2);  % LSF, response function

        % Case 1
%         Prob.dvLocal =  [ 1 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1]};    
%         Prob.si_Sita = {[1]};    
%         Prob.Domain_x_s = {{[-6  6]}};  

        % Case 2
%         Prob.dvLocal =  [ 1 ];  % to be continue
%         Prob.si_Para = {[1,2]};    
%         Prob.si_Sita = {[1,2]};    
%         Prob.Domain_x_s = {{[-6  8],[0.1, 1.5]}};  

        % Case 3
        Prob.dvLocal =  [ 1 2];  % good
        Prob.si_Para = {[1],[1]};    
        Prob.si_Sita = {1,2};    
        Prob.Domain_x_s = {{[-5  5]},{[-3, 3]}};  
        
        Prob.PositiveCheck = 0 ;           % prior distribution of sita
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
end
   
   if NE == 5   % for testing of WLS-BOC , illustration  % 2021-05-13
        Prob.Nx = 2;  % number of variables, x
        Prob.Dist = {'norm','norm'}; % distribution of x
        Prob.Para = {[0, 0.5], [0, 0.5]};
        Prob.Fung  = @(x) exp(0.04* x(:,1)+7) - x(:,2).*exp(0.3* x(:,1).^2+5) ;  % LSF, response function

        % opitimization
        % Case 1: tow dimension sita, two mean
%         Prob.dvLocal =  [1 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1]};    
%         Prob.si_Sita = {[1]};    
%         Prob.Domain_x_s = {{[-2  2]}}; 
%         Prob.PositiveCheck = 0;           % prior distribution of sita
        
        % Case 2: tow dimension sita, 
        Prob.dvLocal =  [ 1 2 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1],[1]};    
        Prob.si_Sita = {[1],[2]};    
        Prob.Domain_x_s = {{[-2  2]},{[-2, 2]}}; 

        % Case 3: one x 
%         Prob.dvLocal =  [ 1 2 ];  % location index of the variable related to sita, index_xr 
%         Prob.Nd = size(Prob.dvDomain,2); % number of sita
%         Prob.si_Para = {[1],[1]};    
%         Prob.si_Sita = {[1],[2]};    
%         Prob.Domain_x_s = {{[-3  3]},{-1,1}}; 

        % Case 3: all 
%         Prob.dvLocal =  [ 1 2 ];  % location index of the variable related to sita, index_xr 
%         Prob.Nd = size(Prob.dvDomain,2); % number of sita
%         Prob.si_Para = {[1],[1]};    
%         Prob.si_Sita = {[1],[2]};    
%         Prob.Domain_x_s = {{[-3  3]},{-1,1}}; 
        
        Prob.PositiveCheck = 0;           % prior distribution of sita
        Prob.Pftol = 0;             % Probability tolerant
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
   end
   
    if NE == 14   % 算例 5 A composite beam 更新修改过的
        % A composite beam (X Du: probabilistic uncertainty analysis by  mean-value FOS)
        Prob.Nx = 19;
        %Prob.Dist = {'norm','norm','norm','norm','norm','norm','norm','norm','norm','norm',     'norm','norm','norm','norm','norm','norm','norm','norm','norm','norm'};
        Prob.Dist = {'norm','norm','norm','norm','norm','norm','norm','norm','norm','norm',     'norm','ev','ev','ev','ev','ev','ev','norm','norm'};
        Prob.mu = [100 200 80 20   200 400 600 800 1000 1200  1400   15 15 15 15 15 15     70 8.75];
        Prob.sd = [5 10 4  1       2   4    6  8   10   12    14     3 3 3 3 3 3          0.7  0.0875 ];
        % Prob.cov = [0.01*ones(1,11),0.2*ones(1,6),0.01*ones(1,2)];
        PL   = @(x) x(:,12).*(x(:,11)-x(:,5))+ x(:,13).*(x(:,11)-x(:,6))+ x(:,14).*(x(:,11)-x(:,7))+ x(:,15).*(x(:,11)-x(:,8))+ x(:,16).*(x(:,11)-x(:,9))+ x(:,17).*(x(:,11)-x(:,10));
        Ymax = @(x) ( 0.5* x(:,1).* x(:,2).^2+ x(:,18)./ x(:,19).* x(:,4).* x(:,3).*( x(:,2)+ x(:,4)))./( x(:,1).* x(:,2)+x(:,18)./ x(:,19).* x(:,4).* x(:,3));
        Iw   = @(x) 1/12*x(:,1).*x(:,2).^3 + x(:,1).*x(:,2).*(Ymax(x)-0.5*x(:,2)).^2+ 1/12*x(:,18)./ x(:,19).* x(:,3).* x(:,4).^3 + x(:,18)./ x(:,19).* x(:,3).* x(:,4).*(0.5*x(:,4)+x(:,2)-Ymax(x)).^2 ;
        Sigma2 =@(x) (PL(x).*x(:,6)./x(:,11)-x(:,12).*(x(:,6)-x(:,5))).*Ymax(x)./ Iw(x) ;
        Sigma3 =@(x) (PL(x).*x(:,7)./x(:,11)-x(:,12).*(x(:,7)-x(:,5))- x(:,13).*(x(:,7)-x(:,6))).*Ymax(x)./ Iw(x) ;
        Sigma4 =@(x) (PL(x).*x(:,8)./x(:,11)-x(:,12).*(x(:,8)-x(:,5))- x(:,13).*(x(:,8)-x(:,6))- x(:,14).*(x(:,8)-x(:,7))).*Ymax(x)./ Iw(x) ;
        Sigma5 =@(x) (PL(x).*x(:,9)./x(:,11)-x(:,12).*(x(:,9)-x(:,5))- x(:,13).*(x(:,9)-x(:,6))- x(:,14).*(x(:,9)-x(:,7))- x(:,15).*(x(:,9)-x(:,8))).*Ymax(x)./ Iw(x) ;
          Prob.Fung  = @(x) 19.8/1000 - max( [Sigma2(x), Sigma3(x), Sigma4(x), Sigma5(x)],[],2 );
        %Prob.Fung=@(x) 19.5/1000-max( [Sigma4(x)],[],2 );
        
        % opitimization
        % Case 1: tow dimension sita, two mean
%         Prob.dvLocal =  [1 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1]};    
%         Prob.si_Sita = {[1]};    
%         Prob.Domain_x_s = {{[90  100]}}; 

        % Case 1: passed
        Prob.dvLocal =  [1 2 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1],[1]};    
        Prob.si_Sita = {[1],[2]};    
        Prob.Domain_x_s = {{[80  120]},{[170  250]}}; 
        
        % Case : not good for sigma 
%         Prob.dvLocal =  [1 2 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1],[1 2]};    
%         Prob.si_Sita = {[1],[2 3]};    
%         Prob.Domain_x_s = {{[80  120]},{[180  220],[2 20]}}; 
        
        % Case : not good for sigma 
%         Prob.dvLocal =  [1 12 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1],[1]};    
%         Prob.si_Sita = {[1],[2]};    
%         Prob.Domain_x_s = {{[80  110]},{[12 18]}}; 

        % Case : not good for sigma 
%         Prob.dvLocal =  [14 15 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1],[1]};    
%         Prob.si_Sita = {[1],[2]};    
%         Prob.Domain_x_s = {{[10  25]},{[10 25]}}; 

        % Prob.dvDomain ={[95 105],[190 210], [75 85],[18 22],[],[],[],[],[],[],    [],[],[],[],[],[],[],[],[]} ;
        Prob.fminobj =@(s) s(:,1)*s(:,2);
        Prob.Pftol = 0.0001;
        Prob.Htol = 100;
        Prob.X0 = Prob.mu;
        % Prob.S0 = Prob.lb;
         Prob.XP = [105   210  85 22 200 400 600 800 1000 1200  1400   15 15 15 15 15 15     70 8.75];
        Prob = ParaStat(Prob,'StoP');
        Prob.Pfmcs = 0.010160 ;
        
        Prob.Sstyle = 'Unif'; % 'Beta';
        % Prob.Sstyle={ 'Unif'};
        Prob.Ssd = [ 5/3  10/3];
        Prob.PositiveCheck = 1;           
        Prob = PrepareProb(Prob);
    end
    
    if NE == 111  %屋架算例,changed, take As and Ac as the design parameters, 2021-09-18
        Prob.Nx = 6;
        Prob.Dist = {'norm','norm','norm','norm','norm','norm'};
        Prob.mu = [ 9.82e-4  0.04 1e11 2e10   20000 12  ];  % [As Ac Es Ec q l ]
        Prob.sd = [ 9.82e-5 0.004 1e10 2e9    2000  1.2 ];
        % Prob.Fung = @(x)  0.05- x(:,6).*x(:,1).^2./2.*(  3.81./(x(:,3).*x(:,5))+1.13./(x(:,2).*x(:,4))   );
        Prob.Fung = @(x)  0.05- x(:,5).*x(:,6).^2./2.*(  3.81./(x(:,2).*x(:,4))+1.13./(x(:,1).*x(:,3))   );
        
        % Case 1: passed
%         Prob.dvLocal =  [1];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1]};    
%         Prob.si_Sita = {[1]};    
%         Prob.Domain_x_s = {{[8e-4 12e-4]}}; 
        
        % Case 2: passed
        Prob.dvLocal =  [1 2 ];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1],[1]};    
        Prob.si_Sita = {[1],[2]};    
        Prob.Domain_x_s = {{[8e-4 12e-4]},{[ 0.03  0.05]}}; 

        Prob.fminobj =@(s) 3.1112.*s(:,1)+1.1859.*s(:,2);  
        Prob.Pftol = 0.0001;
        Prob.Htol = 0.1;    % 相对误差
        Prob.Xtol = 0.02;    % 相对误差
        % Prob.S0 = Prob.ub;   % initial point
        
        Prob.Pfmcs = 0.000670;
        Prob.PositiveCheck = 1;           
        Prob = PrepareProb(Prob);
        Prob.output = 'NE2-truss.output';
    end   
    
   if NE == 123         % front axles with non-normal x(5) and x(6):  2020-12-6
        Prob.Nx = 6; % Pf= 0.1369(MCS,1e5)  %upadeted:Changed to that same as the OPTModel
        Prob.Dist = {'norm','norm','norm','norm','ev','ev'};% [a t b h M T] 
        Prob.Para = {[13,1.3], [15, 1.5], [65, 6.5], [85, 8.5], [3.65 0.27], [3.24, 0.24]};
        Wx   = @(x) x(:,1).*(x(:,4)-2.*x(:,2))./(6.*x(:,4))+x(:,3).*(x(:,4).^3-(x(:,4)-2.*x(:,2)))./(6.*x(:,4));
        Wp =   @(x) 0.8.*x(:,3).*x(:,2).^2+0.4.*(x(:,1).^3.*(x(:,4)-2.*x(:,2))./x(:,2));
        Sigma = @(x) x(:,5)./Wx(x);
        tao = @(x) x(:,6)./Wp(x);
        Prob.Fung  = @(x) 0.68/1000 - sqrt(Sigma(x).^2+3.*tao(x).^2);
        
        % FPF
        Prob.output = 'NE12.output';
        
        % Case 1: tow dimension sita, one x
%         Prob.dvLocal = [2];         
%         Prob.si_Para = {[1, 2]};    
%         Prob.si_Sita = {[1, 2]};    
%         Prob.Domain_x_s = {{[11 13],[0.8  1.6]}};  % Case B_1:  OK for sigma of x_1 
        % Case 2: tow dimension sita, two x
%         Prob.dvLocal = [1,2];       
%         Prob.si_Para = {[2],[1]};    
%         Prob.si_Sita = {[1], [2]};    
%         Prob.Domain_x_s = {{[0.8  1.6]},{[12 14]}};  

        % Case 3: general cases
%         Prob.dvLocal    = [ 1,                    2   ];         % index of x_r in x; which x is related to sita
%         Prob.si_Para    = {[1,           2],     [1]  };         % index of distribution parameters (1, 2 or 3) of x
%         Prob.si_Sita    = {[1,           2],     [3]  };         % define the index of sita_i in the vector of Sita
%         Prob.Domain_x_s = {{[11, 15],[0.8, 1.6]},{[12,18]}};     % design region for each sita_i, cell{} for x_r, [] for component

%           Case 3: general cases
%         Prob.dvLocal    = [ 1,                       6   ];         % index of x_r in x; which x is related to sita
%         Prob.si_Para    = {[1  ],     [1]  };         % index of distribution parameters (1, 2 or 3) of x
%         Prob.si_Sita    = {[1  ],     [2]  };         % define the index of sita_i in the vector of Sita
%         Prob.Domain_x_s = {{[11, 15]},{[2.8,3.8]}};     % design region for each sita_i, cell{} for x_r, [] for component

%         Prob.dvLocal    = [   6   ];         % index of x_r in x; which x is related to sita
%         Prob.si_Para    = {   [1]  };         % index of distribution parameters (1, 2 or 3) of x
%         Prob.si_Sita    = {   [1]  };         % define the index of sita_i in the vector of Sita
%         Prob.Domain_x_s = {{[2.8,3.8]}};     % design region for each sita_i, cell{} for x_r, [] for component
        
        % Case 3: general cases
        Prob.dvLocal    = [ 1,                    2 ,    6   ];         % index of x_r in x; which x is related to sita
        Prob.si_Para    = {[1,           2],     [1],    [1]  };         % index of distribution parameters (1, 2 or 3) of x
        Prob.si_Sita    = {[1,           2],     [3],    [4]  };         % define the index of sita_i in the vector of Sita
        Prob.Domain_x_s = {{[11, 15],[0.8, 1.6]},{[12,18]},{[2.8,3.8]}};     % design region for each sita_i, cell{} for x_r, [] for component
        

        % Prob.dvLocal = unique(x_Local);         % Location (index) of x_r

        % Prob.X0 = Prob.mu;
        % Prob.S0 = [13 1.3 15 ];
        % Prob.X0 = [13  15	65 85 3.65 3.24];% Pf=  6.5023e-03, P=1 
        % Prob.XP = [10.8 11.8  56.1  78.7 3.55 3.38];% Pf=  6.5023e-03, P=1 
        
        Prob.Pfmcs =  0.1372;
        Prob.Pftol = 0;
        % for augmented space (x,sita)
        Prob.Sstyle = 'Unif'; % 'Beta';
        Prob.PositiveCheck = 1;           
        Prob = PrepareProb(Prob);
        % 
   end 
   
    if NE == 124          %汽车前轴结构 log-normal 2021-09-18, used for WLS
        % front axles with lognorm x(5) and x(6): X=lnY~N(mu=P1,simga=P2),
        % mean(Y) = 3.5; CovY = 0.35;
        Prob.Nx = 6; % Pf= 0.1369(MCS,1e5)
        Prob.Dist = {'norm','norm','norm','norm','norm','norm'};
%         Prob.mu = [12  14 65 85 1.2515 1.1302]; % [a t b h M T]=[3.5,3.1] Cov = 0.05 for all
%         Prob.cov = [0.05*ones(1,4)  0.049968  0.049968 ];    %
        Prob.mu = [12  14   65    85   1.2478  1.1264]; % [a t b h M T] Cov = 0.1 for all
        % Prob.cov = [0.1*ones(1,4)  0.09975  0.09975 ];    % 
        Prob.sd = [1.2  1.4  6.5  8.5   0.09975  0.09975] ;
        Wx   = @(x) x(:,1).*(x(:,4)-2.*x(:,2)).^3./(6.*x(:,4))+x(:,3).*(x(:,4).^3-(x(:,4)-2.*x(:,2))).^3./(6.*x(:,4));
        Wp =   @(x) 0.8.*x(:,3).*x(:,2).^2+0.4.*(x(:,1).^3.*(x(:,4)-2.*x(:,2))./x(:,2));
        Sigma = @(x) exp(x(:,5))./Wx(x);    % Note: exp transfromation is used for x(5)
        tao = @(x) exp(x(:,6))./Wp(x);      % Note: exp transfromation is used for x(6)
        Prob.Fung  = @(x) 0.68/1000 - sqrt(Sigma(x).^2+3.*tao(x).^2);
        
        % Case 1: general cases
%         Prob.dvLocal    = [ 1   ];         % index of x_r in x; which x is related to sita
%         Prob.si_Para    = {[1]  };         % index of distribution parameters (1, 2 or 3) of x
%         Prob.si_Sita    = {[1]  };         % define the index of sita_i in the vector of Sita
%         Prob.Domain_x_s = {{[10 14]}};     % design region for each sita_i, cell{} for x_r, [] for component
        
        % Case 2: general cases
        Prob.dvLocal    = [ 1,  2 ,  3 ,4  ];         % index of x_r in x; which x is related to sita
        Prob.si_Para    = {[1],     [1],    [1],    [1]  };         % index of distribution parameters (1, 2 or 3) of x
        Prob.si_Sita    = {[1],     [2],    [3],    [4]  };         % define the index of sita_i in the vector of Sita
        Prob.Domain_x_s = {{[10 14]},{[12 16]},{[60, 70]},{[80 90 ]}};     % design region for each sita_i, cell{} for x_r, [] for component
        
        % opitimization
        Prob.output = 'NE12.output';
        Prob.fminobj =@(s) s(:,1).*s(:,4)+2.*s(:,3).*(s(:,2)-s(:,1));
        % Prob.Cons = @(s) (s(:,1)*s(:,4)+2*(s(:,3)-s(:,1)).*s(:,2))-2350;  % V = ah+2(b-a)t in [2154 2880]；
        % Prob.GradCons = @(s) [s(:,4)-2*s(:,2);  2*(s(:,3)-s(:,1)); 2*s(:,2); s(:,1)];
        Prob.Pftol = 0.001;
        Prob.Htol = 50;
        % Prob.XP = [12.9998   14.9998   61.5944   80.1664 1.2478 1.1264];    % D-L Case B 1/2 ： Pf=0.001336   278*10^7  (default) OK
        Prob.Pfmcs =  0.1372;
        % for augmented space (x,sita)
        Prob.Sstyle = 'Unif'; % 'Beta';
        Prob.Spar1 = [ 1 ];
        Prob.Spar2 = [ 1 ];
        
        Prob.PositiveCheck = 1;           
        Prob = PrepareProb(Prob);
    end
   
   if NE == 1005   % Thermal Stress Analysis of Jet Engine Turbine Blade  % 2021-07-06
        Prob.Nx = 8;  % number of variables, x
        Prob.Dist = {'norm','norm','norm','norm','norm','norm','norm','norm'}; % [E CTE nu p1 p2 kapp Tin Tout]
        Prob.Para = {[225, 23], [13, 1.3], [0.27, 0.027], [500,50],[450,45], [11.5,1.15], [150,15], [1000,100]};
        % Prob.Para = {[230, 11], [13, 0.6], [0.27, 0.0135], [500,50],[450,45], [11.5,0.575], [150,15], [1000,100]};
        Prob.Fung  = @(x) 15E8-TurbineBlade_LSF(x) ;  % LSF, response function

        % opitimization
        % Case 1: tow dimension sita, two mean
%         Prob.dvLocal =  [1 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1]};    
%         Prob.si_Sita = {[1]};    
%         Prob.Domain_x_s = {{[200  250]}}; 
        
        % Case 2: tow dimension sita, two sigma
%         Prob.dvLocal =  [ 1 2 3 4 5 6 7 8 ];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1],[1],[1],[1],[1],[1],[1],[1]};    
%         Prob.si_Sita = {[1],[2],[3],[4],[5],[6],[7],[8]};    
%         Prob.Domain_x_s = {{[200  250]},{[10  15]},{[0.20  0.35]},{[2  10]},{[2  10]},{[8  15]},{[110  195]},{[700  1300]}}; 

        % Case 3: one x 
%         Prob.dvLocal =  [ 1 2 7 8];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1],[1],[1],[1]};    
%         Prob.si_Sita = {[1],[2],[3],[4]};    
%         Prob.Domain_x_s = {{[200  250]},{[10  15]},{[120  160]},{[700  1300]}}; 

        % Case 4:  
        Prob.dvLocal =  [ 1 2];  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1],[1]};    
        Prob.si_Sita = {[1],[2]};    
        % Prob.Domain_x_s = {{[190  280]},{[10  25]}}; 
        Prob.Domain_x_s = {{[190  300]},{[10  25]}}; 
        
        Prob.Pftol = 0;             % Probability tolerant
        Prob.PositiveCheck = 1;           % prior distribution of sita
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob = PrepareProb(Prob);
   end

   if NE == 2001 % For AuIS 
        % [c1,c2,c3,c4,c5,c6, m1,m2,m3,m4,m5,m6, k1,k2,k3,k4,k5,k6, e1,e2,e3,e4,e5,e6 ]
        Prob.S  = 0.05; % m2/s3
        Prob.dt = 0.05; % dT =nt*dt
        Prob.dT = 10;
        Prob.Hi = 0.015*[5.49, 3.81, 3.81, 3.81, 3.81, 3.81]; % interstory drift ratio 2%*Hi
        Prob.threshold   = Prob.Hi;
        Prob.m = 6;
        Prob.Dof = Prob.m;

        Prob.Ns = Prob.Dof*4;
        Prob.Nz = Prob.dT/Prob.dt +1;
        Prob.nx = Prob.Ns + Prob.Nz;
        M = [283, 263, 256, 255, 247, 215 ]; % *10^3
        K = [150, 367, 246, 246, 175, 175 ];
        C = 2*ones(1,6);    % *10^6
        E = 0.05*ones(1,6);
        Prob.mu = [M,K,C,E,                   zeros(1,Prob.Nz)];   
        Prob.sd = [0.1*M,0.1*K,0.1*C,0.1*E,   ones(1,Prob.Nz) ];
        % Prob.cov = [0.1*ones(1,Prob.m ), 0.1*ones(1,Prob.m ), 0.2*ones(1,Prob.m )]; 
        for i=1:Prob.nx
            Prob.Dist{i}= 'norm';
            Prob.Para{i}= [Prob.mu(i), Prob.sd(i)];
        end
        Prob.Fung = @(x) SixStoryWhiteNoise(x);
        % Prob.Fung_AuIS = @(x) SixStory_AuIS(x);
        Prob.output = 'NE2001.output';

        % Case 3: general cases
        Prob.dvLocal    = [ 7, ];         % index of x_r in x; which x is related to sita
        Prob.si_Para    = { [1],  };         % index of distribution parameters (1, 2 or 3) of x
        Prob.si_Sita    = { [1],  };         % define the index of sita_i in the vector of Sita
        Prob.Domain_x_s = {{[100, 200]}};     % design region for each sita_i, cell{} for x_r, [] for component

        Prob.XP = Prob.mu;
        Prob.Sstyle = 'Unif'; % 'Beta';
        Prob = PrepareProb(Prob);
        Prob.xuLocal = [];       
  end
 
   if NE ==2020 % Benchmark Problem 2,Case 2, For AuIS method,and RBDO paper Ex.3  
        % Ref.Ching jianye CMAME (2008)198 52-71
        % [ m1,...,m10, k1,...,k10,c1,...,c10 ]
        % upadated: 2021-08-01
        Prob.Dof   = 10;
        Prob.Ns    = 30; %3*Prob.Dof ;   %number of structural random variables
        Prob.dt    = 0.05;  % Discrete time step
        Prob.dT    = 20;  %0.005, Discrete time step
        Prob.Nz    = Prob.dT/Prob.dt +1;   %30sec and 3001, number of random variables associated with excitation
        Prob.S     = 0.01;
        Prob.Nx    = Prob.Ns + Prob.Nz;
        Prob.EPS     = 0.0; % Parameter of Duffing oscillator
        Prob.threshold   = [0.06 6*ones(1,9)]; % 0.03*ones(1,Prob.Dof);    %threshold level
        M = 10*ones(1,10);        % Note 1e3 is incorporated in to BenchmarkProblem2
        K = 40*ones(1,10); % [40*ones(1,3),36*ones(1,3),32*ones(1,4)];    % 1e6   
        C = 0.04*ones(1,10);  
        Prob.mu = [M,K,C,               zeros(1,Prob.Nz)];   
        Prob.sd = [0.05*M,0.10*K,0.20*C,ones(1,Prob.Nz) ];
        % Prob.cov = [0.05*ones(1,10), 0.10*ones(1,10), 0.20*ones(1,10)];  
        for i=1:Prob.Nx
            Prob.Dist{i}= 'norm';
            Prob.Para{i}= [Prob.mu(i), Prob.sd(i)];
        end
        Prob.output = 'NE20BenchmarkProblem-2.txt';
        % Prob.Fung = @(x) BenchmarkProblem2_AuIS(x);  % probabilistic constraint: P[g(x)<0]<Pftol
        
        % 2021-08-02 one dimension
        Prob.dvLocal =  11;  % location index of the variable related to sita, index_xr 
        Prob.si_Para = {[1]};    
        Prob.si_Sita = {[1]};    
        Prob.Domain_x_s = {  {[20  60]} }; 
        
        % 2021-08-04  two dimensions
%         Prob.dvLocal =  [ 1 11];  % location index of the variable related to sita, index_xr 
%         Prob.si_Para = {[1],[1]};    
%         Prob.si_Sita = {[1],[2]};    
%         Prob.Domain_x_s = {{[5  15]},{[20  60]}}; 
        
        Prob.Pftol = 0;             % Probability tolerant
        Prob.PositiveCheck = 1;         % prior distribution of sita
        Prob.Sstyle = 'Unif';           % prior distribution of sita
        Prob.xuLocal = [];              % 2021-08-12    
        Prob = PrepareProb(Prob);
        
        % Prob.S0 = Prob.lb;                    % Case A
        % Prob.S0 = (Prob.lb + Prob.ub)/2;      % Case B
         Prob.S0 = Prob.ub;                    % Case C
   end  
   
   function  [B, SigmaB]= KLexpansion(Prob)
% K-L展开 (Karhunen-Loeve expansion)
% return B: [namda(i)*Phi(i)]: Prob.Nt* Prob.nKL*Prob.NF
% update: 2020-09-23

    % global Prob
    nKL = Prob.nKL;
       Prob.L  = 1;
    
    t  = 0: Prob.dt : Prob.tmax;    % 行
    C0 = Prob.Covfun(t,Prob.L);     % cell = { [1,nt] }
    
    B=[];SigmaB=[];
    for i=1:Prob.NF % compute the KL seriel for each F(t)
        C = toeplitz(C0{i});                    %generate discrete correlation matrix
        [Psi,Lambda]    = eig(C);                       %calculate eigenvalues/eigenvectors
        [~,ind]         = sort(diag(Lambda),'descend'); %sort values
        Psi             = Psi(:,ind);                   %sort values
        Lambda          = Lambda(ind,ind);              %sort values
        % nKL       = find(cumsum(diag(Lambda))/sum(diag(Lambda))>perc_KL,1,'first');     %determine number of KL terms to be retained
        B0          = Psi(:,1:nKL)*(Lambda(1:nKL,1:nKL))^0.5;   %matrix for generating random samples
        sumB        = sqrt(sum(B0.^2,2));
        
        B = [B,B0];
        SigmaB = [SigmaB, sumB];
    end
   end

function  [Ft]= ZpToFt(Zp,Prob)
% return Ft: [ ] 
% update: 2020-09-23

    [N,~] = size(Zp);
    if N~=1      error('Zp must be inputted as one line! (@ZpToFt) ');   end
    
    % for j = 1:N
       for i=1:Prob.Nt
            B_Ft0 = Prob.B(i,:).*Zp;
            for nF = 1:Prob.NF  % j=1:Prob.NF 
                % n = Prob.nKL*(j-1)+1;  s = Prob.nKL*j;  Z_Ft(i,j) = sum(B_Ft0(:,n:s),2);
                Z_Ft(i,nF) = sum(B_Ft0(:, Prob.nKL*(nF-1)+1 :Prob.nKL*nF),2);
            end
       end
       Ft      = Prob.muF + Z_Ft;   % Ft: [Nt*NF] 
       Ft      = Ft' ;    % Ft: [NF,Nt] 
    end

end
        
        