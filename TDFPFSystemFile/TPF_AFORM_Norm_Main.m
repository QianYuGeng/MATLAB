function [XP,Pf,beta,Ncall] = TPF_AFORM_Norm_Main(Probk,ts)
% [XP,Pf,beta,Ncall] = TPF_AFORM_Norm_Main(Probs,0)
% AFORM for TPf and design point 
% updated: 2021-09-03

    % clear
    methodinfo = 'AFOSM for Pf(t)';
    % NE = 411;  
   %Probk = TDProb(NE);
    % ts = Probk.tmax;
    Nx = Probk.Nx; Dist = Probk.Dist ; 
    %Fungx = Probk.Fungx;
    
    global Ncall
    Ncall = 0;
    
    [Nt,~] = size(ts);
    
    % 设定求解范围
    BetaInit = [0.1, 5];

    MU = Probk.mu;
    % if isfield(Probk, 'cov'),  SIGMA = Probk.cov.*Probk.mu; end
    if isfield(Probk, 'sd'), SIGMA =  Probk.sd; end
    for i=1:Nx
        if MU(i) == 0
           XP(i) = MU(i)+0.01;
        else
           XP(i) = MU(i);
        end
    end

    % XP = MU+1e-3
    beta = BetaInit(1);
    beta0 = BetaInit(2);
    options=optimset('Jacobian','off','MaxIter',50,'Tolfun',1e-8);
    flaginter = 0;
    %迭代求设计点，直到前后两次可靠性指标差值满足一定的要求
    for k = 1:Nt
        tk = ts(k);
        while abs(beta-beta0)>=1e-3 && flaginter<70
            beta0  = beta;
            Gradxi = Gradgx(XP,tk,Probk);                                         %计算各变量偏导数在迭代出的点处的值
            lamda  = -(Gradxi.*SIGMA)./sqrt(sum((Gradxi.^2).*(SIGMA.^2)));         %计算lamda值
            beta   = fsolve(@(x) FunBeta(x,tk,lamda,Probk),beta0,options);
            XP     = MU+SIGMA.*lamda*beta;       %带入lamda后以可靠度指标beida表示的xi*的表达式
            flaginter = flaginter + 1;
        end
        global delta_beta
        delta_beta=beta-beta0;

        %失效概率
        Pf(k) = normcdf(-beta);
        gx = Fungxt(XP,tk,Probk);
        % [Dxi] = difference(XP);                                              %各变量偏导数在设计点处的值
        c0 = gx-sum(Gradxi.*XP);
        mug = c0+sum(Gradxi.*MU);
        sigmag=sum((Gradxi.*SIGMA).^2);
        beta2=mug/sqrt(sigmag);
        pf_mu=-Gradxi/sqrt(sigmag)*normpdf(beta2) ;                                      %灵敏度计算
        pf_sigma=Gradxi.^2.*SIGMA*mug/sqrt(sigmag)^3*normpdf(beta2);
    end
    
    fprintf('\nPf=%12.4e\n',Pf(k));

    % fid=fopen('AFORMresults.txt','at+');
    % outputing(fid,flaginter,Pf,beta,XP,gx,pf_mu,pf_sigma)
    % fclose(fid)
    % open('AFORMresults.txt');
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function    outputing(fid,flaginter,Pf,beta,XP,gxp,pf_mu,pf_sigma)
    fprintf(fid,'\nNo. of iteration: %d',flaginter);
    if flaginter==1
        fprintf(fid,'\nThere may be something wrong with inter num is 1!');
    end
    fprintf(fid,'\nPf=%12.4e',Pf);
    fprintf(fid,'\nbeta= %12.4e',beta);
    fprintf(fid,'\ng(XP)= %12.4e',gxp);
    fprintf(fid,'\nXP = ');
    fprintf(fid,'%12.4e\t',XP);
    fprintf(fid,'\nSensitivity  of Mean and SD are\n');
    fprintf(fid,'%12.4e\t',pf_mu);
    fprintf(fid,'\n');
    fprintf(fid,'%12.4e\t',pf_sigma);
    fprintf(fid,'\n');

function gbeta = FunBeta(beta,tk,cosx,Probk)
% 关于beta 的方程，用于求解
    global Ncall

    MU =  Probk.mu; 
    % if isfield(Prob, 'cov'),  SIGMA = Prob.cov.*Prob.mu; end
    if isfield(Probk, 'sd'), SIGMA =  Probk.sd; end

    x = MU+beta.*cosx.*SIGMA;
    gbeta = Fungxt(x,tk,Probk) ;   % updated: 2021-09-03

    Ncall = Ncall +1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Grad = Gradgx(XP,tk,Probk)     
% 计算指定点处的导数
% Date : 2009-4-30
    global Ncall
    % 求解导数
    NX = length(XP);
    Grad = zeros(size(XP));
    for i=1:NX
        if XP(i) == 0
            h = XP(i)+0.01;
        else
            h = XP(i)*0.01;
        end
        delta = zeros(size(XP));
        delta(i) = h;
        x = XP+delta/2;          g1 =  Fungxt(x,tk,Probk);  % updated: 2021-09-03
        x = XP-delta/2;          g2 =  Fungxt(x,tk,Probk);
        x = XP+delta;            g3 =  Fungxt(x,tk,Probk);
        x = XP-delta;            g4 =  Fungxt(x,tk,Probk);
        Grad(i) = 4.0/(3*h)*(g1-g2)-1.0/(6.0*h)*(g3-g4);
        %h
        %g1
        
        Ncall = Ncall + 2;
    end
          
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gx = Fungxt(x,tk,Probk)
% imput x = [x_s, x_Z ] 
% updated: 2021-09-03

    ts  = 0: Probk.dt : Probk.tmax; 
    IF = find(ts==tk);  % 指定时间处的LSF

    [N,~] =  size(x);
    X =  x(:,1:Probk.Nr); Z = [] ; Ftk = [];

    gx = [];
    for j = 1:N
        if  isfield(Probk,'NF') && Probk.NF>0     
            Z =  x(j,Probk.Nr+1:end); 
            Ft{j} = Probk.ZpToFt(Z(j,:),Probk);       % Ft: [NF,Nt]
%           Ftk = Ft{j}(IF,:);%F维数按列
            Ftk = Ft{j}(:,IF);%F维数按行
        end
        gx(j,:) = Probk.Fungt(X,tk,Ftk) ;%F维数按列

    end
    
    
    
