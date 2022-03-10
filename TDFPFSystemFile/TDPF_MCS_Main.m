function [Pf,CovPf,unit_CovPf] =  TDPF_MCS_Main(Prob,N,ts)
% MCS for Time-dependent Reliability Analysis,
% ts : input time instant sequence 
% updated: 2020-09-21  Version 1, obtain TDPF vs tl function at one time
% updated: 2021-09-01  obtain TDPF vs t function at any time instant
%input ts:[1,Nt]
%output [1,Nt]

tic
    methodinfo = 'MCS for Pf(t)';
    % global Prob  % NE = 53; N = 10000; 
    % disp(nargin); disp(varargin);   % if nargin >= 3     flag = varargin{3};  end
    %Prob = TDProb(NE);
    Nt = length(ts);
    
    tsall = 0:Prob.dt:Prob.tmax;
    IFt = 1:Prob.Nt;
    IF  = IFt(ismember(ts,tsall)==1) ;   % 
    
    xmcs = MCS_Sampling(N,Prob);      % MCS
    X =  xmcs(:,1:Prob.Nr);
%     if isfield(Prob,'NF') && Prob.NF>0,     Z =  xmcs(:,Prob.Nr+1:end); end
    
    % compute the time responses for every X sample, Nt values respect to each Xj 
    gxt = zeros(N,Nt); %Nf = zeros(1,Nt);
%     parfor j = 1:N  % each sample of Z
   parfor j = 1:N  % each sample of Z
        if isfield(Prob,'NF') && Prob.NF>0     
            Z =  xmcs(j,Prob.Nr+1:end); 
            Ft{j} = Prob.ZpToFt(Z,Prob);       % Ft: [NF,Nt]
            Ftk = Ft{j}(:,IF);
        else
            Ftk = 0;
        end
        gxt(j,:) = Prob.Fungt(X(j,:),ts,Ftk);     % gxt: [ N*Nt]
        %gxt(j,:) = Fungx(X(j,:),ts',Ftk');     % gxt: [ N*Nt]
    end
    
    % Time Pf function with regard to tk
    for k = 1:length(ts)
        % compute the time Pf over [0,tk]
        gt      = min(gxt(:,1:k),[],2);            % gt: [N*1](:,k)
        Ifail   = find(gt<=0);
        Nf      = length(Ifail);
        Pf(k)   = Nf/N;
        CovPf(k)    = sqrt((1-Pf(k))/(N-1)/Pf(k)); 
        unit_CovPf  = CovPf*sqrt(N);
    end
    
    % return
    Pf;
toc

%%%%%%%%%%%%%%%
function  XPoint = MCS_Sampling(N,Prob)
% 产生服从特定分布的样本，分布类型由Prob中定义，注意：有些分布只接受一个分布参数 
% date： 2008-12-30
    Nx = Prob.Nx ;
    Dist = Prob.Dist ;
    Para = Prob.Para ;
    for i = 1:Nx    
        % 变量的分布参数个数
        ns = length(Para{i});
        if ns==1
            x(:,i) = random(Dist{i}, Para{i}(1),[N,1]); 
        end
        if ns==2
            x(:,i) = random(Dist{i}, Para{i}(1),Para{i}(2),[N,1]); 
        end
        if ns==3 % 针对gev 分布
            x(:,i) = random(Dist{i}, Para{i}(3), Para{i}(2),Para{i}(1),[N,1]); 
        end
        if ns>=4
            error('you must check the distribution!' ); 
        end
    end
    XPoint = x; 