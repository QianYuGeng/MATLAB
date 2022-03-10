function [ Xfail_dv, fzhz_no, PfMCS, CovPfMCS] = BenchProb2_AuIS(Probs,N)
% [xjIS, PfjIS, PfMCS, CovPfMCS] = BenchProb2_AuIS(Probs,N)
% First excursion probability by Au importance sampling for Benchmark problem 2 (Linear dynamic) 
% input: the relization of structure variables x
% updated: 2020-10-29
% updated: 2021-01-08
% updated: 2021-05-12, change the check x<0, number to N

     % clear;  NE = 2021;     Probs = OPTModel(NE);   N = 100;  
     warning('off');
     Hmu = Probs.mu;  Hsd = Probs.sd; % Probs.cov.*XP; % note that not Prob.XP any more  
     fprintf('The sampling center : ');fprintf('%6.4f\t',Hmu(Probs.dvLocal)); fprintf('\n');
     
    % detail of the structure and excitation
    dt =  Probs.dt; % 0.1; % dT =nt*Dt
    S  =  Probs.S;  % 0.08; % m2/s3
    Nz =  Probs.Nz;
    Dof = Probs.Dof;
    Idv = Probs.dvLocal;
     
    tk = 0:dt:dt*Nz;   % tk = dt:dt:dt*Nz; wrong not from 0
    
    tic
    % rng(96)
    
    % Sampling with Hmu for structural variables x, N=N;
    Nxs = length(Hmu);
    x   = zeros(N,Nxs);
    nsample = 1;
    while nsample <= N   % updated : 2021-05-12
        x(nsample,:) = normrnd(Hmu,Hsd);

        % delete all the minus value
        Xs = x(nsample,1:3*Dof);
        Ifsum = sum(Xs<0,2);
        % x = x(Ifsum==0,:);
        % [N,~] = size(x);
        
        if Ifsum==0
            nsample = nsample+1;
        end
    end
    
    % construct H(x) for each xi i =1,...,N
    for p=1:N
        if Probs.Ns == 3*Dof % asign the structural variables
            M   = diag(x(p,1:Dof)*1e3);  % kg
            ki  = x(p,Dof+1:2*Dof)*1e6;    %% N/m  
            etai= x(p,2*Dof+1:3*Dof);
        elseif  Probs.Ns == 0 % deterministic structure
            M   = diag(Probs.mu(1:Dof)*1e3);  % kg
            ki  = Probs.mu(Dof+1:2*Dof)*1e6;    %% N/m  
            etai = Probs.mu(2*Dof+1:3*Dof);
        else
            fprintf('check the number of the variabls.');
        end
        
        % K
        for i=1:Dof
            for j=1:Dof
                if i==j && j<Dof,     K(i,j) = ki(i)+ki(i+1);
                elseif i==j && j==Dof,  K(i,j) = ki(i);
                end
                if i>1 && j==i-1,  K(i,j) = -ki(i);  end
                if i<Dof && j==i+1,   K(i,j) = -ki(j);   end
            end
        end

        %[Q,W]= eig(K,M); % A*V=B*V*D: [V,D] = eig(A,B) KQ=w^2MQ    Q'*M*Q = E
        %W=sqrt(W);
        
        % C
        for i =1:Dof
            ci(i) = 2*etai(i)*sqrt(M(i,i)*ki(i));
        end
        for i=1:Dof
            for j=1:Dof
                if i==j && j<Dof,     C(i,j) = ci(i)+ci(i+1);
                elseif i==j && j==Dof,  C(i,j) = ci(i);
                end
                if i>1 && j==i-1,  C(i,j) = -ci(i);  end
                if i<Dof && j==i+1,   C(i,j) = -ci(j);   end
            end
        end
        
        % impluse response of the filter
        ome_1g  = 15;       %rad/s
        z_1g    = 0.8; 
        ome_2g  = 0.3;      %rad/s
        z_2g    = 0.995;

        Af   = [...
           0              1                 0           0;
        -ome_1g^2   -2*z_1g*ome_1g          0           0;
           0              0                 0           1;
        ome_1g^2    2*z_1g*ome_1g     -ome_2g^2   -2*z_2g*ome_2g
        ];
        Bf   = [0; 1; 0; 0];
        Cf   = [ome_1g^2 2*z_1g*ome_1g -ome_2g^2 -2*z_2g*ome_2g];
        Df   = 0;

        sysf = ss(Af,Bf,Cf,Df); % put the system in state space formulation 

        % impulse response of filter
        afilter = impulse(sysf,tk)';  % a   = lsim(sys,w,t,x0,'zoh')';

        % Time-Domain method: solving state-space modal
        Ass = [zeros(Dof),eye(Dof); -M^(-1)*K, -M^(-1)*C];
        Bss = [zeros(Dof,1); ones(Dof,1)];
        Css = [eye(Dof),zeros(Dof);zeros(Dof),zeros(Dof)];%eye(2*Dof);
        %Css =  eye(2*Dof);
        Dss = 0;
        sysY = ss(Ass,Bss,Css,Dss);

        % impulse response of augmented system
        Yimp =  lsim(sysY,afilter,tk,zeros(2*Dof,1),'zoh')';
        
        % Yoimp = (impulse(sysY,tk))'; % impulse response of orignal system
        % for i=1:Dof     Yimp(i,:) = conv(Yoimp(i,:),afilter);       end
        
        % impulse response of interstory drift ratio
        Ydrift = zeros(Dof,Nz);
        Ydrift(1,:) = Yimp(1,1:Nz);  % first floor
        for i =2:Dof                
            Ydrift(i,:) = Yimp(i,1:Nz)-Yimp(i-1,1:Nz);
        end
        
        et   =  [(0:dt:2)*0.5,ones(1,8/dt-1),exp(-0.1*(0:dt:Probs.dT-10))]; % envelop function        
        
        % calculate by AuIS
        Pfj(p,:) = IS_Au_LinearSys(Probs,Ydrift,et);   % fz/hz
        fprintf('\nPfj(fz/hz) = %12.4e\n',Pfj);

    end
    
    % calculate for MCS
     PfMCS = mean(Pfj);% this is the Pf calcultate by MCS when xp and Hsd
    VarPfMCS = (mean(Pfj.^2)-PfMCS^2)/N;
    CovPfMCS = sqrt(VarPfMCS)/PfMCS;
    fprintf(1,'\nPfmcsxp(Cov) = %12.4e( %12.4e) N=%d\n',PfMCS,CovPfMCS,N);

    % save for GIS calculating    
    Xfail_dv = x(:,1:3*Dof); % xjIS = x(:,Idv); 
    fzhz_no  = Pfj;     % PfjIS = Pfj;
    % fprintf(1,'\nPfiIS = %12.4e\n',PfjIS);
    toc

 function [Pfj] = IS_Au_LinearSys(Probs,Ydrift,et)
% Au IS for Z. 
% updated : 2021-01-09

        dt =  Probs.dt; % 0.1; % dT =nt*Dt
        S  =  Probs.S;  % 0.08; % m2/s3
        Nz =  Probs.Nz;
        m  =  Probs.Dof;
        bik = diag(Probs.threshold)*ones(m,Nz); % interstory drift ratio 0.2%*Hi

        % impulse response 
        giks = cell(1,m);
        for i=1:m
            giks{i} = zeros(Nz); 
            for k=1:Nz
                for s=1:k
                    giks{i}(k,s) = Ydrift(i,k-s+1);
                end
            end
        end
                
        
        % Mu_g and Sigma_g of each LSF
        % Sigma2ik = zeros(m,2*Nz-1);
        % for i=1:m  % Eq. 11 Au
        %    Sigma2ik(i,:) = conv(Ydrift(i,:).^2,ones(size(Ydrift(i,:))))*2*pi*S*dt;
        % end

        Sigma2ik = [];
        for i=1:m % First story % Dof
            Sigma2ik(i,:) = conv(Ydrift(i,:).^2,et.^2)*2*pi*S*dt;
            for j = 1:length(Sigma2ik(i,:))
                if Sigma2ik(i,j)==0, Sigma2ik(i,j)= 1E-10; end
            end
        end

        % Sigma2ik = cumsum(Ydrift.^2,2)*2*pi*S*dt;  % updated 2020-11-1
                 
        Sigmaik = sqrt(Sigma2ik(:,1:Nz));
        Betaik  = bik./Sigmaik;
        PFik    = 2*normcdf(-Betaik);
        % for case beta>15 which result in zero in f_beta, assign a small value to it 2020-10-25 Yuan
        if all(PFik<1E-50)
            PFik = 1E-50*ones(size(PFik)); fprintf('All beta are very small: %e \n',PFik(end));
        end
                    
                  
        % The design point
        ZXPik   = cell(m,Nz);       UXP     = cell(m,Nz);
        for i=1:m
            for k=1:Nz
                % design pionts
                ZXPik{i,k}=zeros(1,Nz);
                for s=1:k
                    ZXPik{i,k}(s) = et(s)*sqrt(2*pi*S*dt)*giks{i}(k,s)/Sigma2ik(i,k)*bik(i,k);
                    UXP{i,k} = ZXPik{i,k}/Betaik(i,k);
                end
            end
        end

        PFsum = sum(sum(PFik,2),1);
        Wik = PFik/PFsum;
    
               
        % caculate the cdf of wik for IK sampling
        CdfWik = cumsum(reshape(Wik',1,[])); % updated 2020-11-1
        %plot([1:j],CdfWik);
        if isnan(CdfWik(end)), CdfWik(end)=1; fprintf('Cdf of wi not equal to 1 : %e \n',CdfWik(end)); end

        % choose IK
        u=unifrnd(0,1,1,1);
        Ci = find(CdfWik>u, 1 );
        Iik = ceil(Ci/Nz); % determin the number of rows 
        Kik = Ci- (Iik-1)*Nz; % determin the number of colums 
        % tranform
        U = unifrnd(0,1,1,2);
        if Betaik(Iik,Kik)<0            % updated: 2020-10-24
            alpa = norminv(U(1)+(1-U(1)).*normcdf(Betaik(Iik,Kik)));
        else
            alpa = -norminv((1-U(1)).*normcdf(-Betaik(Iik,Kik)));
        end
        if isinf(alpa), alpa= Betaik(Iik,Kik)+1; fprintf('alpa is too big! : %e \n',alpa); end
             
        
        Z = normrnd(0,1,1,Nz);
        if U(2)<=1/2 
            Zr = Z+(alpa-Z*UXP{Iik,Kik}')*UXP{Iik,Kik};
        else
            Zr = -Z-(alpa-Z*UXP{Iik,Kik}')*UXP{Iik,Kik};
        end
        
                 
        % compute the IFik(Z)
        Y = zeros(m,Nz);
        for i=1:m
            for k=1:Nz
                Y(i,k) = sum(giks{i}(k,1:k).*Zr(1:k)*sqrt(2*pi*S*dt).*et(1:k),2);  % 2021-01-10
            end
        end
        
                    
        Nfailj = sum(sum(  abs(Y)>bik  ));

        if Nfailj==0  && alpa > 4 %  
            fprintf('\nNfailj=0!!! alpa = %6.4f  So set Pfj (p)=0',alpa);  Pfj = 0; 
        else   
            Pfj    = PFsum*(1/Nfailj); % Pfj = f(z)/H(z), is the weights of Z part, not include fxhx 2020
        end
        
       % fprintf('\nPfj = %12.4e\t ',Pfj);
    
% end
    


