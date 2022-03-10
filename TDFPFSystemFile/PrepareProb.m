function Prob = PrepareProb(Prob)
% parepare the parameter for Prob: Prob.S0, xuLocal, Prob.Nd,Prob.lb,Prob.dvDomain
% updeated: 2020-12-23
% updeated: 2021-05-29  2021-06-22

    % transform 
    if isfield(Prob,'Para')
        Prob = ParaStat(Prob,'PtoS');
    else
        Prob = ParaStat(Prob,'StoP');
    end

    % prepare the index of x_u
    Prob.nxr = length(Prob.dvLocal);
    ID = 1:Prob.Nx;
    if ~isfield(Prob,'xuLocal')  % 2021-08-10, 注：高维问题里面需定义Prob.xuLocal=[]，因其在加权中fxhx_no无需计算
        Prob.xuLocal = ID(ismember(ID,Prob.dvLocal)==0);        % the index of variables unrelated with sita
    end

    % prepare the 'Prob.dvDomain' and 'Prob.S0'
    Prob.dvDomain = {}; Prob.S0= [];   % change dsDomain to dvDomain
    p = 1;
    for i = 1:Prob.nxr   % number of x_r
        ixr = Prob.dvLocal(i);  % e.g. 1  OLD name: dvLocal
        isP = Prob.si_Para{i};  % e.g. [1, 2]
        Prob.S0 =  [Prob.S0, Prob.Para{ixr}(:,isP)];
        for j = 1:length(Prob.si_Para{i})  % number of si in each x_r
            Prob.dvDomain{p} = Prob.Domain_x_s{i}{j};   % e.g. {[11, 13],[0.8, 1.6]}
            % Prob.dvDomain{p} = cell2mat(Prob.Domain_x_s{i}(j));
            p = p+1;
        end
    end

    % prepare the 'Prob.lb' and 'Prob.ub'
    Prob.Nd = size(Prob.dvDomain,2); % number of sita
    for i = 1:Prob.Nd
        Prob.lb(i) = Prob.dvDomain{i}(1);
        Prob.ub(i) = Prob.dvDomain{i}(2);
    end

end

    
function Prob = ParaStat(Prob,Flag)
% 变量分布参数与统计参数的相互转换
% Transformation between distribution parameters and statistical parameters of variables
% Flag: 'StoP','PtoS'
% updated： 2020-11-30

    Nx = Prob.Nx;
    Dist = Prob.Dist;
    if strcmp(Flag,'StoP')
        if  isfield(Prob,'mu') && isfield(Prob,'sd')
            mu = Prob.mu ;  sd = Prob.sd;
            for i=1:Nx
                mux = mu(i); sdx = sd(i); 
                if strcmp(Dist{i}, 'unif')      Para{i} = [mux-sqrt(3)*sdx  mux+sqrt(3)*sdx]; 
                elseif strcmp(Dist{i}, 'exp')   Para{i} = [mux ];  
                elseif strcmp(Dist{i}, 'norm')  Para{i} = [mux  sdx ]; 
                elseif strcmp(Dist{i}, 'logn')  Para{i} = [log(mux)-1/2*log(1+(sdx/mux)^2)   sqrt(log((sdx/mux)^2+1)) ]; % 
                elseif strcmp(Dist{i}, 'rayl')  Para{i} = [mux/sqrt(pi/2) ]; 
                elseif strcmp(Dist{i}, 'gam')   Para{i} = [(mux/sdx)^2 sdx^2/mux]; 
                elseif strcmp(Dist{i}, 'ev')    Para{i} = [mux+0.577215*sqrt(6)*sdx/pi   sqrt(6)*sdx/pi];  
                else error('the distribution is not include in the ParaDist !');
                end
            end
            Prob.Para = Para;
        else
            fprintf('Error: There are no mu or sd!');
        end

    elseif strcmp(Flag,'PtoS')
        if isfield(Prob,'Para')
            Para = Prob.Para ;
            mu = zeros(1,Nx); sd = zeros(1,Nx);
            for i=1:Nx
                ab = Para{i};  % ab: [p1 p2]
                if strcmp(Dist{i}, 'unif'),  [mu(i), sd(i)] = unifstat(ab(1),ab(2)); 
                elseif strcmp(Dist{i}, 'exp'),   [mu(i), sd(i)] = expstat(ab(1)); 
                elseif strcmp(Dist{i}, 'norm'),  [mu(i), sd(i)] = normstat(ab(1),ab(2));
                elseif strcmp(Dist{i}, 'logn'),  [mu(i), sd(i)] = lognstat(ab(1),ab(2)); 
                elseif strcmp(Dist{i}, 'wbl'),   [mu(i), sd(i)] = wblstat(ab(1),ab(2));
                elseif strcmp(Dist{i}, 'rayl'),  [mu(i), sd(i)] = raylstat(ab(1));
                elseif strcmp(Dist{i}, 'gam'),   [mu(i), sd(i)] = gamstat(ab(1),ab(2));
                elseif strcmp(Dist{i}, 'ev'),    [mu(i), sd(i)] = evstat(ab(1),ab(2)); 
                else  error('the distribution is not include in the StatDist !');
                end 
            end
            sd=sqrt(sd);
            Prob.mu = mu ; 
            Prob.sd = sd ;
        else
            fprintf('Error: There are no Para !');
        end
    else
        fprintf('Error: No trans of Para&Stat has take place.')
    end

end

