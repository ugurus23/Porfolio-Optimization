% Load .csv data from yahoo finance
CME1 = readmatrix('CME.csv');    % Jan 14 to Jan 22 monthly data
CME2 = readmatrix('CME-2.csv');   % Feb 14 to Feb 22 monthly data
BR1 = readmatrix('BR.csv');      % Jan 14 to Jan 22 monthly data
BR2 = readmatrix('BR-3.csv');     % Feb 14 to Feb 22 monthly data
CBOE1 = readmatrix('CBOE.csv');  % Jan 14 to Jan 22 monthly data
CBOE2 = readmatrix('CBOE-2.csv'); % Feb 14 to Feb 22 monthly data
ICE1 = readmatrix('ICE.csv');    % Jan 14 to Jan 22 monthly data
ICE2 = readmatrix('ICE-2.csv');   % Feb 14 to Feb 22 monthly data
ACN1 = readmatrix('ACN.csv');    % Jan 14 to Jan 22 monthly data
ACN2 = readmatrix('ACN-2.csv');   % Feb 14 to Feb 22 monthly data

% Return (r_it) = (Month 2 Close - Month 1 Close)/Month 1 Close
CME_rtn = (CME2(:, 6) - CME1(:,6))./CME1(:, 6);      % mo. return
BR_rtn  = (BR2(:, 6) - BR1(:,6))./BR1(:, 6);         % mo. return
CBO_rtn = (CBOE2(:, 6) - CBOE1(:,6))./CBOE1(:, 6);   % mo. return
ICE_rtn = (ICE2(:, 6) - ICE1(:,6))./ICE1(:, 6);      % mo. return
ACN_rtn = (ACN2(:, 6) - ACN1(:,6))./ACN1(:, 6);      % mo. return

% Calculate arithmetic means (rbar_i)
CME_rbar = mean(CME_rtn)
BR_rbar  = mean(BR_rtn)
CBO_rbar = mean(CBO_rtn)
ICE_rbar = mean(ICE_rtn)
ACN_rbar = mean(ACN_rtn)

% Calculate geometric Means (mu_i)
CME_rtn_t = 1;      % initialize
BR_rtn_t = 1;      % intiialize
CBO_rtn_t = 1;      % initialize
ICE_rtn_t = 1;      % intiialize
ACN_rtn_t = 1;      % initialize

% loop through to calculate total product of returns
for t = 1:97        %T = 97 data points
    CME_rtn_t = (1 + CME_rtn(t,:))*CME_rtn_t;   % cumulative
    BR_rtn_t = (1 + BR_rtn(t,:))*BR_rtn_t;      % cumulative
    CBO_rtn_t = (1 + CBO_rtn(t,:))*CBO_rtn_t;   % cumulative
    ICE_rtn_t = (1 + ICE_rtn(t,:))*ICE_rtn_t;   % cumulative
    ACN_rtn_t = (1 + ACN_rtn(t,:))*ACN_rtn_t;   % cumulative
end

% take the 'T'th root to get final geometric mean (mu_i)
CME_mu = CME_rtn_t^(1/97) - 1
BR_mu = BR_rtn_t^(1/97) - 1
CBO_mu = CBO_rtn_t^(1/97) - 1
ICE_mu = ICE_rtn_t^(1/97) - 1
ACN_mu = ACN_rtn_t^(1/97) - 1

% Calculate standard deviations
CME_std = std(CME_rtn)
BR_std = std(BR_rtn)
CBO_std = std(CBO_rtn)
ICE_std = std(ICE_rtn)
ACN_std = std(ACN_rtn)

SPY1 = readmatrix('SPY-6.csv')    
SPY2 = readmatrix('SPY-7.csv')
SPY_rtn = (SPY2(:, 6) - SPY1(:,6))./SPY1(:, 6); 
SPY_rbar = mean(SPY_rtn)

EEMV1 = readmatrix('EEMV-3.csv')    
EEMV2 = readmatrix('EEMV-4.csv')
EEMV_rtn = (EEMV2(:, 6) - EEMV1(:,6))./EEMV1(:, 6); 
EEMV_rbar = mean(EEMV_rtn)

GOVT1 = readmatrix('GOVT-4.csv')    
GOVT2 = readmatrix('GOVT-5.csv')
GOVT_rtn = (GOVT2(:, 6) - GOVT1(:,6))./GOVT1(:, 6); 
GOVT_rbar = mean(GOVT_rtn)

% Calculate geometric Means (mu_i)
SPY_rtn_t = 1;      % initialize
GOVT_rtn_t = 1;      % intiialize
EEMV_rtn_t = 1;      % initialize

% loop through to calculate total product of returns
for t = 1:97       %T = 98 data points
    SPY_rtn_t = (1 + SPY_rtn(t,:))*SPY_rtn_t;   % cumulative
    GOVT_rtn_t = (1 + GOVT_rtn(t,:))*GOVT_rtn_t;   % cumulative
    EEMV_rtn_t = (1 + EEMV_rtn(t,:))*EEMV_rtn_t;   % cumulative
end

% take the 'T'th root to get final geometric mean (mu_i)
SPY_mu = SPY_rtn_t^(1/97) - 1
GOVT_mu = GOVT_rtn_t^(1/97) - 1
EEMV_mu = EEMV_rtn_t^(1/97) - 1

SPY_std = std(SPY_rtn)
GOVT_std = std(GOVT_rtn)
EEMV_std = std(EEMV_rtn)
% Calculate covariances (SPY: i=1, GOV: i=2, EEM: i=3)
rtn_matrix = [SPY_rtn GOVT_rtn EEMV_rtn];     % i
format shortE
cov_matrix = cov(rtn_matrix,1)

R = linspace(0.03, 0.18, 15);    % expected return data points
H = cov_matrix;                  % covariance matrix
f = zeros(3,1);                  % [0 0 0]T

Aeq = [-SPY_mu -GOVT_mu -EEMV_mu; ones(1,3)]  %Equal. constr. matrix

% initiate vector, matrix to store values from loop:
p_var = zeros(size(R,1),1);     % portfolio variance (obj function)
w = zeros(3, size(R,1));        % asset weights (dec. variable)
lb = [-inf; -inf; -inf];
ub = [inf; inf; inf];

% lb = [0; 0; 0]

% loop through i from 1 to 15
for i = 1:size(R,2)
    beq = [-R(i); 1];                % Equality constraint vector
    % quadratic program
    
    [w(:,i), p_var(i)] = quadprog(H, f, [], [], Aeq, beq,lb,ub);
end

% Results table of: 'R', 'p_var', w_1, w_2, w_3:
T = array2table([R',p_var', w']);
T.Properties.VariableNames(1:5) = {'exp. rtn', 'port. var', 'w_1', 'w_2', 'w_3'}

display(T)

% Calculate covariances 
% (SPY=1, GOV=2, EEM=3,CME=4, BR=5, CBO=6, ICE=7, ACN=8)

rtn_matrix2 = [SPY_rtn GOVT_rtn EEMV_rtn CME_rtn BR_rtn CBO_rtn ICE_rtn ACN_rtn];
format shortE
cov_matrix2 = cov(rtn_matrix2,1)



% Generating efficient frontier using all 8 assets with short selling:

R2 = linspace(0.03, 0.18, 15);     % expected return data points
H2 = cov_matrix2;                  % covariance matrix
f2 = zeros(8,1);                    % [0 0 0 0 0 0 0 0]T
 
Aeq2 = [SPY_mu GOVT_mu EEMV_mu CME_mu BR_mu CBO_mu ICE_mu ACN_mu; ones(1,8)];  %Equality constraint matrix

% initiate vector, matrix to store values from loop:
p_var2 = zeros(size(R2,1),1);  % portfolio variance (obj function)
w2 = zeros(8, size(R2,1));     % asset weights (dec. variable)

% loop through i from 1 to 15
for i = 1:size(R2,2)
    beq2 = [R2(i); 1];          % Equality constraint vector
    % quadratic program
    [w2(:,i), p_var2(i)] = quadprog(H2, f2, [], [], Aeq2, beq2,[],[]);
end

% Results table of: 'R', 'p_var', w_i
T2 = array2table([R2',p_var2', w2']);
T2.Properties.VariableNames(1:10) = {'exp return', 'port. var','w_1', 'w_2', 'w_3', 'w_4', 'w_5', 'w_6', 'w_7', 'w_8'}
display(T2)

% Plot efficient frontier
figure('Name','Efficient Frontier - With Shorting');
plot(p_var2, R2, 'b-*');
hold on
plot(p_var, R, 'r--o');
legend('8-asset portfolio', '3-asset portfolio')
title('Figure 3 - Efficient Frontier with 8 assets');
xlabel('portfolio variance');
ylabel('expected return');