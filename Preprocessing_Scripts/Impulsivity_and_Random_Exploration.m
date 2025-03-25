%% Test the relationship between random exploration Psychometric SLOPE AND impulsivity scores

% load impulsivity dataset

% impulsivity_table = readtable('Impulsivity_Scores_compatable.csv');
% 
% total_impulsivity_scores = double(impulsivity_table.Var33);

load('Impulsiveness_corrected.mat')

total_impulsivity_scores = impulsivity_corrected(:,end);


%% create sub_scale scores

attention_items = [5 6 9 11 20 24 26 28]+2;
motor_items = [2 3 4 16 17 19 21 22 23 25 30]+2;
nonplanning_items = [1 7 8 10 12 13 14 15 18 27 29]+2; 

attention_subscores = sum(impulsivity_corrected(:,attention_items),2);

motor_subscores = sum(impulsivity_corrected(:,motor_items),2);

nonplanning_subscores = sum(impulsivity_corrected(:,nonplanning_items),2);




%%

% load individual random epxloration metrics

full_random_exp_slopes = readtable('Replication_Random_Exploration_Slope_All.csv');


% data wrangling for short horizon 
Ind1_random_exp_short = double(full_random_exp_slopes.Ind1_rand_short_slope);
Ind2_random_exp_short = double(full_random_exp_slopes.Ind2_rand_short_slope);

short_random_exp = zeros(size(Ind1_random_exp_short,1) + size(Ind2_random_exp_short,1),1);

short_random_exp(1:2:length(short_random_exp)-1) = Ind1_random_exp_short;
short_random_exp(2:2:length(short_random_exp)) = Ind2_random_exp_short;

% correlation btw short random and impulsivity

% Compute Pearson correlation
[r, p, rlo, rup] = corrcoef(short_random_exp, total_impulsivity_scores);

% Extract relevant statistics
r_value = r(1,2);          % Correlation coefficient
p_value = p(1,2);          % P-value
ci_lower = rlo(1,2);       % Lower bound of 95% CI
ci_upper = rup(1,2);       % Upper bound of 95% CI

% Display results
fprintf('Pearson correlation: r = %.3f, 95%% CI [%.3f, %.3f], p = %.3g\n', ...
        r_value, ci_lower, ci_upper, p_value);
    
    
scatter(short_random_exp,total_impulsivity_scores)


[coef p_val] = corrcoef(short_random_exp, attention_subscores)

[coef p_val] =  corrcoef(short_random_exp, motor_subscores)

[coef p_val] = corrcoef(short_random_exp, nonplanning_subscores)



% data wrangling for long horizon 


Ind1_random_exp_long = double(full_random_exp_slopes.Ind1_rand_long_slope);
Ind2_random_exp_long = double(full_random_exp_slopes.Ind2_rand_long_slope);

long_random_exp = zeros(size(Ind1_random_exp_long,1) + size(Ind2_random_exp_long,1),1);

long_random_exp(1:2:length(long_random_exp)-1) = Ind1_random_exp_long;
long_random_exp(2:2:length(long_random_exp)) = Ind2_random_exp_long;



% correlation btw long random and impulsivity

% Compute Pearson correlation
[r, p, rlo, rup] = corrcoef(long_random_exp, total_impulsivity_scores);

% Extract relevant statistics
r_value = r(1,2);          % Correlation coefficient
p_value = p(1,2);          % P-value
ci_lower = rlo(1,2);       % Lower bound of 95% CI
ci_upper = rup(1,2);       % Upper bound of 95% CI

% Display results
fprintf('Pearson correlation: r = %.3f, 95%% CI [%.3f, %.3f], p = %.3g\n', ...
        r_value, ci_lower, ci_upper, p_value);
    
scatter(long_random_exp,total_impulsivity_scores)


[coef p_val] = corrcoef(long_random_exp, attention_subscores)

[coef p_val] =  corrcoef(long_random_exp, motor_subscores)

[coef p_val] = corrcoef(long_random_exp, nonplanning_subscores)


%% Test the relationship btw random exploration RATE and Impulsivity Scores

% load random exploration rate first row is dyad's score

load('Replication_Random_Exploration_Rate.mat')


random_exp(1:3:end,:) = []; 



% calculate correlations btw short and long horizons (2nd column is short and 3rd column is long horizon)


% correlation btw short horizon random exploration rate and impulsivity

% Compute Pearson correlation
[r, p, rlo, rup] = corrcoef(random_exp(:,2), total_impulsivity_scores);

% Extract relevant statistics
r_value = r(1,2);          % Correlation coefficient
p_value = p(1,2);          % P-value
ci_lower = rlo(1,2);       % Lower bound of 95% CI
ci_upper = rup(1,2);       % Upper bound of 95% CI

% Display results
fprintf('Pearson correlation: r = %.3f, 95%% CI [%.3f, %.3f], p = %.3g\n', ...
        r_value, ci_lower, ci_upper, p_value);
    
scatter(random_exp(:,2),total_impulsivity_scores)




[coef p_val] = corrcoef(random_exp(:,2), attention_subscores)

[coef p_val] =  corrcoef(random_exp(:,2), motor_subscores)

[coef p_val] = corrcoef(random_exp(:,2), nonplanning_subscores)









% correlation btw long horizon random exploration rate and impulsivity


% Compute Pearson correlation
[r, p, rlo, rup] = corrcoef(random_exp(:,3), total_impulsivity_scores);

% Extract relevant statistics
r_value = r(1,2);          % Correlation coefficient
p_value = p(1,2);          % P-value
ci_lower = rlo(1,2);       % Lower bound of 95% CI
ci_upper = rup(1,2);       % Upper bound of 95% CI

% Display results
fprintf('Pearson correlation: r = %.3f, 95%% CI [%.3f, %.3f], p = %.3g\n', ...
        r_value, ci_lower, ci_upper, p_value);
    
scatter(random_exp(:,3),total_impulsivity_scores,'filled')






[coef p_val] = corrcoef(random_exp(:,3), attention_subscores)

[coef p_val] =  corrcoef(random_exp(:,3), motor_subscores)

[coef p_val] = corrcoef(random_exp(:,3), nonplanning_subscores)






%% Test the relationship between directed exploration AND impulsivity scores

% load impulsivity dataset

load('Impulsiveness_corrected.mat')

total_impulsivity_scores = impulsivity_corrected(:,end);

% load individual directed epxloration metrics

full_directed_exp_PSE = readtable('Replication_Directed_Exploration_PSE_All.csv');


% data wrangling for short horizon 
Ind1_directed_exp_short = double(full_directed_exp_PSE.Ind1_directed_short_PSE);
Ind2_directed_exp_short = double(full_directed_exp_PSE.Ind2_directed_short_PSE);

short_directed_exp = zeros(size(Ind1_directed_exp_short,1) + size(Ind2_directed_exp_short,1),1);

short_directed_exp(1:2:length(short_directed_exp)-1) = Ind1_directed_exp_short;
short_directed_exp(2:2:length(short_directed_exp)) = Ind2_directed_exp_short;

% correlation btw short directed and impulsivity

% Compute Pearson correlation
[r, p, rlo, rup] = corrcoef(short_directed_exp, total_impulsivity_scores);

% Extract relevant statistics
r_value = r(1,2);          % Correlation coefficient
p_value = p(1,2);          % P-value
ci_lower = rlo(1,2);       % Lower bound of 95% CI
ci_upper = rup(1,2);       % Upper bound of 95% CI

% Display results
fprintf('Pearson correlation: r = %.3f, 95%% CI [%.3f, %.3f], p = %.3g\n', ...
        r_value, ci_lower, ci_upper, p_value);
    
scatter(total_impulsivity_scores,short_directed_exp)


% data wrangling for long horizon 


Ind1_directed_exp_long = double(full_directed_exp_PSE.Ind1_directed_long_PSE);
Ind2_directed_exp_long = double(full_directed_exp_PSE.Ind2_directed_long_PSE);

long_directed_exp = zeros(size(Ind1_directed_exp_long,1) + size(Ind2_directed_exp_long,1),1);

long_directed_exp(1:2:length(long_directed_exp)-1) = Ind1_directed_exp_long;
long_directed_exp(2:2:length(long_directed_exp)) = Ind2_directed_exp_long;



% correlation btw long directed and impulsivity

% Compute Pearson correlation
[r, p, rlo, rup] = corrcoef(long_directed_exp, total_impulsivity_scores);

% Extract relevant statistics
r_value = r(1,2);          % Correlation coefficient
p_value = p(1,2);          % P-value
ci_lower = rlo(1,2);       % Lower bound of 95% CI
ci_upper = rup(1,2);       % Upper bound of 95% CI

% Display results
fprintf('Pearson correlation: r = %.3f, 95%% CI [%.3f, %.3f], p = %.3g\n', ...
        r_value, ci_lower, ci_upper, p_value);
    
scatter(long_directed_exp,total_impulsivity_scores)








%% correlations with accuracy - SHORT

% load impulsivity dataset

load('Impulsiveness_corrected.mat')

total_impulsivity_scores = impulsivity_corrected(:,end);

% load individual random epxloration metrics

short_acc = readtable('Replication_Short_Acc.csv');


% data wrangling for short horizon 
Ind1_short_acc = double(short_acc.Ind1_Short_Acc);
Ind2_short_acc = double(short_acc.Ind2_Short_Acc);

short_acc = zeros(size(Ind1_short_acc,1) + size(Ind2_short_acc,1),1);

short_acc(1:2:length(short_acc)-1) = Ind1_short_acc;
short_acc(2:2:length(short_acc)) = Ind2_short_acc;


% correlation btw short directed and impulsivity

[coef,p] = corr(short_acc,total_impulsivity_scores)

scatter(short_acc,total_impulsivity_scores)

%% correlations with accuracy - LONG
% load impulsivity dataset

load('Impulsiveness_corrected.mat')

total_impulsivity_scores = impulsivity_corrected(:,end);

% load individual random epxloration metrics

long_acc = readtable('Replication_Long_Acc.csv');
long_acc = table2array(long_acc)

% exclude dyad scores
long_acc(:,3:3:end) = [];

ind1_all_scores = long_acc(:,1:2:end-1);
ind2_all_scores = long_acc(:,2:2:end);





% data wrangling for short horizon 
Ind1_long_mean_acc = mean(ind1_all_scores,2);
Ind2_long_mean_acc =mean(ind2_all_scores,2);

long_acc = zeros(size(Ind1_long_mean_acc,1) + size(Ind2_long_mean_acc,1),1);

long_acc(1:2:length(long_acc)-1) = Ind1_long_mean_acc;
long_acc(2:2:length(long_acc)) = Ind2_long_mean_acc;


% correlation btw short directed and impulsivity

[coef,p] = corr(long_acc,total_impulsivity_scores)

scatter(long_acc,total_impulsivity_scores)


% correlation per each choice sequence
% load impulsivity dataset
load('Impulsiveness_corrected.mat')

total_impulsivity_scores = impulsivity_corrected(:,end);

% load individual random epxloration metrics

long_acc = readtable('Replication_Long_Acc.csv');
long_acc = table2array(long_acc)

% exclude dyad scores
long_acc(:,3:3:end) = [];

% have all choice sequence
ind1_all_scores = long_acc(:,1:2:end-1);
ind2_all_scores = long_acc(:,2:2:end);

long_acc = zeros(34,6);

long_acc(1:2:length(long_acc)-1,:) = ind1_all_scores;
long_acc(2:2:length(long_acc),:) = ind2_all_scores;



for i = 1:size(long_acc,2);
    
    % correlation btw short directed and impulsivity

[coefs(i), p_values(i)] = corr(long_acc(:,i),total_impulsivity_scores);


end
