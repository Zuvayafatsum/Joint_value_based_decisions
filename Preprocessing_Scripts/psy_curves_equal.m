% Equal condition


% Psychometric curve for equalvs. mean differences

clc
clear all

directory = 'D:\Program Files\MATLAB\Joint_Horizon\main_study\Collective\';
data_files = dir('D:\Program Files\MATLAB\Joint_Horizon\main_study\Collective\');

data_file_names = {data_files.name}';
data_file_names(1:2,:) = [];


 

gavg_tenhor = []; 
gavg_fivehor = []; 
data_type = {}; 
meta_data_s_e_c = [];
higher_info_choice_only_equal=[];
data_clean = [];

%% SHORT EQUAL Collective
for t = 1: length(data_file_names);

    data_ID_c = string(data_file_names(t));

    
   load('D:\Program Files\MATLAB\Joint_Horizon\main_study\Collective\' + data_ID_c);



% Determines higher information for short trials ONLY

%%% Filter for horizon 5 & 10 only

answer = [game.mean]'
trial_type = [game.gameLength]'


short_trials = find(trial_type == 5)
long_trials = find(trial_type == 10)
short_only_data = game(short_trials)
long_only_data = game(long_trials)

forced_short = {short_only_data.nforced}.'
forced_long = {long_only_data.nforced}.'
condition_unequal_majority = zeros(size(forced_short, 1), 1);
forced_short =cell2mat(forced_short)
    
    for i = 1:size(forced_short, 1)
        countOfOnes = sum(forced_short(i, :) == 1)
        if countOfOnes == 2;
            condition_unequal_majority(i) = 0;
        elseif countOfOnes == 1;
            condition_unequal_majority(i) = 2;
        else
            condition_unequal_majority(i) = 1;
        end
    end
    
  equal_mask = find(condition_unequal_majority == 0);
  equal_trials_only_short = short_only_data(equal_mask) %checked
    
  short_trials_answers = {equal_trials_only_short.a}.'
  short_trials_answers = cell2mat(short_trials_answers)

    % Accuracy is choosing the option on the right

  short_responses = short_trials_answers(:,5)

  accuracy_short = zeros(length(short_responses),1);
    
  for i = 1 : length(short_responses);
    
   if short_responses(i) == 2 % 2 is the option on the right
        accuracy_short(i) = 1;
   else 
        accuracy_short(i) = 0;
   end
 end
    
    accuracy_fivehorizon = sum(accuracy_short)/length(accuracy_short)



means = [equal_trials_only_short.mean]'
mean_difference = [means(:,2)-means(:,1)]



data_clean(:,1) = 2*ones(40,1) %right
data_clean(:,2) = accuracy_short
data_clean(:,3) = mean_difference




data_ID_c = char(data_ID_c)

data_clean(:,6)=str2double(data_ID_c(5:6))


   if data_ID_c(8) == 'C' | data_ID_c(9) == 'C' 
       data_clean(:,4) = 1
   else
       data_clean(:,4) = 0;
   end

   if data_clean(:,4) == 0;
   data_clean(:,5) = str2double(data_ID_c(10:11))
   else 
   data_clean(:,5) = 0
   end




meta_data_s_e_c = [meta_data_s_e_c ; data_clean];

data_clean = [];

higher_info_choice_only_equal=[];

end




data = [meta_data_s_e_c(:,3), meta_data_s_e_c(:,2)]

% Unique stimulus intensity levels
s_e_c_intensity_levels = unique(data(:, 1));

% Calculate the probability of correct responses for each intensity level
s_e_c_prob_correct = arrayfun(@(x) mean(data(data(:, 1) == x, 2)), s_e_c_intensity_levels);

s_e_c_std_error = arrayfun(@(x) std(data(data(:, 1) == x, 2))/sqrt(sum(data(:, 1) == x)), s_e_c_intensity_levels);



%% SHORT EQUAL INDIVIDUAL 


directory = 'D:\Program Files\MATLAB\Joint_Horizon\main_study\Individual\';
data_files = dir('D:\Program Files\MATLAB\Joint_Horizon\main_study\Individual\');

data_file_names = {data_files.name}';
data_file_names(1:2,:) = [];


%% 

gavg_tenhor = []; 
gavg_fivehor = []; 
data_type = {}; 
meta_data_s_e_i = [];
data = [];
data_clean = [] ;


for t = 1: length(data_file_names);

    data_ID_i = string(data_file_names(t));

    
   load('D:\Program Files\MATLAB\Joint_Horizon\main_study\Individual\' + data_ID_i);


% Determines higher information for short trials ONLY

%%% Filter for horizon 5 & 10 only

answer = [game.mean]'
trial_type = [game.gameLength]'


short_trials = find(trial_type == 5)
long_trials = find(trial_type == 10)
short_only_data = game(short_trials)
long_only_data = game(long_trials)

forced_short = {short_only_data.nforced}.'
forced_long = {long_only_data.nforced}.'
condition_unequal_majority = zeros(size(forced_short, 1), 1);
forced_short =cell2mat(forced_short)
    
    for i = 1:size(forced_short, 1)
        countOfOnes = sum(forced_short(i, :) == 1)
        if countOfOnes == 2;
            condition_unequal_majority(i) = 0;
        elseif countOfOnes == 1;
            condition_unequal_majority(i) = 2;
        else
            condition_unequal_majority(i) = 1;
        end
    end
    
  equal_mask = find(condition_unequal_majority == 0);
  equal_trials_only_short = short_only_data(equal_mask) %checked
    
  short_trials_answers = {equal_trials_only_short.a}.'
  short_trials_answers = cell2mat(short_trials_answers)

    % Accuracy is choosing the option on the right

  short_responses = short_trials_answers(:,5)

  accuracy_short = zeros(length(short_responses),1);
    
  for i = 1 : length(short_responses);
    
   if short_responses(i) == 2 % 2 is the option on the right
        accuracy_short(i) = 1;
   else 
        accuracy_short(i) = 0;
   end
 end

accuracy_fivehorizon = sum(accuracy_short)/length(accuracy_short)

% % Sanity check - completed
% sanitycheck(:,1) = short_responses
% sanitycheck(:,2) = higher_info_choice_only_unequal
% sanitycheck(:,3) = accuracy_short


   
means = [equal_trials_only_short.mean]'
mean_difference = [means(:,2)-means(:,1)]

%Short horizon
data_clean(:,1)= 2*ones(40,1)
data_clean(:,2)= accuracy_short %Correct or incorrect
data_clean(:,3) = mean_difference


data_ID_i = char(data_ID_i)

data_clean(:,6)=str2double(data_ID_i(5:6))


   if data_ID_i(8) == 'C' | data_ID_i(9) == 'C' 
       data_clean(:,4) = 1
   else
       data_clean(:,4) = 0;
   end

   if data_clean(:,4) == 0;
   data_clean(:,5) = str2double(data_ID_i(10:11))
   else 
   data_clean(:,5) = 0
   end



meta_data_s_e_i = [meta_data_s_e_i; data_clean];

data_clean = [];
end




data = [meta_data_s_e_i(:,3), meta_data_s_e_i(:,2)]

% Unique stimulus intensity levels
s_e_i_intensity_levels = unique(data(:, 1));

% Calculate the probability of correct responses for each intensity level
s_e_i_prob_correct = arrayfun(@(x) mean(data(data(:, 1) == x, 2)), s_e_i_intensity_levels);
s_e_i_std_error = arrayfun(@(x) std(data(data(:, 1) == x, 2))/sqrt(sum(data(:, 1) == x)), s_e_i_intensity_levels);




%% Long EQUAL Collective
clc



gavg_tenhor = []; 
gavg_fivehor = []; 
data_type = {}; 
meta_data_l_e_c = [];
data = [];
data_clean = [] ;

directory = 'D:\Program Files\MATLAB\Joint_Horizon\main_study\Collective\';
data_files = dir('D:\Program Files\MATLAB\Joint_Horizon\main_study\Collective\');

data_file_names = {data_files.name}';
data_file_names(1:2,:) = [];


for t = 1: length(data_file_names);

    data_ID_c = string(data_file_names(t));

    
   load('D:\Program Files\MATLAB\Joint_Horizon\main_study\Collective\' + data_ID_c);


   

% Determines higher information for short trials ONLY

%%% Filter for horizon 5 & 10 only

answer = [game.mean]'
trial_type = [game.gameLength]'


short_trials = find(trial_type == 5)
long_trials = find(trial_type == 10)
short_only_data = game(short_trials)
long_only_data = game(long_trials)

forced_short = {short_only_data.nforced}.'
forced_long = {long_only_data.nforced}.'

%Horizon 10 - "accuracy" is probability of choosing more informative option

condition_unequal_majority = zeros(size(forced_long, 1), 1);
more_informative_long= zeros(size(forced_long, 1), 1);
forced_long =cell2mat(forced_long)

for i = 1:size(forced_long, 1)
    countOfOnes = sum(forced_long(i, :) == 1)
    if countOfOnes == 2;
        condition_unequal_majority(i) = 0;
        more_informative_long(i) = 0;
    elseif countOfOnes == 1;
        condition_unequal_majority(i) = 2;
        more_informative_long(i)=1
    else
        condition_unequal_majority(i) = 1;
        more_informative_long(i)=2;
    end
end

more_informative_long = more_informative_long'

equal_mask = find(condition_unequal_majority == 0);
equal_trials_only_long = long_only_data(equal_mask)






%Take last element of each array in answers and store them in a separate
%variable, then compare with ground truth

long_trials_answers = {equal_trials_only_long.a}.'
%long_trials_answers = {unequal_trials_only_long.a}.'
long_trials_answers = cell2mat(long_trials_answers)
%long_trials_answers = cell2mat(long_trials_answers)



%For horizon 5 only

long_responses = long_trials_answers(:,5)

accuracy_long = zeros(length(long_responses),1);

for i = 1 : length(long_responses);
    
   if long_responses(i) == 2 % 2 is the option on the right
        accuracy_long(i) = 1;
   else 
        accuracy_long(i) = 0;
   end
 end
accuracy_tenhorizon = sum(accuracy_long)/length(accuracy_long)

% % Sanity check - completed
% sanitycheck(:,1) = long_responses
% sanitycheck(:,2) = higher_info_choice_only_unequal
% sanitycheck(:,3) = accuracy_long






means = [equal_trials_only_long.mean]'
mean_difference = [means(:,2)-means(:,1)]

%long horizon
data_clean(:,1)= 2*ones(40,1)
data_clean(:,2)= accuracy_long %Correct or incorrect
data_clean(:,3) = mean_difference





data_ID_c = char(data_ID_c)

data_clean(:,6)=str2double(data_ID_c(5:6))


   if data_ID_c(8) == 'C' | data_ID_c(9) == 'C' 
       data_clean(:,4) = 1
   else
       data_clean(:,4) = 0;
   end

   if data_clean(:,4) == 0;
   data_clean(:,5) = str2double(data_ID_c(10:11))
   else 
   data_clean(:,5) = 0
   end





meta_data_l_e_c = [meta_data_l_e_c; data_clean];

data_clean = [];
end




data = [meta_data_l_e_c(:,3), meta_data_l_e_c(:,2)]

% Unique stimulus intensity levels
l_e_c_intensity_levels = unique(data(:, 1));

% Calculate the probability of correct responses for each intensity level
l_e_c_prob_correct = arrayfun(@(x) mean(data(data(:, 1) == x, 2)), l_e_c_intensity_levels);
l_e_c_std_error = arrayfun(@(x) std(data(data(:, 1) == x, 2))/sqrt(sum(data(:, 1) == x)), l_e_c_intensity_levels);




%% LONG equal INDIVIDUAL 


directory = 'D:\Program Files\MATLAB\Joint_Horizon\main_study\Individual\';
data_files = dir('D:\Program Files\MATLAB\Joint_Horizon\main_study\Individual\');

data_file_names = {data_files.name}';
data_file_names(1:2,:) = [];


gavg_tenhor = []; 
gavg_fivehor = []; 
data_type = {}; 
meta_data_l_e_i= [];
data = [];
data_clean = [] ;

%% 

gavg_tenhor = []; 
gavg_fivehor = []; 
data_type = {}; 
meta_data_l_e_i = [];

for t = 1: length(data_file_names);

    data_ID_i = string(data_file_names(t));

    
   load('D:\Program Files\MATLAB\Joint_Horizon\main_study\Individual\' + data_ID_i);


% Determines higher information for long trials ONLY

%%% Filter for horizon 5 & 10 only

answer = [game.mean]'
trial_type = [game.gameLength]'


short_trials = find(trial_type == 5)
long_trials = find(trial_type == 10)
short_only_data = game(short_trials)
long_only_data = game(long_trials)

forced_short = {short_only_data.nforced}.'
forced_long = {long_only_data.nforced}.'

%Horizon 5 - "accuracy" is probability of choosing more informative option

condition_unequal_majority = zeros(size(forced_long, 1), 1);
more_informative_long = zeros(size(forced_long, 1), 1);
forced_long =cell2mat(forced_long)

for i = 1:size(forced_long, 1)
    countOfOnes = sum(forced_long(i, :) == 1)
    if countOfOnes == 2;
        condition_unequal_majority(i) = 0;
        more_informative_long(i) = 0;
    elseif countOfOnes == 1;
        condition_unequal_majority(i) = 2;
        more_informative_long(i)=1
    else
        condition_unequal_majority(i) = 1;
        more_informative_long(i)=2;
    end
end

more_informative_long = more_informative_long'

equal_mask = find(condition_unequal_majority == 0);
equal_trials_only_long = long_only_data(equal_mask)
% higher_info_choice_only_equal= more_informative_long(more_informative_long ~= 0)
% higher_info_choice_only_equal= higher_info_choice_only_unequal'





%Take last element of each array in answers and store them in a separate
%variable, then compare with ground truth

long_trials_answers = {equal_trials_only_long.a}.'
%long_trials_answers = {equal_trials_only_long.a}.'
long_trials_answers = cell2mat(long_trials_answers)
%long_trials_answers = cell2mat(long_trials_answers)



%For horizon 5 only

long_responses = long_trials_answers(:,5)

accuracy_long = zeros(length(long_responses),1);

% higher_info_choice_only_equal= higher_info_choice_only_unequal'
for i = 1 : length(long_responses);
    
   if long_responses(i) == 2 % 2 is the option on the right
        accuracy_long(i) = 1;
   else 
        accuracy_long(i) = 0;
   end
 end

accuracy_tenhorizon = sum(accuracy_long)/length(accuracy_long)

% % Sanity check - completed
% sanitycheck(:,1) = long_responses
% sanitycheck(:,2) = higher_info_choice_only_unequal
% sanitycheck(:,3) = accuracy_long






means = [equal_trials_only_long.mean]'
mean_difference = [means(:,2)-means(:,1)]

%long horizon
data_clean(:,1)= 2*ones(40,1)
data_clean(:,2)= accuracy_long %Correct or incorrect
data_clean(:,3) = mean_difference




data_ID_i = char(data_ID_i)

data_clean(:,6)=str2double(data_ID_i(5:6))


   if data_ID_i(8) == 'C' | data_ID_i(9) == 'C' 
       data_clean(:,4) = 1
   else
       data_clean(:,4) = 0;
   end

   if data_clean(:,4) == 0;
   data_clean(:,5) = str2double(data_ID_i(10:11))
   else 
   data_clean(:,5) = 0
   end




meta_data_l_e_i = [meta_data_l_e_i; data_clean];

data_clean = [];
end




data = [meta_data_l_e_i(:,3), meta_data_l_e_i(:,2)]

% Unique stimulus intensity levels
l_e_i_intensity_levels = unique(data(:, 1));

% Calculate the probability of correct responses for each intensity level
l_e_i_prob_correct = arrayfun(@(x) mean(data(data(:, 1) == x, 2)), l_e_i_intensity_levels);
l_e_i_std_error = arrayfun(@(x) std(data(data(:, 1) == x, 2))/sqrt(sum(data(:, 1) == x)), l_e_i_intensity_levels);



%% PLOPTTING

green = [39 83 22 ]/ 256;
purple = [112, 48, 160]/256;
capsize = 24;
linewidth = 8;
dotsize = 256;
%% plot for col vs ind only for equal and short horizon

errorbar(s_e_i_intensity_levels, s_e_i_prob_correct, s_e_i_std_error, 'o', 'MarkerSize', 12, 'MarkerEdgeColor', green, 'MarkerFaceColor', green,'LineWidth',2,'Color',green)
hold on;

% Bahador's code
bhat = glmfit(s_e_i_intensity_levels, [s_e_i_prob_correct ones(size(s_e_i_prob_correct))], 'binomial', 'link', 'probit');

estimatedMean = -bhat(1)/bhat(2);
estimatedSD = 1/bhat(2);

% plotting 
hold on
C = 1.3 .* (min(s_e_i_intensity_levels) : 0.01 : max(s_e_i_intensity_levels));
Ind_Short = plot(C, cdf('norm', C, estimatedMean, estimatedSD),'LineWidth', 6, 'color', green, 'DisplayName', 'Individual Short');

hold on

errorbar(s_e_c_intensity_levels, s_e_c_prob_correct, s_e_c_std_error, 'o', 'MarkerSize', 12, 'MarkerEdgeColor', purple, 'MarkerFaceColor', purple,'LineWidth',2,'Color',purple)
hold on;

% Bahador's code
bhat = glmfit(s_e_c_intensity_levels, [s_e_c_prob_correct ones(size(s_e_c_prob_correct))], 'binomial', 'link', 'probit');

estimatedMean = -bhat(1)/bhat(2);
estimatedSD = 1/bhat(2);

% plotting 
hold on
D = 1.3 .* (min(s_e_c_intensity_levels) : 0.01 : max(s_e_c_intensity_levels));
Coll_Short = plot(D, cdf('norm', D, estimatedMean, estimatedSD), 'LineWidth', 6, 'color',purple, 'DisplayName', 'Collective Short');



xlabel('')
ylabel('')
title('')
set(gca, 'FontSize', 28, 'FontWeight', 'bold'); % Change the font size of axis labels
set(gca, 'LineWidth', 2); % Change the width of the axis lines
set(gcf, 'Color', 'w')
set(gca, 'TickDir', 'out')
% legend([Ind_Short, Coll_Short], 'Location', 'southeast')
% hLegend = findobj(gcf, 'Type', 'Legend');
% set(hLegend, 'FontSize', 32, 'FontWeight', 'bold');% Remove top and right ticks
ax = gca;
ax.XAxis.TickLength = [0.015 0];  % Bottom ticks only
ax.YAxis.TickLength = [0.015 0];  % Left ticks only
ax.Box = 'off';  % Turn off the box
ax.TickLength = [0.015 0];  % Set the tick length
% legend boxoff
axis square

set(gca, 'FontSize', 36, "FontWeight","bold"); % Change the font size of axis labels
set(gca, 'LineWidth', 8); % Change the width of the axis lines
xlabel('')
ylabel('')
title('')
set(gcf, 'Color', 'w')
set(gca, 'TickDir', 'out')



%% figure 2, equal information , ind vs collective for long horizon
green = [39 83 22 ]/ 256;
purple = [112, 48, 160]/256;
capsize = 24;
linewidth = 8;
dotsize = 256;

% same for long col vs ind

errorbar(l_e_i_intensity_levels, l_e_i_prob_correct, l_e_i_std_error, 'o', 'MarkerSize', 12, 'MarkerEdgeColor', green, 'MarkerFaceColor', green,'LineWidth',2,'Color',green)
hold on;

% Bahador's code
bhat = glmfit(l_e_i_intensity_levels, [l_e_i_prob_correct ones(size(l_e_i_prob_correct))], 'binomial', 'link', 'probit');

estimatedMean = -bhat(1)/bhat(2);
estimatedSD = 1/bhat(2);

% plotting 
hold on
E = 1.3 .* (min(l_e_i_intensity_levels) : 0.01 : max(l_e_i_intensity_levels));
Ind_Long = plot(E, cdf('norm', E, estimatedMean, estimatedSD), 'Color',green, 'LineWidth', 6, 'color', green, 'DisplayName', 'Individual Long');

hold on

errorbar(l_e_c_intensity_levels, l_e_c_prob_correct, l_e_c_std_error, 'o', 'MarkerSize', 12, 'MarkerEdgeColor', purple, 'MarkerFaceColor', purple,'LineWidth',2,'Color',purple)
hold on;

% Bahador's code
bhat = glmfit(l_e_c_intensity_levels, [l_e_c_prob_correct ones(size(l_e_c_prob_correct))], 'binomial', 'link', 'probit');

estimatedMean = -bhat(1)/bhat(2);
estimatedSD = 1/bhat(2);

% plotting 
hold on
F = 1.3 .* (min(l_e_c_intensity_levels) : 0.01 : max(l_e_c_intensity_levels));
Coll_Long = plot(F, cdf('norm', F, estimatedMean, estimatedSD), 'LineWidth', 6, 'color', purple, 'DisplayName', 'Collective Long');

hold on 
xlabel('')
ylabel('')
title('')
set(gca, 'FontSize', 28, 'FontWeight', 'bold'); % Change the font size of axis labels
set(gca, 'LineWidth', 2); % Change the width of the axis lines
set(gcf, 'Color', 'w')
set(gca, 'TickDir', 'out')
% legend([Ind_Long, Coll_Long], 'Location', 'southeast')
% hLegend = findobj(gcf, 'Type', 'Legend');
% set(hLegend, 'FontSize', 32, 'FontWeight', 'bold');% Remove top and right ticks
ax = gca;
ax.XAxis.TickLength = [0.015 0];  % Bottom ticks only
ax.YAxis.TickLength = [0.015 0];  % Left ticks only
ax.Box = 'off';  % Turn off the box
ax.TickLength = [0.015 0];  % Set the tick length
% legend boxoff

axis square

set(gca, 'FontSize', 36, "FontWeight","bold"); % Change the font size of axis labels
set(gca, 'LineWidth', 8); % Change the width of the axis lines
xlabel('')
ylabel('')
title('')
set(gcf, 'Color', 'w')
set(gca, 'TickDir', 'out')

%% Across horizons
% Combine short and long conditions



short_combined = [s_e_i_prob_correct s_e_c_prob_correct];
long_combined = [l_e_i_prob_correct l_e_c_prob_correct]


mean_short_horizon = mean(short_combined,2)

mean_long_horizon = mean(long_combined,2)


se_short_horizon_combined = [s_e_i_std_error s_e_c_std_error]

se_long_horizon_combined = [l_e_i_std_error l_e_c_std_error]



se_short_horizon = mean(se_short_horizon_combined,2)

se_long_horizon = mean(se_long_horizon_combined,2)



%% plot for long vs short across conditions

errorbar(s_e_i_intensity_levels, mean_short_horizon, se_short_horizon, 'o', 'MarkerSize', 12, 'MarkerEdgeColor', 'black', 'MarkerFaceColor', 'black','LineWidth',2,'Color','black')
hold on;

% Bahador's code
bhat = glmfit(s_e_i_intensity_levels, [mean_short_horizon ones(size(mean_short_horizon))], 'binomial', 'link', 'probit');

estimatedMean = -bhat(1)/bhat(2);
estimatedSD = 1/bhat(2);

% plotting 
hold on
C = 1.3 .* (min(s_e_i_intensity_levels) : 0.01 : max(s_e_i_intensity_levels));
Short = plot(C, cdf('norm', C, estimatedMean, estimatedSD), 'b-', 'LineWidth', 6, 'color', 'black', 'DisplayName', 'Short Horizon');

hold on

errorbar(s_e_c_intensity_levels, mean_long_horizon, se_long_horizon, 'o', 'MarkerSize', 12, 'MarkerEdgeColor', [0.8 0.39 0.039], 'MarkerFaceColor', [0.8 0.39 0.039],'LineWidth',2,'Color',[0.8 0.39 0.039])
hold on;

% Bahador's code
bhat = glmfit(s_e_c_intensity_levels, [mean_long_horizon ones(size(mean_long_horizon))], 'binomial', 'link', 'probit');

estimatedMean = -bhat(1)/bhat(2);
estimatedSD = 1/bhat(2);

% plotting 
hold on
D = 1.3 .* (min(s_e_c_intensity_levels) : 0.01 : max(s_e_c_intensity_levels));
Long = plot(D, cdf('norm', D, estimatedMean, estimatedSD), 'LineWidth', 6, 'color', [0.8 0.39 0.039], 'DisplayName', 'Long Horizon');

xlabel('Mean Difference (Right - Left)')
ylabel('P(Choosing Right)')
title('Random Exploration - Short VS Long')
set(gca, 'FontSize', 28, 'FontWeight', 'bold'); % Change the font size of axis labels
set(gca, 'LineWidth', 2); % Change the width of the axis lines
set(gcf, 'Color', 'w')
set(gca, 'TickDir', 'out')
legend([Short, Long], 'Location', 'southeast')
hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend, 'FontSize', 32, 'FontWeight', 'bold');% Remove top and right ticks
ax = gca;
ax.XAxis.TickLength = [0.015 0];  % Bottom ticks only
ax.YAxis.TickLength = [0.015 0];  % Left ticks only
ax.Box = 'off';  % Turn off the box
ax.TickLength = [0.015 0];  % Set the tick length

legend boxoff




