
clc
clear all

directory = 'D:\Program Files\MATLAB\Joint_Horizon\main_study\Mixed\';
data_files = dir('D:\Program Files\MATLAB\Joint_Horizon\main_study\Mixed\');

data_file_names = {data_files.name}';
data_file_names(1:2,:) = [];


%% 

gavg_tenhor = []; 
gavg_fivehor = []; 
data_type = {}; 

for t = 1: length(data_file_names);

    data_ID = string(data_file_names(t));

    
   load('D:\Program Files\MATLAB\Joint_Horizon\main_study\Mixed\' + data_ID);



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

%Horizon 5 - "accuracy" is probability of choosing more informative option

condition_unequal_majority = zeros(size(forced_short, 1), 1);
more_informative = zeros(size(forced_short, 1), 1);
forced_short =cell2mat(forced_short)

for i = 1:size(forced_short, 1)
    countOfOnes = sum(forced_short(i, :) == 1)
    if countOfOnes == 2;
        condition_unequal_majority(i) = 0;
        more_informative_short(i) = 0;
    elseif countOfOnes == 1;
        condition_unequal_majority(i) = 2;
        more_informative_short(i)=1
    else
        condition_unequal_majority(i) = 1;
        more_informative_short(i)=2;
    end
end

more_informative_short = more_informative_short'

unequal_mask = find(more_informative_short ~= 0);
unequal_trials_only_short = short_only_data(unequal_mask)
higher_info_choice_only_unequal = more_informative_short(unequal_mask)
higher_info_choice_only_unequal = higher_info_choice_only_unequal'
%unequal_trials_only_long = long_only_data(unequal_mask)






%Take last element of each array in answers and store them in a separate
%variable, then compare with ground truth

short_trials_answers = {unequal_trials_only_short.a}.'
%long_trials_answers = {unequal_trials_only_long.a}.'
short_trials_answers = cell2mat(short_trials_answers)
%long_trials_answers = cell2mat(long_trials_answers)



%For horizon 5 only

short_responses = short_trials_answers(:,5)

accuracy_short = zeros(length(short_responses),1);

higher_info_choice_only_unequal = higher_info_choice_only_unequal'
for i = 1 : length(short_responses);

if short_responses(i) == higher_info_choice_only_unequal(i)
    accuracy_short(i) = 1;
else 
    accuracy_short(i) = 0;
end
end

accuracy_fivehorizon = sum(accuracy_short)/length(accuracy_short)

% Sanity check - completed
sanitycheck(:,1) = short_responses
sanitycheck(:,2) = higher_info_choice_only_unequal
sanitycheck(:,3) = accuracy_short


%For long trials -  only first free choice!!

% Remove equal trials

condition_unequal_majority = zeros(size(forced_long, 1), 1);
more_informative = zeros(size(forced_long, 1), 1);
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

unequal_mask = find(more_informative_long ~= 0);
unequal_trials_only_long = long_only_data(unequal_mask)
higher_info_choice_only_unequal = more_informative_long(unequal_mask)
higher_info_choice_only_unequal = higher_info_choice_only_unequal'
unequal_trials_only_long = long_only_data(unequal_mask)




%Take last element of each array in answers and store them in a separate
%variable, then compare with ground truth

long_trials_answers = {unequal_trials_only_long.a}.'
long_trials_answers = cell2mat(long_trials_answers)



%For horizon 5 only

long_responses = long_trials_answers(:,5)

accuracy_long = zeros(length(long_responses),1);

higher_info_choice_only_unequal = higher_info_choice_only_unequal'
for i = 1 : length(long_responses);

if long_responses(i) == higher_info_choice_only_unequal(i)
    accuracy_long(i) = 1;
else 
    accuracy_long(i) = 0;
end
end

accuracy_tenhorizon = sum(accuracy_long)/length(accuracy_long)

%Sanity check - completed

sanitycheck(:,1) = long_responses
sanitycheck(:,2) = higher_info_choice_only_unequal
sanitycheck(:,3) = accuracy_long


gavg_tenhor(t,1) = accuracy_tenhorizon; 
gavg_fivehor(t,1)= accuracy_fivehorizon; 

   data_ID = char(data_ID)

   if data_ID(8) == 'C' | data_ID(9) == 'C' 
       data_group(t,1) = 1
   else
       data_group(t,1) = 0;
   end
   
end


%% create df for common use 
t = repmat(1:17, 3, 1);
t = t(:);
t(1:3:end) = 0;
% Data group 1 is collective, 0 is individual
directed_exp(:,1) =data_group;
directed_exp(:,2) =gavg_fivehor; 
directed_exp(:,3) =gavg_tenhor;
directed_exp(:,4) = t;
directed_exp(:,5) = 1:length(directed_exp(:,1));




colnames = {'Dyad_or_Ind' 'Directed_Short' 'Directed_Long' 'Dyad_ID','Ind_ID'};
Replication_Directed_Exploration = table(directed_exp(:,1), directed_exp(:,2),directed_exp(:,3),directed_exp(:,4),directed_exp(:,5), 'VariableNames', colnames);



% Initialize the result vector to store the averages for unique numbers
unique_numbers = unique(t(t > 0));  % Extract the unique non-zero numbers (1 to 17)
short_directed_ind_mean = zeros(size(unique_numbers));  % Preallocate for speed


% Loop through each unique number (1 to 17)
for i = 1:length(unique_numbers)
    num = unique_numbers(i);
    
    % Find the indices where 't' equals the current number
    indices = find(t == num);
    
    % Get the corresponding values from the 'values' vector
    corresponding_values = directed_exp(indices,2);
    
    % Calculate the average of these corresponding values
    short_directed_ind_mean(i) = mean(corresponding_values);
end


long_directed_ind_mean = zeros(size(unique_numbers));  % Preallocate for speed


% Loop through each unique number (1 to 17)
for i = 1:length(unique_numbers)
    num = unique_numbers(i);
    
    % Find the indices where 't' equals the current number
    indices = find(t == num);
    
    % Get the corresponding values from the 'values' vector
    corresponding_values = directed_exp(indices,3);
    
    % Calculate the average of these corresponding values
    long_directed_ind_mean(i) = mean(corresponding_values);
end



col_idx = find(directed_exp(:,1) == 1);
short_directed_col = directed_exp(col_idx,2)

long_directed_col =  directed_exp(col_idx,3)


% long vs short for individuals

[h, p, ci, stats] = ttest(short_directed_ind_mean, long_directed_ind_mean)
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(short_directed_ind_mean),mean(long_directed_ind_mean)])


% long vs short for collectives 

[h, p, ci, stats] = ttest(short_directed_col, long_directed_col)
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(short_directed_col),mean(long_directed_col)])


ind_idx = find(directed_exp(:,1) == 0); 

short_ind_all = directed_exp(ind_idx,2)
long_ind_all =directed_exp(ind_idx,3)

[h, p, ci, stats] = ttest(short_ind_all, long_ind_all)
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(short_ind_all),mean(long_ind_all)])


Summ_Replication_Directed_Exploration = table(short_directed_ind_mean,long_directed_ind_mean,short_directed_col,long_directed_col)


% pooled testing 

pooled_short = [short_directed_ind_mean; short_directed_col]

pooled_long = [long_directed_ind_mean; long_directed_col]


[h, p, ci, stats] = ttest(pooled_short, pooled_long)
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(pooled_short),mean(pooled_long)])




colnames = {'short_directed_ind' 'short_directed_col' 'long_directed_ind' 'long_directed_col'};
Directed_Exploration_AVG = table(short_directed_ind_mean, short_directed_col,long_directed_ind_mean,long_directed_col, 'VariableNames', colnames);




%% PLOTTING

capsize = 24;
linewidth = 8;
dotsize = 256;

col_idx = find(directed_exp(:,1) == 1)
ind_idx = find(directed_exp(:,1) == 0)

col_data = directed_exp(col_idx,:);
ind_data = directed_exp(ind_idx,:);

scatter(1:2,[mean(ind_data(:,2),1) mean(col_data(:,2),1)],dotsize,"black","filled")
hold on
scatter(1:2,[mean(ind_data(:,3),1) mean(col_data(:,3),1)],dotsize,[0.8 0.39 0.039],"filled")

set(gca, 'FontSize', 36, "FontWeight","bold"); % Change the font size of axis labels
set(gca, 'LineWidth', 8); % Change the width of the axis lines


%Error bars
eb_fivehor_ind = std(ind_data(:,2),1) / sqrt(length(ind_data(:,2)));
eb_fivehor_col = std(col_data(:,2),1) / sqrt(length(col_data(:,2)));
eb_tenhor_ind= std(ind_data(:,3),1) / sqrt(length(ind_data(:,3)));
eb_tenhor_col = std(col_data(:,2),1) / sqrt(length(col_data(:,2)));

errorbar(1,mean(ind_data(:,2),1),eb_fivehor_ind,'Color','black','LineWidth',linewidth,'CapSize',capsize)
errorbar(2,mean(col_data(:,2),1),eb_fivehor_col,'Color','black','LineWidth',linewidth,'CapSize',capsize)
errorbar(1,mean(ind_data(:,3),1),eb_tenhor_ind,'Color',[0.8 0.39 0.039],'LineWidth',linewidth,'CapSize',capsize)
errorbar(2,mean(col_data(:,3),1),eb_tenhor_col,'Color',[0.8 0.39 0.039],'LineWidth',linewidth,'CapSize',capsize)


set(gcf, 'Color', 'w')


%General stuff
xlim([0.5 2.5])
ylim([0.35 0.7])
xticks([1 2])
xticklabels({'',''}) %??

yticks([0.4 0.5 0.6 0.7])


ylabel("")
xlabel("")
title("")
ax = gca;
ax.XAxis.TickLength = [0.015 0];  % Bottom ticks only
ax.YAxis.TickLength = [0.015 0];  % Left ticks only
ax.Box = 'off';  % Turn off the box
ax.TickLength = [0.015 0];  % Set the tick length
axis square
%%

capsize = 24;
linewidth = 8;
dotsize = 256;

col_idx = find(directed_exp(:,1) == 1);
ind_idx = find(directed_exp(:,1) == 0);

col_data = directed_exp(col_idx,:);
ind_data = directed_exp(ind_idx,:);

% Means for the data points
mean_ind_2 = mean(ind_data(:,2),1);
mean_col_2 = mean(col_data(:,2),1);
mean_ind_3 = mean(ind_data(:,3),1);
mean_col_3 = mean(col_data(:,3),1);

% Scatter plot with dots
scatter(1:2, [mean_ind_2 mean_col_2], dotsize, "black", "filled");
hold on;
scatter(1:2, [mean_ind_3 mean_col_3], dotsize, [0.8 0.39 0.039], "filled");

% Connect the dots with lines in the same color
plot(1:2, [mean_ind_2 mean_col_2], 'Color', 'black', 'LineWidth', linewidth); % Line for black dots
plot(1:2, [mean_ind_3 mean_col_3], 'Color', [0.8 0.39 0.039], 'LineWidth', linewidth); % Line for colored dots

set(gca, 'FontSize', 36, "FontWeight", "bold"); % Change the font size of axis labels
set(gca, 'LineWidth', 8); % Change the width of the axis lines

% Error bars
eb_fivehor_ind = std(ind_data(:,2),1) / sqrt(length(ind_data(:,2)));
eb_fivehor_col = std(col_data(:,2),1) / sqrt(length(col_data(:,2)));
eb_tenhor_ind = std(ind_data(:,3),1) / sqrt(length(ind_data(:,3)));
eb_tenhor_col = std(col_data(:,2),1) / sqrt(length(col_data(:,2)));

errorbar(1, mean(ind_data(:,2),1), eb_fivehor_ind, 'Color', 'black', 'LineWidth', linewidth, 'CapSize', capsize);
errorbar(2, mean(col_data(:,2),1), eb_fivehor_col, 'Color', 'black', 'LineWidth', linewidth, 'CapSize', capsize);
errorbar(1, mean(ind_data(:,3),1), eb_tenhor_ind, 'Color', [0.8 0.39 0.039], 'LineWidth', linewidth, 'CapSize', capsize);
errorbar(2, mean(col_data(:,3),1), eb_tenhor_col, 'Color', [0.8 0.39 0.039], 'LineWidth', linewidth, 'CapSize', capsize);

set(gcf, 'Color', 'w');

% General settings
xlim([0.5 2.5]);
ylim([0.35 0.7]);
xticks([1 2]);
xticklabels({'', ''}); % Empty xticks

yticks([0.4 0.5 0.6 0.7]);

ylabel("");
xlabel("");
title("");
ax = gca;
ax.XAxis.TickLength = [0.015 0];  % Bottom ticks only
ax.YAxis.TickLength = [0.015 0];  % Left ticks only
ax.Box = 'off';  % Turn off the box
ax.TickLength = [0.015 0];  % Set the tick length
axis square;

