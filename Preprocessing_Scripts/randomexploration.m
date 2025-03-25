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


answer = [game.mean]'
trial_type = [game.gameLength]'

%Filter for horizon 5 & 10 only

short_trials = find(trial_type == 5)
long_trials = find(trial_type == 10)
short_only_data = game(short_trials)
long_only_data = game(long_trials)


%Take last element of each array in answers and store them in a separate
%variable, then compare with ground truth

short_trials_answers = {short_only_data.a}.'
long_trials_answers = {long_only_data.a}.'
short_trials_answers = cell2mat(short_trials_answers)
long_trials_answers = cell2mat(long_trials_answers)

forced_short = {short_only_data.nforced}.'
forced_long = {long_only_data.nforced}.'

%Filter for equal trials only

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



% Determines "incorrect" answers (LOWER mean) for short and equal trials
short_only_equal = short_only_data(equal_mask)
answershort = [short_only_equal.mean]'
answer_comparison = answershort %Has mean choices for equal short trials
incorrect_answer = zeros(size(game, 1), 1);
for i = 1:size(answer_comparison, 1)
    if answer_comparison(i, 1) < answer_comparison(i, 2)
        incorrect_answer(i) = 1;
    else
        incorrect_answer(i) = 2;
    end
end
incorrect_answer = incorrect_answer'





%For horizon 5 only - "accuracy" is the probability of choosing lower mean

short_responses = short_trials_answers(:,5)

accuracy_short = zeros(length(short_responses),1);

for i = 1 : length(short_responses);

if short_responses(i) == incorrect_answer(i)
    accuracy_short(i) = 1;
else 
    accuracy_short(i) = 0;
end
end

accuracy_fivehorizon = sum(accuracy_short)/length(accuracy_short)

sanity_check(:,1) = short_responses
sanity_check(:,2) = incorrect_answer
sanity_check(:,3) = accuracy_short % check complete :)




%For horizon 10 only - determine lower mean choices for long and equal
%trials

%Filter for equal trials only - long

condition_unequal_majority = zeros(size(forced_long, 1), 1);
forced_long =cell2mat(forced_long)

for i = 1:size(forced_long, 1)
    countOfOnes = sum(forced_long(i, :) == 1)
    if countOfOnes == 2;
        condition_unequal_majority(i) = 0;
    elseif countOfOnes == 1;
        condition_unequal_majority(i) = 2;
    else
        condition_unequal_majority(i) = 1;
    end
end

equal_mask = find(condition_unequal_majority == 0);
long_only_equal = long_only_data(equal_mask)
long_trials_answers = {long_only_equal.a}.'
long_trials_answers = cell2mat(long_trials_answers)



% Determines "incorrect" answers (LOWER mean) for long and equal trials

answerlong = [long_only_equal.mean]'
long_answer_comparison = answerlong %Has mean choices for equal short trials
long_incorrect_answer = zeros(size(game, 1), 1);
for i = 1:size(long_answer_comparison, 1)
    if long_answer_comparison(i, 1) < long_answer_comparison(i, 2)
        long_incorrect_answer(i) = 1;
    else
        long_incorrect_answer(i) = 2;
    end
end
long_incorrect_answer = long_incorrect_answer'

%Compare choices with lower mean options

long_responses = long_trials_answers(:,5)
accuracy_long = zeros(length(long_responses),1);


for i = 1 : length(long_responses);

if long_responses(i) == long_incorrect_answer(i)
    accuracy_long(i) = 1;
else 
    accuracy_long(i) = 0;
end
end

accuracy_tenhorizon = sum(accuracy_long)/length(accuracy_long)


sanity_check(:,1) = long_responses
sanity_check(:,2) = long_incorrect_answer
sanity_check(:,3) = accuracy_long % check complete :)



gavg_tenhor(t,1) = accuracy_tenhorizon; 
gavg_fivehor(t,1)= accuracy_fivehorizon; 

   data_ID = char(data_ID)

   if data_ID(8) == 'C' | data_ID(9) == 'C' 
       data_group(t,1) = 1
   else
       data_group(t,1) = 0;
   end
   
end

random_exp(:,1) =data_group;
random_exp(:,2) =gavg_fivehor; 
random_exp(:,3) =gavg_tenhor;


%% 
capsize = 24;
linewidth = 8;
dotsize = 256;




col_idx = find(random_exp(:,1) == 1)
ind_idx = find(random_exp(:,1) == 0)

col_data = random_exp(col_idx,:);
ind_data = random_exp(ind_idx,:);


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

errorbar(1,mean(ind_data(:,2),1),eb_fivehor_ind,'Color','black','LineWidth',linewidth,"CapSize",capsize)
errorbar(2,mean(col_data(:,2),1),eb_fivehor_col,'Color','black','LineWidth',linewidth,"CapSize",capsize)
errorbar(1,mean(ind_data(:,3),1),eb_tenhor_ind,'Color',[0.8 0.39 0.039],'LineWidth',linewidth,"CapSize",capsize)
errorbar(2,mean(col_data(:,3),1),eb_tenhor_col,'Color',[0.8 0.39 0.039],'LineWidth',linewidth,"CapSize",capsize)

set(gcf, 'Color', 'w')

xlim([0.5 2.5])
ylim([0.1 0.3])
yticks([0 0.1 0.15 0.2 0.25 0.3])


xticks([1 2])
xticklabels({'',''}) %??
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



