directory = 'D:\Program Files\MATLAB\Joint_Horizon\pilot_data\Individual\';
data_files = dir('D:\Program Files\MATLAB\Joint_Horizon\pilot_data\Individual\');

data_file_names = {data_files.name}';
data_file_names(1:2,:) = [];



%% 

grandavg_short_ind = []; 
grandavg_long_ind = []; 

for t = 1: length(data_file_names);

    data_ID = string(data_file_names(t));

   load('D:\Program Files\MATLAB\Joint_Horizon\pilot_data\Individual\' + data_ID);



answer = [game.mean]'
trial_type = [game.gameLength]'
short_trials_answers = [game.a]'


% Determines correct answers (higher mean) for all trials

answer_comparison = answer
correct_answer = zeros(size(game, 1), 1);
for i = 1:size(answer_comparison, 1)
    if answer_comparison(i, 1) > answer_comparison(i, 2)
        correct_answer(i) = 1;
    else
        correct_answer(i) = 2;
    end
end
correct_answer = correct_answer'

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

%For horizon 5 only

short_responses = short_trials_answers(:,5)

correct_short = correct_answer(short_trials)

accuracy_short = zeros(length(short_trials),1);

for i = 1 : length(short_trials);

if short_responses(i) == correct_short(i)
    accuracy_short(i) = 1;
else 
    accuracy_short(i) = 0;
end
end

accuracy_fivehorizon = sum(accuracy_short)/length(accuracy_short)



%For horizon 10 only

long_responses = long_trials_answers(:,5:10)

correct_long = correct_answer(long_trials)

accuracy_long = zeros(length(long_trials),6)

for k = 1:6

    for i = 1: length(long_trials);

if long_responses(i,k) == correct_long(i)
    accuracy_long(i,k) = 1;
else 
    accuracy_long(i,k) = 0;    
end

    end
end

accuracy_tenhorizon = sum(accuracy_long)/length(accuracy_long)




grandavg_short_ind(t,1) = accuracy_fivehorizon;
grandavg_long_ind(t,1:6) = accuracy_tenhorizon;


accuracy_tenhorizon = [];
accuracy_fivehorizon = [];
end

short_se_ind = std(grandavg_short_ind,1) / sqrt(length(grandavg_short_ind));
long_se_ind = std(grandavg_long_ind,1) / sqrt(length(grandavg_long_ind));



%% Calculate accuracy over trials for collective conditions

directory = 'D:\Program Files\MATLAB\Joint_Horizon\pilot_data\Collective\';
data_files = dir('D:\Program Files\MATLAB\Joint_Horizon\pilot_data\Collective\');

data_file_names = {data_files.name}';
data_file_names(1:2,:) = [];

%% Collective accuracy 

grandavg_short_col = []; 
grandavg_long_col = []; 

for t = 1: length(data_file_names);

    data_ID = string(data_file_names(t));

   load('D:\Program Files\MATLAB\Joint_Horizon\pilot_data\Collective\' + data_ID);



answer = [game.mean]'
trial_type = [game.gameLength]'
short_trials_answers = [game.a]'


% Determines correct answers (higher mean) for all trials

answer_comparison = answer
correct_answer = zeros(size(game, 1), 1);
for i = 1:size(answer_comparison, 1)
    if answer_comparison(i, 1) > answer_comparison(i, 2)
        correct_answer(i) = 1;
    else
        correct_answer(i) = 2;
    end
end
correct_answer = correct_answer'

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

%For horizon 5 only

short_responses = short_trials_answers(:,5)

correct_short = correct_answer(short_trials)

accuracy_short = zeros(length(short_trials),1);

for i = 1 : length(short_trials);

if short_responses(i) == correct_short(i)
    accuracy_short(i) = 1;
else 
    accuracy_short(i) = 0;
end
end

accuracy_fivehorizon = sum(accuracy_short)/length(accuracy_short)



%For horizon 10 only

long_responses = long_trials_answers(:,5:10)

correct_long = correct_answer(long_trials)

accuracy_long = zeros(length(long_trials),6)

for k = 1:6

    for i = 1: length(long_trials);

if long_responses(i,k) == correct_long(i)
    accuracy_long(i,k) = 1;
else 
    accuracy_long(i,k) = 0;    
end

    end
end

accuracy_tenhorizon = sum(accuracy_long)/length(accuracy_long)




grandavg_short_col(t,1) = accuracy_fivehorizon;
grandavg_long_col(t,1:6) = accuracy_tenhorizon;


accuracy_tenhorizon = [];
accuracy_fivehorizon = [];
end

short_se_col = std(grandavg_short_col,1) / sqrt(length(grandavg_short_col));
long_se_col = std(grandavg_long_col,1) / sqrt(length(grandavg_long_col));




%% Plot the Combıned Graphs

green = [39 83 22 ]/ 256;
purple = [112, 48, 160]/256;
capsize = 24;
linewidth = 8;
dotsize = 256;


% Plot a la Fig1C of the Horizon Task paper
figure(1)



tenplot_col = scatter(1.2:6.5, mean(grandavg_long_col,1),180,purple,"filled",'DisplayName', 'Collective')
hold on
fiveplot_col = scatter(1.2,mean(grandavg_short_col,1),dotsize,purple,"filled")
errorbar(1.2:6.5,mean(grandavg_long_col),long_se_col,'Color',purple,'LineWidth',4,"CapSize",18)
errorbar(1.2,mean(grandavg_short_col),short_se_col,'Color',purple,'LineWidth',linewidth,"CapSize",capsize)


hold on

set(gca, 'FontSize', 36, "FontWeight","bold"); % Change the font size of axis labels
set(gca, 'LineWidth', 8); % Change the width of the axis lines






tenplot_ind = scatter(1:6, mean(grandavg_long_ind,1),180,green,"filled",'DisplayName','Individual')
hold on
fiveplot_ind = scatter(1,mean(grandavg_short_ind,1),dotsize, green,"filled")
errorbar(1:6,mean(grandavg_long_ind),long_se_ind,'Color',green,'LineWidth',4,"CapSize",18)
errorbar(1,mean(grandavg_short_ind),short_se_ind,'Color',green,'LineWidth',linewidth,"CapSize",capsize)


xlabel('')
ylabel('')
title('')
set(gcf, 'Color', 'w')
set(gca, 'TickDir', 'out')

ax = gca;
ax.XAxis.TickLength = [0.015 0];  % Bottom ticks only
ax.YAxis.TickLength = [0.015 0];  % Left ticks only
ax.Box = 'off';  % Turn off the box
ax.TickLength = [0.015 0];  % Set the tick length

axis square

ylim([0.6 1])
xlim([0 6.5])

yticks(0.5:0.1:1)
xticks(1:1:6)

ylim([0.6 0.95])
xlim([0 6.5])


%% Create csv dataframes 

% Short horizon accuracies:

Ind1_idx = 1:2:size(grandavg_long_ind,1);

Ind1_Short_Acc = grandavg_short_ind(Ind1_idx)


Ind2_idx = 2:2:size(grandavg_long_ind,1);

Ind2_Short_Acc = grandavg_short_ind(Ind2_idx)


Dyad_Short_Acc = grandavg_short_col;


Discovery_Short_Acc = table(Ind1_Short_Acc,Ind2_Short_Acc,Dyad_Short_Acc) 


writetable(Discovery_Short_Acc,"Discovery_Short_Acc.csv")


% Long horion accuracies


Ind1_idx = 1:2:size(grandavg_long_ind,1);

Ind1_Long_Acc_Choice1 = grandavg_long_ind(Ind1_idx,1)
Ind1_Long_Acc_Choice2 = grandavg_long_ind(Ind1_idx,2)
Ind1_Long_Acc_Choice3 = grandavg_long_ind(Ind1_idx,3)
Ind1_Long_Acc_Choice4 = grandavg_long_ind(Ind1_idx,4)
Ind1_Long_Acc_Choice5 = grandavg_long_ind(Ind1_idx,5)
Ind1_Long_Acc_Choice6 = grandavg_long_ind(Ind1_idx,6)



Ind2_idx = 2:2:size(grandavg_long_ind,1);

Ind2_Long_Acc_Choice1 = grandavg_long_ind(Ind2_idx,1)
Ind2_Long_Acc_Choice2 = grandavg_long_ind(Ind2_idx,2)
Ind2_Long_Acc_Choice3 = grandavg_long_ind(Ind2_idx,3)
Ind2_Long_Acc_Choice4 = grandavg_long_ind(Ind2_idx,4)
Ind2_Long_Acc_Choice5 = grandavg_long_ind(Ind2_idx,5)
Ind2_Long_Acc_Choice6 = grandavg_long_ind(Ind2_idx,6)




Dyad_Long_Acc_Choice1= grandavg_long_col(:,1);
Dyad_Long_Acc_Choice2= grandavg_long_col(:,2);
Dyad_Long_Acc_Choice3= grandavg_long_col(:,3);
Dyad_Long_Acc_Choice4= grandavg_long_col(:,4);
Dyad_Long_Acc_Choice5= grandavg_long_col(:,5);
Dyad_Long_Acc_Choice6= grandavg_long_col(:,6);



Discovery_Long_Acc = table(Ind1_Long_Acc_Choice1,Ind2_Long_Acc_Choice1,Dyad_Long_Acc_Choice1,...
    Ind1_Long_Acc_Choice2,Ind2_Long_Acc_Choice2,Dyad_Long_Acc_Choice2,...
    Ind1_Long_Acc_Choice3,Ind2_Long_Acc_Choice3,Dyad_Long_Acc_Choice3,...
    Ind1_Long_Acc_Choice4,Ind2_Long_Acc_Choice4,Dyad_Long_Acc_Choice4,...
    Ind1_Long_Acc_Choice5,Ind2_Long_Acc_Choice5,Dyad_Long_Acc_Choice5,...
    Ind1_Long_Acc_Choice6,Ind2_Long_Acc_Choice6,Dyad_Long_Acc_Choice6)


writetable(Discovery_Long_Acc,"Discovery_Long_Acc.csv")



