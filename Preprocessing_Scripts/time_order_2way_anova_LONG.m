%% Runs a 2-way mixed factor ANOVA to find main effect of TIME, main effect of ORDER and their interactions
clc
clear all

%Prepare data first for the bigger set = replication study

% Load the binary experimental condition
load('block_type_order_replication.mat'); % block_type_order is a 17x1 column 
% 1 = First block individual then collective
% 2 first block collective then individual

replication_block_type = block_type_order;


% Convert `long_data` to array
long_data_replication = readtable('All_Accuracy_Replication.csv');
long_data_replication = table2array(long_data_replication);

% Calculate individual and dyad accuracies
Ind_Choice1_replication = (long_data_replication(:,1) + long_data_replication(:,2)) / 2;
Ind_Choice2_replication = (long_data_replication(:,4) + long_data_replication(:,5)) / 2;
Ind_Choice3_replication = (long_data_replication(:,7) + long_data_replication(:,8)) / 2;
Ind_Choice4_replication = (long_data_replication(:,10) + long_data_replication(:,11)) / 2;
Ind_Choice5_replication = (long_data_replication(:,13) + long_data_replication(:,14)) / 2;
Ind_Choice6_replication = (long_data_replication(:,16) + long_data_replication(:,17)) / 2;

Ind_Acc_replication = [Ind_Choice1_replication Ind_Choice2_replication Ind_Choice3_replication Ind_Choice4_replication Ind_Choice5_replication Ind_Choice6_replication];
Dyad_Acc_replication = [long_data_replication(:,3) long_data_replication(:,6) long_data_replication(:,9) long_data_replication(:,12) long_data_replication(:,15) long_data_replication(:,18)];

% Combine data into one matrix
data_replication = [mean(Ind_Acc_replication,2), mean(Dyad_Acc_replication,2)];  % Now this is a 17x12 matrix (6 measurements per condition)
% first column IND second column is COL

%% Do the same for pilot/discovery study

%Prepare data first for the bigger set = discovery study

% Load the binary experimental condition
load('block_type_order_discovery.mat'); % block_type_order is a 17x1 column

discovery_block_type = block_type_order_discovery;


% Convert `long_data` to array
long_data_discovery = readtable('Discovery_Long_Acc.csv');
long_data_discovery = table2array(long_data_discovery);

% Calculate individual and dyad accuracies
Ind_Choice1_discovery = (long_data_discovery(:,1) + long_data_discovery(:,2)) / 2;
Ind_Choice2_discovery = (long_data_discovery(:,4) + long_data_discovery(:,5)) / 2;
Ind_Choice3_discovery = (long_data_discovery(:,7) + long_data_discovery(:,8)) / 2;
Ind_Choice4_discovery = (long_data_discovery(:,10) + long_data_discovery(:,11)) / 2;
Ind_Choice5_discovery = (long_data_discovery(:,13) + long_data_discovery(:,14)) / 2;
Ind_Choice6_discovery = (long_data_discovery(:,16) + long_data_discovery(:,17)) / 2;

Ind_Acc_discovery = [Ind_Choice1_discovery Ind_Choice2_discovery Ind_Choice3_discovery Ind_Choice4_discovery Ind_Choice5_discovery Ind_Choice6_discovery];
Dyad_Acc_discovery = [long_data_discovery(:,3) long_data_discovery(:,6) long_data_discovery(:,9) long_data_discovery(:,12) long_data_discovery(:,15) long_data_discovery(:,18)];

% Combine data into one matrix
data_discovery = [mean(Ind_Acc_discovery,2), mean(Dyad_Acc_discovery,2)];  % Now this is a 17x12 matrix (6 measurements per condition)


%% combine datasets for analysis 
data = [data_replication; data_discovery];
block_type_order = [replication_block_type ; discovery_block_type];


%% 
new_datastr = zeros(26,3);

for i = 1:26;

if block_type_order(i,1) == 1;
    new_datastr(i,1) = data(i,1);
    new_datastr(i,2) = data(i,2);
else
    new_datastr(i,1) = data(i,2);
    new_datastr(i,2) = data(i,1);
end
end

new_datastr(:,3) = block_type_order;


% Load your dataset (assuming it's a 26x3 matrix named 'data')
% Column 1: Performance at Time 1
% Column 2: Performance at Time 2
% Column 3: Order (1 = Individual first, 2 = Social first)

% Extract relevant variables
performance = new_datastr(:,1:2); % Repeated measures (within-subjects factor: Time)
order = new_datastr(:,3);         % Between-subjects factor (Order)


% Convert data into a table format suitable for repeated-measures ANOVA
numSubjects = size(data,1);
Subject = (1:numSubjects)'; % Subject IDs
Order = categorical(order); % Convert Order to categorical variable

% Create a table with repeated measures (wide format)
T = table(Subject, Order, performance(:,1), performance(:,2), ...
          'VariableNames', {'Subject', 'Order', 'Time1', 'Time2'});

% Define repeated-measures factor (Time: 1st session vs. 2nd session)
WithinFactor = table([1; 2], 'VariableNames', {'Time'}); 

% Define repeated-measures model
rm = fitrm(T, 'Time1,Time2 ~ Order', 'WithinDesign', WithinFactor);

% Run the repeated-measures ANOVA (for within-subjects factor: Time)
ranovaTable = ranova(rm);
disp('Repeated-Measures ANOVA Results (Main Effect of Time and Interaction):');
disp(ranovaTable);

% Extract and display the main effect of Time
mainEffectTime = ranovaTable(1,:); % The first row corresponds to the main effect of Time
disp('Main Effect of Time:');
disp(mainEffectTime);

% Run the between-subjects ANOVA (for Order effect)
anovaTable = anova(rm);
disp('Between-Subjects ANOVA Results (Main Effect of Order):');
disp(anovaTable);

% reporting og the results
% Extract values for the main effect of Time
timeEffect = ranovaTable(1,:); % First row = Main Effect of Time
F_time = timeEffect.F; % F-value
df_time1 = timeEffect.DF; % Numerator df
df_time2 = ranovaTable.DF(end); % Correct denominator df (from the residual/error term)
p_time = timeEffect.pValue; % p-value
eta2_time = timeEffect.SumSq / sum(ranovaTable.SumSq); % Partial eta squared

% Display in APA format
fprintf('Main effect of Time: F(%d, %d) = %.3f, p = %.3f, ?² = %.3f\n', ...
    df_time1, df_time2, F_time, p_time, eta2_time);

% Extract values for the main effect of Order (Between-Subjects)
orderEffect = anovaTable(2,:); % Second row corresponds to the Order effect
F_order = orderEffect.F; % F-value from the Order effect
df_order1 = orderEffect.DF; % Numerator df (for Order)
df_order2 = anovaTable.DF(3); % Denominator df (from Error term, the third row)
p_order = orderEffect.pValue; % p-value for Order
eta2_order = orderEffect.SumSq / sum(anovaTable.SumSq); % Partial eta squared for Order

% Display in APA format
fprintf('Main effect of Order: F(%d, %d) = %.3f, p = %.3f, ?² = %.3f\n', ...
    df_order1, df_order2, F_order, p_order, eta2_order);


% Extract values for the Interaction (Time × Order)
interactionEffect = ranovaTable(2,:);
F_interaction = interactionEffect.F;
df_interaction1 = interactionEffect.DF; % Numerator df
df_interaction2 = ranovaTable.DF(end); % Correct denominator df
p_interaction = interactionEffect.pValue;
eta2_interaction = interactionEffect.SumSq / sum(ranovaTable.SumSq);

% Display in APA format
fprintf('Interaction effect: F(%d, %d) = %.3f, p = %.3f, ?² = %.3f\n', ...
    df_interaction1, df_interaction2, F_interaction, p_interaction, eta2_interaction)



% Post-hoc tests (if needed)
disp('Post-hoc pairwise comparisons:');
multcompare(rm, 'Time') % Test within-subject Time effect
multcompare(rm, 'Order') % Test between-subject Order effect
multcompare(rm, 'Time', 'By', 'Order') % Test interaction effects

%% post-hoc comparisons

% Run pairwise comparisons using multcompare
comparisons = multcompare(rm, 'Time', 'By', 'Order');

% Extract relevant values from the multcompare output
differences = table2array(comparisons(:, 4)); % Mean difference between groups
stdErr =table2array( comparisons(:, 5)); % Standard error
p_values = table2array(comparisons(:, 6)); % p-values
lower_CI = table2array(comparisons(:, 7)); % Lower bound of the confidence interval
upper_CI = table2array(comparisons(:, 8)); % Upper bound of the confidence interval
time1 = table2array(comparisons(:, 2)); % Time point 1
time2 = table2array(comparisons(:, 3)); % Time point 2

% Calculate t-values based on the difference and standard error
t_values = differences ./ stdErr;

% Degrees of freedom for each comparison (assuming same df as in ANOVA output)
df = 24; % Based on the number of subjects (26), so df = 26 - 2 = 24

% Format the results in APA style
for i = 1:length(t_values)
    fprintf('Comparison %d: t(%d) = %.3f, p = %.3f, CI = [%.3f, %.3f]\n', ...
        i, df, t_values(i), p_values(i), lower_CI(i), upper_CI(i));
end




%% further inspect interaction effect: compare ind vs collective in time 1 and time 2
% for time1

ind_time1 = find(new_datastr(:,3) == 1);
col_time1 = find(new_datastr(:,3) == 2);

time1_ind_acc = new_datastr(ind_time1,1);
time1_col_acc = new_datastr(col_time1,1);

[h,p,ci,stats] = ttest2(time1_ind_acc,time1_col_acc)

% for time 2
ind_time2 = find(new_datastr(:,3) == 2);
col_time2 = find(new_datastr(:,3) == 1);

time2_ind_acc = new_datastr(ind_time2,2);
time2_col_acc = new_datastr(col_time2,2);

[h,p,ci,stats] = ttest2(time2_ind_acc,time2_col_acc)

% Calculate means and standard deviations
mean1 = mean(time2_ind_acc);
std1 = std(time2_ind_acc);
mean2 = mean(time2_col_acc);
std2 = std(time2_col_acc);

% Display results in a structured format
fprintf('Group\tMean\tSD\tn\nt(df)\tp\n');
fprintf('Group 1\t%.3f\t%.3f\t%d\n', mean1, std1, length(time2_ind_acc));
fprintf('Group 2\t%.3f\t%.3f\t%d\n', mean2, std2, length(time2_col_acc));
fprintf('t(%d) = %.3f, p = %.4f\n', stats.df, stats.tstat, p);

%% Visualize the effects

%% visualize the effect for the average of long horizon choice sequence

% Sample data (replace these with your actual data)
accuracy = data;  % 26x2 matrix (accuracy percentages)



% Extracting data based on order
indiv_first = block_type_order == 1;  % Subjects who did individual first
soc_first = block_type_order == 2;    % Subjects who did social first

% Compute means and standard errors
mean_1st_session_indiv = mean(accuracy(indiv_first,1)); % Individual first session
mean_2nd_session_soc = mean(accuracy(indiv_first,2));   % Social second session

mean_1st_session_soc = mean(accuracy(soc_first,2)); % Social first session
mean_2nd_session_indiv = mean(accuracy(soc_first,1)); % Individual second session

se_1st_session_indiv = std(accuracy(indiv_first,1)) / sqrt(sum(indiv_first));
se_2nd_session_soc = std(accuracy(indiv_first,2)) / sqrt(sum(indiv_first));

se_1st_session_soc = std(accuracy(soc_first,2)) / sqrt(sum(soc_first));
se_2nd_session_indiv = std(accuracy(soc_first,1)) / sqrt(sum(soc_first));

% X positions for plot (1st and 2nd session)
x_positions = [1, 2];
x_positions_jittered = [1.05, 1.95];

% Define colors
indiv_color = [0, 0.5, 0];  % Dark green (Individual)
soc_color = [0.5, 0, 0.5];  % Purple (Social)

% Create figure
figure; hold on;
axis square
set(gcf, 'Color', 'w')
set(gca, 'TickDir', 'out')

ax = gca;
ax.XAxis.TickLength = [0.015 0];  % Bottom ticks only
ax.YAxis.TickLength = [0.015 0];  % Left ticks only
ax.Box = 'off';  % Turn off the box
ax.TickLength = [0.015 0];  % Set the tick length
set(gca, 'FontSize', 24, "FontWeight","bold"); % Change the font size of axis labels
set(gca, 'LineWidth', 4); % Change the width of the axis lines

% Connect Individual (1st session) ? Social (2nd session) [Indiv First Group]
plot(x_positions, 100*[mean_1st_session_indiv, mean_2nd_session_soc], '-o', ...
    'Color', [0 0 0 0.5], 'MarkerFaceColor', 'k', ...
   'LineWidth', 2, 'MarkerSize', 12);

% Connect Social (1st session) ? Individual (2nd session) [Soc First Group]
plot(x_positions_jittered, 100*[mean_1st_session_soc, mean_2nd_session_indiv], '-o', ...
    'Color', [0 0 0 0.5], 'MarkerFaceColor', 'k', ...
  'LineWidth', 2, 'MarkerSize', 12);

% Add error bars
errorbar(x_positions, 100*[mean_1st_session_indiv, mean_2nd_session_soc], ...
   100*[se_1st_session_indiv, se_2nd_session_soc], 'k', 'LineStyle', 'none', 'LineWidth', 2);

errorbar(x_positions_jittered, 100*[mean_1st_session_soc, mean_2nd_session_indiv], ...
    100*[se_1st_session_soc, se_2nd_session_indiv], 'k', 'LineStyle', 'none', 'LineWidth', 2);

% Customize plot
xticks([1 2]);
xticklabels({'', ''});
ylabel('');
axis square
set(gca, 'FontSize', 24, 'LineWidth', 2);
xlim([0.8 2.2]);
ylim([75 100])
hold off;

