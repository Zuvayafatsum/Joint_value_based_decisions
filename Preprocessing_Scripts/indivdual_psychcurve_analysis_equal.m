%Run psy_curves_equal first!!

unique_dyads = unique(meta_data_l_e_c(:,6));


for i = 1:length(unique_dyads);

dyad_idx = find (meta_data_l_e_c(:,6) == i);

temp_data = meta_data_l_e_c(dyad_idx,:)

data = [temp_data(:,3), temp_data(:,2)]

% Unique stimulus intensity levels
intensity_levels = unique(data(:, 1));

% Calculate the probability of correct responses for each intensity level
prob_correct = arrayfun(@(x) mean(data(data(:, 1) == x, 2)), intensity_levels);


bhat = glmfit(intensity_levels,[prob_correct ones(size(prob_correct))],'binomial','link','probit')

estimatedMean = -bhat(1)/bhat(2);
estimatedSD   = 1/bhat(2);
slope_l_e_c =bhat(2)

l_e_c_metrics(i,1) = estimatedMean
l_e_c_metrics (i,2)= estimatedSD

l_e_c_metrics (i,3) = i; %dyad nr
l_e_c_metrics(i,4) = slope_l_e_c


end

unique_dyads = []

%Short horizon
unique_dyads = unique(meta_data_s_e_c(:,6));


for i = 1:length(unique_dyads);

dyad_idx = find (meta_data_s_e_c(:,6) == i);

temp_data = meta_data_s_e_c(dyad_idx,:)

data = [temp_data(:,3), temp_data(:,2)]

% Unique stimulus intensity levels
intensity_levels = unique(data(:, 1));

% Calculate the probability of correct responses for each intensity level
prob_correct = arrayfun(@(x) mean(data(data(:, 1) == x, 2)), intensity_levels);


bhat = glmfit(intensity_levels,[prob_correct ones(size(prob_correct))],'binomial','link','probit')

estimatedMean = -bhat(1)/bhat(2);
estimatedSD   = 1/bhat(2);
slope_s_e_c =bhat(2)

s_e_c_metrics(i,1) = estimatedMean
s_e_c_metrics (i,2)= estimatedSD

s_e_c_metrics (i,3) = i; %dyad nr
s_e_c_metrics(i,4) = slope_s_e_c


end

% Individuals


%Long horizon
unique_dyads = unique(meta_data_l_e_i(:,6));


for i = 1:length(unique_dyads);

dyad_idx = find (meta_data_l_e_i(:,6) == i);


temp_data = meta_data_l_e_i(dyad_idx,:)
ind_1 = find(temp_data(:,5)==1)
ind_2 = find(temp_data(:,5)==2)


temp_data_1 = temp_data(ind_1,:)
temp_data_2 = temp_data(ind_2,:)


data = [temp_data_1(:,3), temp_data_1(:,2),temp_data_2(:,3), temp_data_2(:,2)]

% Unique stimulus intensity levels
intensity_levels = unique(data(:, 1));

% Calculate the probability of correct responses for each intensity level
prob_correct_1 = arrayfun(@(x) mean(data(data(:, 1) == x, 2)), intensity_levels);
prob_correct_2 = arrayfun(@(x) mean(data(data(:, 1) == x, 4)), intensity_levels);


bhat_1 = glmfit(intensity_levels,[prob_correct_1 ones(size(prob_correct_1))],'binomial','link','probit')
bhat_2 = glmfit(intensity_levels,[prob_correct_2 ones(size(prob_correct_2))],'binomial','link','probit')


estimatedMean_1 = -bhat_1(1)/bhat_1(2);
estimatedSD_1   = 1/bhat_1(2);

estimatedMean_2 = -bhat_2(1)/bhat_2(2);
estimatedSD_2   = 1/bhat_2(2);

estimatedMean_avg = mean([estimatedMean_1 estimatedMean_2])
estimatedSD_avg = mean([estimatedSD_1 estimatedSD_2])

slope_l_e_i_1 =bhat_1(2)
slope_l_e_i_2 =bhat_2(2)
slope_avg= mean([slope_l_e_i_1 slope_l_e_i_2])



l_e_i_metrics(i,1) = estimatedMean_1
l_e_i_metrics (i,2)= estimatedSD_1
l_e_i_metrics (i,3)= 1
l_e_i_metrics(i,4) = estimatedMean_2
l_e_i_metrics (i,5)= estimatedSD_2
l_e_i_metrics (i,6)= 2


l_e_i_metrics(i,7) = estimatedMean_avg
l_e_i_metrics(i,8) = estimatedSD_avg



l_e_i_metrics (i,9) = i
l_e_i_metrics(i,10) = slope_l_e_i_1
l_e_i_metrics(i,11) = slope_l_e_i_2
l_e_i_metrics(i,12) = slope_avg

end



%Short horizon
unique_dyads = unique(meta_data_s_e_i(:,6));


for i = 1:length(unique_dyads);

dyad_idx = find (meta_data_s_e_i(:,6) == i);


temp_data = meta_data_s_e_i(dyad_idx,:)
ind_1 = find(temp_data(:,5)==1)
ind_2 = find(temp_data(:,5)==2)


temp_data_1 = temp_data(ind_1,:)
temp_data_2 = temp_data(ind_2,:)


data = [temp_data_1(:,3), temp_data_1(:,2),temp_data_2(:,3), temp_data_2(:,2)]

% Unique stimulus intensity levels
intensity_levels = unique(data(:, 1));

% Calculate the probability of correct responses for each intensity level
prob_correct_1 = arrayfun(@(x) mean(data(data(:, 1) == x, 2)), intensity_levels);
prob_correct_2 = arrayfun(@(x) mean(data(data(:, 1) == x, 4)), intensity_levels);


bhat_1 = glmfit(intensity_levels,[prob_correct_1 ones(size(prob_correct_1))],'binomial','link','probit')
bhat_2 = glmfit(intensity_levels,[prob_correct_2 ones(size(prob_correct_2))],'binomial','link','probit')


estimatedMean_1 = -bhat_1(1)/bhat_1(2);
estimatedSD_1   = 1/bhat_1(2);

estimatedMean_2 = -bhat_2(1)/bhat_2(2);
estimatedSD_2   = 1/bhat_2(2);

estimatedMean_avg = mean([estimatedMean_1 estimatedMean_2])
estimatedSD_avg = mean([estimatedSD_1 estimatedSD_2])
slope_s_e_i_1 =bhat_1(2)
slope_s_e_i_2 =bhat_2(2)
slope_avg= mean([slope_s_e_i_1 slope_s_e_i_2])



s_e_i_metrics(i,1) = estimatedMean_1
s_e_i_metrics (i,2)= estimatedSD_1
s_e_i_metrics (i,3)= 1
s_e_i_metrics(i,4) = estimatedMean_2
s_e_i_metrics (i,5)= estimatedSD_2
s_e_i_metrics (i,6)= 2
s_e_i_metrics(i,7) = estimatedMean_avg
s_e_i_metrics(i,8) = estimatedSD_avg



s_e_i_metrics (i,9) = i; %dyad nr
s_e_i_metrics(i,10) = slope_s_e_i_1
s_e_i_metrics(i,11) = slope_s_e_i_2
s_e_i_metrics(i,12) = slope_avg

end

%Statistics

% T test for short ind vs col estimatedmean
[h, p, ci, stats] = ttest(s_e_i_metrics(:,7), s_e_c_metrics(:,1))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(s_e_c_metrics(:,1))])

%Long ind vs col estimated mean
[h, p, ci, stats] = ttest(l_e_i_metrics(:,7), l_e_c_metrics(:,1))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(l_e_i_metrics(:,7)),mean(l_e_c_metrics(:,1))])


% long ind vs short ind estimated mean
[h, p, ci, stats] = ttest(l_e_i_metrics(:,7), s_e_i_metrics(:,7))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(l_e_i_metrics(:,7)),mean(s_e_i_metrics(:,7))])


% long col vs ind estimated mean

[h, p, ci, stats] = ttest(l_e_c_metrics(:,1), s_e_c_metrics(:,1))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(l_e_c_metrics(:,1)),mean(s_e_c_metrics(:,1))])


%Statistics - SD

% T test for short ind vs col
[h, p, ci, stats] = ttest(s_e_i_metrics(:,8), s_e_c_metrics(:,2))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(s_e_i_metrics(:,8)),mean(s_e_c_metrics(:,2))])




%Long
[h, p, ci, stats] = ttest(l_e_i_metrics(:,8), l_e_c_metrics(:,2))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(l_e_i_metrics(:,8)),mean(l_e_c_metrics(:,2))])


%short vs long ind

[h, p, ci, stats] = ttest(l_e_i_metrics(:,8), s_e_i_metrics(:,8))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(l_e_i_metrics(:,8)),mean(s_e_i_metrics(:,8))])



%short vs long coll
[h, p, ci, stats] = ttest(l_e_c_metrics(:,2), s_e_c_metrics(:,2))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(l_e_c_metrics(:,2)),mean(s_e_c_metrics(:,2))])



%Statistics - slope


% T test for short ind vs short col
[h, p, ci, stats] = ttest(s_e_i_metrics(:,12), s_e_c_metrics(:,4))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(s_e_i_metrics(:,12)),mean(s_e_c_metrics(:,4))])




%Long individual vs long collective
[h, p, ci, stats] = ttest(l_e_i_metrics(:,12), l_e_c_metrics(:,4))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(l_e_i_metrics(:,12)),mean(l_e_c_metrics(:,4))])


%short vs long ind

[h, p, ci, stats] = ttest(l_e_i_metrics(:,12), s_e_i_metrics(:,12))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(l_e_i_metrics(:,12)),mean(s_e_i_metrics(:,12))])



%short vs long coll
[h, p, ci, stats] = ttest(l_e_c_metrics(:,4), s_e_c_metrics(:,4))
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(l_e_c_metrics(:,4)),mean(s_e_c_metrics(:,4))])



% across horizons

short_horizon_general = [s_e_c_metrics(:,4) ; s_e_i_metrics(:,12)]; 


long_horizon_general = [l_e_c_metrics(:,4) ; l_e_i_metrics(:,12)]; 

%short vs long over context

[h, p, ci, stats] = ttest(short_horizon_general, long_horizon_general)
disp(['p-value: ', num2str(p)]);
disp(['Confidence Interval: [', num2str(ci(1)), ', ', num2str(ci(2)), ']']);
disp(['Test Statistic: ', num2str(stats.tstat)]);
disp([mean(short_horizon_general),mean(long_horizon_general)])


%% create csv dataframes that includes slope values for Ind1, Ind2 and the dyad; for both horizons.

% Slope for short horizon

Ind1_rand_short_slope = s_e_i_metrics(:,10);

Ind2_rand_short_slope = s_e_i_metrics(:,11);

Dyad_rand_short_slope = s_e_c_metrics(:,4);


Random_Exploration_Short_Horizon_Slope = table(Ind1_rand_short_slope,Ind2_rand_short_slope,Dyad_rand_short_slope)

writetable(Random_Exploration_Short_Horizon_Slope, "Random_Exploration_Short_Horizon_Slope.csv")

% Slope for long horizon

Ind1_rand_long_slope = l_e_i_metrics(:,10);

Ind2_rand_long_slope = l_e_i_metrics(:,11);

Dyad_rand_long_slope = l_e_c_metrics(:,4);


Random_Exploration_Long_Horizon_Slope = table(Ind1_rand_long_slope,Ind2_rand_long_slope,Dyad_rand_long_slope)

writetable(Random_Exploration_Long_Horizon_Slope, "Random_Exploration_Long_Horizon_Slope.csv")


% combined table


Replication_Random_Exploration_Slope_All =  table(Ind1_rand_short_slope,Ind2_rand_short_slope,Dyad_rand_short_slope,Ind1_rand_long_slope,Ind2_rand_long_slope,Dyad_rand_long_slope)


writetable(Replication_Random_Exploration_Slope_All, "Replication_Random_Exploration_Slope_All.csv")



