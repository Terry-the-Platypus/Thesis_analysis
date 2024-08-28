clear all
%Specify the folder containing the .xls files to be used for the analysis
folder_path = 'C:\Users\Teren\OneDrive\Bureaublad\RM CN yr2\Internship and Thesis\analysis\';

%Get list of the files in the folder, used .csv as that is how its stored
files = dir(fullfile(folder_path, '*.csv'));

for a = 1:length(files)
    %Construct full file path of curent xls file
    file_path = fullfile(folder_path, files(a).name);

    %Read data from xls file into a table
    data = readtable(file_path);

    %Store data table along with file name
    data_tables{a, 1} = files(a).name;
    data_tables{a, 2} = data;
end



%Initialize a cell array to store tables and thier coresponding file names,
%as well as removed data
filtered_data_tables = cell(length(files), 4); % 4 columns {file name, filtered data}
removed_data_tables = cell(length(files), 5); % 5 columns {file name, removed data}

%Loop through each table stored in second column of data_table
for a = 1:size(data_tables, 1)
    %Retreive the table from second column
    table_data = data_tables{a, 2};

    %Find the correct rejections for each file
    correct_rejections = table_data(isnan(table_data.RT) & table_data.target == 0, :);

    %Find the missed targets for each file
    missed_targets = table_data(isnan(table_data.RT) & table_data.target == 1, :);

    %Filter the data for RTs within the 200 to 1000 ms time window
    time_filtered_data = table_data(table_data.RT >= 200 & table_data.RT <=1000, :);

    %Store the removed data based on the time filter above
    RT_too_fast = table_data(table_data.RT <= 200, :);
    RT_too_slow = table_data(table_data.RT >= 1000, :);

    %Select correct presses with RT values within IQR
    %Calc the IQR for the RT column
    rt_column = time_filtered_data.RT;
    iqr_rt = iqr(rt_column);

    %Compute lower and upper bounds based on IQR
    lower_bound = quantile(rt_column, 0.25) - 1.5 * iqr_rt;
    upper_bound = quantile(rt_column, 0.75) + 1.5 * iqr_rt;

    %Filter data to inlcude only RTs w/in IQR
    IQR_filtered_data= time_filtered_data(rt_column >= lower_bound & rt_column <= upper_bound, :);

    %Create table containing the lower and upper bounds of the RT data
    Q1_IQR_data = time_filtered_data(rt_column < lower_bound, :);
    Q4_IQR_data = time_filtered_data(rt_column > upper_bound, :);
    Q1_and_Q4_IQR = vertcat(Q1_IQR_data, Q4_IQR_data);

    %Filter the data table such that we get all the correct button presses
    %made
    correct_presses = IQR_filtered_data(~isnan(IQR_filtered_data.RT) & IQR_filtered_data.target == 1, :);

    %Find the false alarms (button pressed with no target) for each file
    false_alarms = IQR_filtered_data(~isnan(IQR_filtered_data.RT) & IQR_filtered_data.target == 0, :);


    %Store filtered data table along with file name
    filtered_data_tables{a, 1} = files(a).name;
    filtered_data_tables{a, 2} = correct_rejections; %Correct Rejections
    filtered_data_tables{a, 3} = missed_targets; %Misses
    filtered_data_tables{a, 4} = time_filtered_data;
    filtered_data_tables{a, 5} = IQR_filtered_data;
    filtered_data_tables{a, 6} = correct_presses; %Hits
    filtered_data_tables{a, 7} = false_alarms; %FAs

    %Store th remoed data along with file name
    removed_data_tables{a, 1} = files(a).name;
    removed_data_tables{a, 2} = RT_too_slow;
    removed_data_tables{a, 3} = RT_too_fast;
    removed_data_tables{a, 4} = Q1_and_Q4_IQR;
   end






%Initialize array to store the d' scores in
trial_types = {};
d_prime_scores = {};

%Calculate the d' value for each of the files
%Loop through each file
for b = 1:size(filtered_data_tables, 1)
    %Retreive the number of hits, misses, false alarms, and correct
    %rejections
    hits = size(filtered_data_tables{b, 6});
    misses = size(filtered_data_tables{b, 3});
    FAs = size(filtered_data_tables{b, 7});
    cor_rejections = size(filtered_data_tables{b, 2});

    %Calculate the proportions needed
    signal_trials = hits(1) + misses(1);
    noise_trials = FAs(1) + cor_rejections(1);
    total_trials = signal_trials + noise_trials;

    %Calculate the proportions of signal and noise trials out of the total
    %trials in order to apply the loglinear method
    prop_signal = signal_trials/total_trials;
    prop_noise = 1 - prop_signal;

    %Store the number of signal and noise trials
    trial_types{b, 1} = filtered_data_tables{b, 1};
    trial_types{b, 2} = signal_trials;
    trial_types{b, 3} = noise_trials;
    trial_types{b, 4} = total_trials;
    trial_types{b, 5} = prop_signal;
    trial_types{b, 6} = prop_noise;

    %Proportions of hits and false alarms, with the additions based on the
    %loglinear method (Hautus, 1995)
    p_hits = (hits(1) + prop_signal) / (signal_trials + (2 * prop_signal));
    p_false_alarms = (FAs(1) + prop_noise) / (noise_trials + (2 * prop_noise));

    %Calculate z-scores
    z_hits = norminv(p_hits);
    z_false_alarms = norminv(p_false_alarms);

    %Calculate d'
    d_prime = z_hits - z_false_alarms;

    d_prime_scores{b, 1} = filtered_data_tables{b, 1};
    d_prime_scores{b, 2} = d_prime;
end

%Save the d prime scores as an xls file
filename = 'D_prime_scores.xls';
writecell(d_prime_scores, filename);


%Split the data based on the different conditions
%Sort files as TMS or SHAM stimulation
%Read TMSorSHAM file
TMSorSHAM_file = dir(fullfile(folder_path, "TMSorSHAM.xlsx"));
TMSorSHAM_file_path = fullfile(folder_path, TMSorSHAM_file.name);

%Load order from TMSorSHAM xls file
TMSorSHAM_data = readtable(TMSorSHAM_file_path);

% Initialize cell arrays to store TMS and SHAM data
tms_pre_A_d_prime_data = {};
tms_post_A_d_prime_data = {};
tms_pre_N_d_prime_data = {};
tms_post_N_d_prime_data = {};
sham_pre_A_d_prime_data = {};
sham_post_A_d_prime_data = {};
sham_pre_N_d_prime_data = {};
sham_post_N_d_prime_data = {};

% Loop through each file name in filtered_data_tables
for b = 1:size(d_prime_scores, 1)
    % Extract information from the current file name
    file_name = d_prime_scores{b, 1};
    file_name_parts = strsplit(file_name, '_');
    PP_number = str2double(file_name_parts{1}(3:end)); % Extract the PP number from the file name
    session = str2double(file_name_parts{2}); % Extract session number from file name
    pre_or_post = str2double(file_name_parts{3}); % Extract pre or post stimulation task
    face_type = file_name_parts{5}(1); % Extract if angry or neutral faces presented
    
    % Find matching row in TMSorSHAM_data based on x value
    match_row = find(TMSorSHAM_data{:, 1} == PP_number, 1);
    % Check if a matching row is found
    if ~isempty(match_row)
        % Determine y values based on conditions from TMSorSHAM_data
        if TMSorSHAM_data{match_row, 2} == 1
            % If condition is 1, select y = 1 as TMS and y = 2 as SHAM
            tms_session = 1;
            sham_session = 2;
        else
            % If condition is 0, select y = 1 as SHAM and y = 2 as TMS
            tms_session = 2;
            sham_session = 1;
        end
        % Determine if file should be stored as angry or neutral face in
        % the pre-stimulation or post-stimulation for either TMS or SHAM
        % stimulation
        if pre_or_post == 1
            if session == tms_session
                if face_type == 'A'
                    tms_pre_A_d_prime_data = [tms_pre_A_d_prime_data; d_prime_scores(b, :)];
                elseif face_type == 'N'
                    tms_pre_N_d_prime_data = [tms_pre_N_d_prime_data; d_prime_scores(b, :)];
                end
            elseif session == sham_session
                if face_type =='A'
                    sham_pre_A_d_prime_data = [sham_pre_A_d_prime_data; d_prime_scores(b, :)];
                elseif face_type == 'N'
                    sham_pre_N_d_prime_data = [sham_pre_N_d_prime_data; d_prime_scores(b, :)];
                end
            end
        elseif pre_or_post == 2
            if session == tms_session
                if face_type == 'A'
                    tms_post_A_d_prime_data = [tms_post_A_d_prime_data; d_prime_scores(b, :)];
                elseif face_type == 'N'
                    tms_post_N_d_prime_data = [tms_post_N_d_prime_data; d_prime_scores(b, :)];
                end
            elseif session == sham_session
                if face_type == 'A'
                    sham_post_A_d_prime_data = [sham_post_A_d_prime_data; d_prime_scores(b, :)];
                elseif face_type =='N'
                    sham_post_N_d_prime_data = [sham_post_N_d_prime_data; d_prime_scores(b, :)];
                end
            end
        end
    end
end

%Create averages for the accuracies of the participants
% Define the condition labels
conditions = {'tms_pre_A_d_prime_data', 'tms_post_A_d_prime_data', 'tms_pre_N_d_prime_data', 'tms_post_N_d_prime_data', ...
    'sham_pre_A_d_prime_data', 'sham_post_A_d_prime_data', 'sham_pre_N_d_prime_data', 'sham_post_N_d_prime_data'};

% Define condition data
condition_data = {tms_pre_A_d_prime_data, tms_post_A_d_prime_data, tms_pre_N_d_prime_data, tms_post_N_d_prime_data, ...
    sham_pre_A_d_prime_data, sham_post_A_d_prime_data, sham_pre_N_d_prime_data, sham_post_N_d_prime_data};

%Save each of the created arrays as an xls file, to do so loop through the
%different conditions as they were saved
for idx = 1:numel(conditions)
    %Retreive the name for the files
    filename = [conditions{idx} '.xls'];

    %write the data table into an xls file
    writecell(condition_data{idx}, filename);
end

% Initialize a structure to store the mean RTs for each participant and condition
mean_d_primes = struct();

% Initialize a cell array to store tables for each condition
condition_tables = cell(1, numel(conditions));

% Iterate through each condition
for condition_idx = 1:numel(conditions)
    % Get the current condition
    condition = conditions{condition_idx};

    % Initialize a structure to store mean RTs for the current condition
    mean_d_primes.(condition) = struct();

    % Iterate through each data element for the current condition
    for data_idx = 1:size(condition_data{1, condition_idx})
        % Extract data for the current condition and data index
        data = condition_data{1, condition_idx}{data_idx, 2};

        % Extract the participant number from the first column for data
        % index
        file_name = condition_data{1, condition_idx}{data_idx, 1};
        file_name_parts = strsplit(file_name, '_');

        % Extract participant number from the first column
        participant_number = str2double(file_name_parts{1}(3:end));

        % Access the tables in the fourth column and calculate mean from RTs in the second column
        mean_d_prime = mean(data);

        % Store the mean d prime for the current participant and condition
        participant_key = ['PP', num2str(participant_number)]; % Prefix with PP
        if isfield(mean_d_primes.(condition), participant_key)
            mean_d_primes.(condition).(participant_key) = ...
                [mean_d_primes.(condition).(participant_key), mean_d_prime];
        else
            mean_d_primes.(condition).(participant_key) = mean_d_prime;
        end
    end
    
    % Get the fieldnames of the mean d prime for the current condition
    participant_numbers = fieldnames(mean_d_primes.(condition));
    
    % Iterate over each participant number
    for c = 1:size(participant_numbers)
        participant_keys = participant_numbers{c};

        % Get the means stored for the participant and condition
        participant_means = mean_d_primes.(condition).(participant_keys);

        % Calculate the mean of the three means
        overall_mean = mean(participant_means);

        % Store the overall mean for the participant and condition
        mean_d_primes.(condition).(participant_keys) = overall_mean;
    end

    % Get structure with means per participant for the current condition
    condition_d_prime_means = mean_d_primes.(condition);

    % Convert the structure to a table
    condition_d_prime_table = struct2table(condition_d_prime_means);

    % Transpose the table to flip rows and columns
    condition_d_prime_table = rows2vars(condition_d_prime_table);

    % Rename first column to participant number
    condition_d_prime_table.Properties.VariableNames{1} = 'Participant_Number';
    
    % Rename second column to the current condition name
    condition_d_prime_table.Properties.VariableNames{2} = condition;

    % Store the table for the current condition
    condition_tables{condition_idx} = condition_d_prime_table;
end

% Combine tables for all conditions into one table
combined_d_prime_table = condition_tables{1};
for condition_idx = 2:numel(conditions)
    combined_d_prime_table = outerjoin(combined_d_prime_table, condition_tables{condition_idx}, 'MergeKeys', true);
end

% Extract numeric part of participant number
combined_d_prime_table.Participant_Number = extractAfter(combined_d_prime_table.Participant_Number, 2);
combined_d_prime_table.Participant_Number = str2double(combined_d_prime_table.Participant_Number);

% Sort the table based on the numeric participant number
combined_d_prime_table = sortrows(combined_d_prime_table, 'Participant_Number');

% Rename the first variable of the TMSorSHAM_data table to match the variable name in combined_table
TMSorSHAM_data.Properties.VariableNames{1} = 'Participant_Number';

% Replace participant numbers with order as found in the TMSorSHAM data
% table created earlier
[~, idx] = ismember(combined_d_prime_table(:, 1), TMSorSHAM_data(:, 1));
combined_d_prime_table(:, 1) = TMSorSHAM_data(idx, 2);

%Read BDI_Scores file
BDI_Scores_file = dir(fullfile(folder_path, "BDI_Scores.xlsx"));
BDI_Scores_file_path = fullfile(folder_path, BDI_Scores_file.name);

%Load BDI_Scores from xls file
BDI_Scores_data = readtable(BDI_Scores_file_path);

% Extract the second column of BDI scores
BDI_scores = BDI_Scores_data{:, 2};

% Initialize BDI_Group variable
BDI_Group = zeros(size(BDI_scores));

% Set values above 10 to 1 (depressed)
BDI_Group(BDI_scores > 13) = 1;

% Set values below 10 to 0 (healthy)
BDI_Group(BDI_scores <= 13) = 0;

% Store existing columns of the combined table
existing_columns = combined_d_prime_table(:, 1:end);

% Add BDI scores as the first column
combined_d_prime_table = table(BDI_scores, 'VariableNames', {'BDI_Scores'});
combined_d_prime_table = table(BDI_Group, 'VariableNames', {'BDI_Group'});

% Concatenate existing columns to the combined table
combined_d_prime_table = [combined_d_prime_table, existing_columns];

% Rename the Participant Number variable name to order
combined_d_prime_table.Properties.VariableNames{2} = 'Order stim';

% Read %RMT file
percent_RMT_file = dir(fullfile(folder_path, "stim_intensity.xlsx"));
percent_RMT_file_path = fullfile(folder_path, percent_RMT_file.name);

% Load the %RMT from file
percent_RMT_data = readtable(percent_RMT_file_path);

% Extract the 3rd column with the %RMT
percent_RMT = percent_RMT_data{:, 3};

% Store existing columns of the combined table
existing_columns = combined_d_prime_table(:, 1:end);

% Add %RMT as the first column
combined_d_prime_table_BDI = table(percent_RMT, 'VariableNames', {'percent_RMT'});

% Concatenate existing columns to the combined table
combined_d_prime_table_BDI = [combined_d_prime_table_BDI, existing_columns];

% Extract the 4th column with the machine output
Machine_output = percent_RMT_data{:, 4};

% Store existing columns of the combined table
existing_columns = combined_d_prime_table(:, 1:end);

% Add Machine output as first column
combined_d_prime_table = table(Machine_output, 'VariableNames', {'Machine_output'});

% Concatenate existing columns to the combined table
combined_d_prime_table = [combined_d_prime_table, existing_columns];

% Save the table as an xls file
filename = 'combined_d_prime_table_BDI.xls'; %Name of the xls file
writetable(combined_d_prime_table_BDI, filename); %Write the table to an xls file


%Initialize arrays to store the means and standard errors
Angry_d_prime = [];
Angry_d_prime_errors = [];
Neutral_d_prime = [];
Neutral_d_prime_errors = [];

%Loop through each of the files
for d = 1:size(filtered_data_tables, 1)
    %Split file name to get the face info of the file
    file_name = filtered_data_tables{d, 1};
    file_name_parts = strsplit(file_name, '_');
    face_type = file_name_parts{5}(1);

    %Extract RT data from filtered_data_tables (IQR data)
    d_prime_data = d_prime_scores{d, 2};

    %Calculate the mean and standard error
    mean_d_prime = mean(d_prime_data);
    std_d_prime_error = std(d_prime_data)/sqrt(length(d_prime_data));

    %Add mean and std error to the arrays    
    if face_type == 'A'
        Angry_d_prime = [Angry_d_prime; mean_d_prime];
        Angry_d_prime_errors = [Angry_d_prime_errors; std_d_prime_error];
    elseif face_type == 'N'
        Neutral_d_prime = [Neutral_d_prime; mean_d_prime];
        Neutral_d_prime_errors = [Neutral_d_prime_errors; std_d_prime_error];
    end
end

% Calculate overall mean and standard error for Angry and Neutral faces
overall_Angry_d_prime = mean(Angry_d_prime);
overall_Neutral_d_prime = mean(Neutral_d_prime);
overall_Angry_std_error_d_prime = std(Angry_d_prime) / sqrt(length(Angry_d_prime));
overall_Neutral_std_error_d_prime = std(Neutral_d_prime) / sqrt(length(Neutral_d_prime));

% Create a scatter plot of the means with error bars
figure;
hold on;
scatter(size(overall_Angry_d_prime), overall_Angry_d_prime, 'b');
scatter(2*size(overall_Neutral_d_prime), overall_Neutral_d_prime, 'r')
errorbar([1, 2], [overall_Angry_d_prime, overall_Neutral_d_prime], [overall_Angry_std_error_d_prime, overall_Neutral_std_error_d_prime], 'k', 'LineWidth', 1.5, 'LineStyle', 'none');
plot([1, 2], [overall_Angry_d_prime, overall_Neutral_d_prime], 'k--', 'LineWidth', 1);
xlim([0.5, 2.5]);
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Angry', 'Neutral'});
xlabel('Face Type');
ylabel('Mean d prime');
title('Overall Mean d prime for Angry and Neutral Faces');
hold off
