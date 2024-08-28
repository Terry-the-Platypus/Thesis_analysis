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

%Code for cleaning up the RT data for all PP and find the mean RTs per file
%name

%Initialize a cell array to store tables and thier coresponding file names,
%as well as removed data
filtered_data_tables = cell(length(files), 4); % 4 columns {file name, filtered data}
removed_data_tables = cell(length(files), 4); % 4 columns {file name, removed data}

%Loop through each table stored in second column of data_table
for a = 1:size(data_tables, 1)
    %Retreive the table from second column
    table_data = data_tables{a, 2};

    %Filter the data table such that only trials with a numeric value for
    %RT and a target trial
    filtered_data = table_data(~isnan(table_data.RT) | table_data.target == 1, :);

    %Create variable containing the non-targets as well as the targets for
    %which there was no rt value
    non_targets = table_data(table_data.target == 0, :);
    no_rt_targets = table_data(isnan(table_data.RT) & table_data.target == 1, :);

    %Remove all RT values <200ms and >1000ms
    filtered_data_RT = filtered_data(filtered_data.RT >= 200 & filtered_data.RT <= 1000, :);
    
    %Create variable for responses that were too fast
    RT_data_too_fast = filtered_data(filtered_data.RT <= 200, :);

    %Create variable for responses that were too slow
    RT_data_too_slow = filtered_data(filtered_data.RT >= 1000, :);

    %Calc the IQR for the RT column
    rt_column = filtered_data_RT.RT;
    iqr_rt = iqr(rt_column);

    %Compute lower and upper bounds based on IQR
    lower_bound = quantile(rt_column, 0.25) - 1.5 * iqr_rt;
    upper_bound = quantile(rt_column, 0.75) + 1.5 * iqr_rt;

    %Filter data to inlcude only RTs w/in IQR
    filtered_data_IQR = filtered_data_RT(rt_column >= lower_bound & rt_column <= upper_bound, :);

    %Create table containing the lower and upper bounds of the RT data
    lower_IQR_data = filtered_data_RT(rt_column < lower_bound, :);
    upper_IQR_data = filtered_data_RT(rt_column > upper_bound, :);
    lower_and_upper_IQR_data = vertcat(lower_IQR_data, upper_IQR_data);

    %Store filtered data table along with file name
    filtered_data_tables{a, 1} = files(a).name;
    filtered_data_tables{a, 2} = filtered_data;
    filtered_data_tables{a, 3} = filtered_data_RT;
    filtered_data_tables{a, 4} = filtered_data_IQR;

    %Save the data to be used for further analysis

    %Store removed data tables along with file name
    removed_data_tables{a, 1} = files(a).name;
    removed_data_tables{a, 2} = non_targets;
    removed_data_tables{a, 3} = no_rt_targets;
    removed_data_tables{a, 4} = RT_data_too_slow;
    removed_data_tables{a, 5} = RT_data_too_fast;
    removed_data_tables{a, 6} = lower_and_upper_IQR_data;
end




%In order to make a histogram of the removed trials per participant
%nummber, first create a counts array
removal_counts = zeros(60, 4);

%Loop through each file name in removed_data_tables
for b = 1:size(removed_data_tables, 1)
    %Extract info from current file name
    file_name = removed_data_tables{b, 1};
    file_name_parts = strsplit(file_name, '_');
    PP_number = str2double(file_name_parts{1}(3:end));

    %Iterate through each removal type (columns 3 to 6)
    for removal_type = 3:6
        %Count the number of trials removed for current removal type
        removal_counts(PP_number, removal_type-2) = size(removed_data_tables{PP_number, removal_type}, 1);
    end
end

% Plot histogram
figure;
bar(removal_counts);
xlabel('Participant Number');
ylabel('Number of Trials Removed');
title('Histogram of Trials Removed by Participant Number and Removal Type');
legend({'No RT Targets', 'Too Slow RTs', 'Too Fast RTs', 'Data Outside IQR'}, 'Location', 'best');

%In order to make a histogram of the included trials per participant
%nummber, first create a counts array
inclusion_counts = zeros(60, 4);

%Loop through each file name in filtered_data_tables
for b = 1:size(filtered_data_tables, 1)
    %Extract info from current file name
    file_name = filtered_data_tables{b, 1};
    file_name_parts = strsplit(file_name, '_');
    PP_number = str2double(file_name_parts{1}(3:end));

    %Iterate through each filter type (columns 2 to 4)
    for filter_type =2:4
        %Count the number of trials included for current filter type
        inclusion_counts(PP_number, filter_type-1) = size(filtered_data_tables{PP_number, filter_type}, 1);
    end
end

% Plot histogram
figure;
bar(inclusion_counts);
xlabel('Participant Number');
ylabel('Number of Trials Included');
title('Histogram of Trials Included by Participant Number and Filter Type');
legend({'Target or RT value', 'RT filter (fast and slow)', 'Data Within IQR'}, 'Location', 'best');


%Sort files as TMS or SHAM stimulation
%Read TMSorSHAM file
TMSorSHAM_file = dir(fullfile(folder_path, "TMSorSHAM.xlsx"));
TMSorSHAM_file_path = fullfile(folder_path, TMSorSHAM_file.name);

%Load order from TMSorSHAM xls file
TMSorSHAM_data = readtable(TMSorSHAM_file_path);

% Initialize cell arrays to store TMS and SHAM data
tms_pre_A_data = {};
tms_post_A_data = {};
tms_pre_N_data = {};
tms_post_N_data = {};
sham_pre_A_data = {};
sham_post_A_data = {};
sham_pre_N_data = {};
sham_post_N_data = {};

% Loop through each file name in filtered_data_tables
for b = 1:size(filtered_data_tables, 1)
    % Extract information from the current file name
    file_name = filtered_data_tables{b, 1};
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
                    tms_pre_A_data = [tms_pre_A_data; filtered_data_tables(b, :)];
                elseif face_type == 'N'
                    tms_pre_N_data = [tms_pre_N_data; filtered_data_tables(b, :)];
                end
            elseif session == sham_session
                if face_type =='A'
                    sham_pre_A_data = [sham_pre_A_data; filtered_data_tables(b, :)];
                elseif face_type == 'N'
                    sham_pre_N_data = [sham_pre_N_data; filtered_data_tables(b, :)];
                end
            end
        elseif pre_or_post == 2
            if session == tms_session
                if face_type == 'A'
                    tms_post_A_data = [tms_post_A_data; filtered_data_tables(b, :)];
                elseif face_type == 'N'
                    tms_post_N_data = [tms_post_N_data; filtered_data_tables(b, :)];
                end
            elseif session == sham_session
                if face_type == 'A'
                    sham_post_A_data = [sham_post_A_data; filtered_data_tables(b, :)];
                elseif face_type =='N'
                    sham_post_N_data = [sham_post_N_data; filtered_data_tables(b, :)];
                end
            end
        end
    end
end


% Define the condition labels
conditions = {'tms_pre_A_data', 'tms_post_A_data', 'tms_pre_N_data', 'tms_post_N_data', ...
              'sham_pre_A_data', 'sham_post_A_data', 'sham_pre_N_data', 'sham_post_N_data'};

% Define condition data
condition_data = {tms_pre_A_data, tms_post_A_data, tms_pre_N_data, tms_post_N_data, ...
                  sham_pre_A_data, sham_post_A_data, sham_pre_N_data, sham_post_N_data};

% Initialize a structure to store the mean RTs for each participant and condition
mean_RTs = struct();
raw_RT = struct();

% Initialize a cell array to store tables for each condition
condition_tables = cell(1, numel(conditions));

% Iterate through each condition
for condition_idx = 1:numel(conditions)
    % Get the current condition
    condition = conditions{condition_idx};
    
    % Initialize a structure to store mean RTs for the current condition
    mean_RTs.(condition) = struct();
    raw_RT.(condition) = struct();
    
    % Iterate through each data element for the current condition
    for data_idx = 1:size(condition_data{1, condition_idx})
        % Extract data for the current condition and data index
        data = condition_data{1, condition_idx}{data_idx, 4};
        % Extract the participant number from the first column for data
        % index
        file_name = condition_data{1, condition_idx}{data_idx, 1};
        file_name_parts = strsplit(file_name, '_');

        % Extract participant number from the first column
        participant_number = str2double(file_name_parts{1}(3:end));
        
        % Access the tables in the fourth column and calculate mean from RTs in the second column
        RTs = data.RT;

        % Access the tables in the fourth column and extract RTs from the second column
        RTs = data.RT;
        
        % Store raw RTs for the current participant and condition
        participant_key = ['PP', num2str(participant_number)]; % Prefix with PP
        if isfield(raw_RT.(condition), participant_key)
            raw_RT.(condition).(participant_key) = ...
                [raw_RT.(condition).(participant_key), {RTs}];
        else
            raw_RT.(condition).(participant_key) = {RTs};
        end


        mean_RT = mean(RTs);
        
        % Store the mean RT for the current participant and condition
        participant_key = ['PP', num2str(participant_number)];% Prefix with PP
        if isfield(mean_RTs.(condition), num2str(participant_key))
            mean_RTs.(condition).(num2str(participant_key)) = ...
                [mean_RTs.(condition).(num2str(participant_key)), mean_RT];
        else
            mean_RTs.(condition).(num2str(participant_key)) = mean_RT;
        end
    end

    % Get the fieldnames of the mean RTs for the current condition
    participant_numbers = fieldnames(mean_RTs.(condition));

    % Iterate over each participant number
    for c = 1:size(participant_numbers)
        participant_keys = participant_numbers{c};
        
        % Get the means stored for the participant and condition
        participant_means = mean_RTs.(condition).(participant_keys);
            
        % Calculate the mean of the three means
        overall_mean = mean(participant_means);
            
        % Store the overall mean for the participant and condition
        mean_RTs.(condition).(participant_keys) = overall_mean;
    end

    % Get structure with means per participant for the current condition
    condition_means = mean_RTs.(condition);

    % Convert the structure to a table
    condition_table = struct2table(condition_means);

    % Transpose the table to flip rows and columns
    condition_table = rows2vars(condition_table);

    % Rename first column to participant number
    condition_table.Properties.VariableNames{1} = 'Participant_Number';

    % Rename second column to the current condition name
    condition_table.Properties.VariableNames{2} = condition;

    % Store the table for the current condition
    condition_tables{condition_idx} = condition_table;
end

% Save raw_RT structure for later use in RStudio
save('raw_RT.mat', 'raw_RT');

% Combine tables for all conditions into one table
combined_table = condition_tables{1};
for condition_idx = 2:numel(conditions)
    combined_table = outerjoin(combined_table, condition_tables{condition_idx}, 'MergeKeys', true);
end

% Extract numeric part of participant number
combined_table.Participant_Number = extractAfter(combined_table.Participant_Number, 2);
combined_table.Participant_Number = str2double(combined_table.Participant_Number);

% Sort the table based on the numeric participant number
combined_table = sortrows(combined_table, 'Participant_Number');

% Rename the first variable of the TMSorSHAM_data table to match the variable name in combined_table
TMSorSHAM_data.Properties.VariableNames{1} = 'Participant_Number';

% Replace participant numbers with order as found in the TMSorSHAM data
% table created earlier
[~, idx] = ismember(combined_table(:, 1), TMSorSHAM_data(:, 1));
combined_table(:, 1) = TMSorSHAM_data(idx, 2);

% Read BDI_Scores file
BDI_Scores_file = dir(fullfile(folder_path, "BDI_Scores.xlsx"));
BDI_Scores_file_path = fullfile(folder_path, BDI_Scores_file.name);

% Load BDI_Scores from xls file
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
existing_columns = combined_table(:, 1:end);

% Add BDI scores as the first column
combined_table = table(BDI_Group, 'VariableNames', {'BDI_Group'});

% Concatenate existing columns to the combined table
combined_table = [combined_table, existing_columns];

% Rename the Participant Number variable name to order
combined_table.Properties.VariableNames{2} = 'Order_stim';

% Read %RMT file
percent_RMT_file = dir(fullfile(folder_path, "stim_intensity.xlsx"));
percent_RMT_file_path = fullfile(folder_path, percent_RMT_file.name);

% Load the %RMT from file
percent_RMT_data = readtable(percent_RMT_file_path);

% Extract the 3rd column with the %RMT
percent_RMT = percent_RMT_data{:, 3};

% Store existing columns of the combined table
existing_columns = combined_table(:, 1:end);

% Add %RMT as the first column
combined_table_BDI = table(percent_RMT, 'VariableNames', {'percent_RMT'});

% Concatenate existing columns to the combined table
combined_table_BDI = [combined_table_BDI, existing_columns];

% Extract the 4th column with the machine output
Machine_output = percent_RMT_data{:, 4};

% Store existing columns of the combined table
existing_columns = combined_table(:, 1:end);

% Add Machine output as first column
combined_table = table(Machine_output, 'VariableNames', {'Machine_output'});

% Concatenate existing columns to the combined table
combined_table = [combined_table, existing_columns];

% Save the table as an xls file
filename = 'combined_RT_table_BDI.xls'; %Name of the xls file
writetable(combined_table_BDI, filename); %Write the table to an xls file

%Initialize arrays to store the means and standard errors
Angry_means = [];
Angry_errors = [];
Neutral_means = [];
Neutral_errors = [];

%Loop through each of the files
for d = 1:size(filtered_data_tables, 1)
    %Split file name to get the face info of the file
    file_name = filtered_data_tables{d, 1};
    file_name_parts = strsplit(file_name, '_');
    face_type = file_name_parts{5}(1);

    %Extract RT data from filtered_data_tables (IQR data)
    RT_data = filtered_data_tables{d, 4}.RT;

    %Calculate the mean and standard error
    mean_RT = mean(RT_data);
    std_error = std(RT_data)/sqrt(length(RT_data));

    %Add mean and std error to the arrays    
    if face_type == 'A'
        Angry_means = [Angry_means; mean_RT];
        Angry_errors = [Angry_errors; std_error];
    elseif face_type == 'N'
        Neutral_means = [Neutral_means; mean_RT];
        Neutral_errors = [Neutral_errors; std_error];
    end
end

% Calculate overall mean and standard error for Angry and Neutral faces
overall_Angry_mean = mean(Angry_means);
overall_Neutral_mean = mean(Neutral_means);
overall_Angry_std_error = std(Angry_means) / sqrt(length(Angry_means));
overall_Neutral_std_error = std(Neutral_means) / sqrt(length(Neutral_means));

% Create a scatter plot of the means with error bars
figure;
hold on;
scatter(size(overall_Angry_mean), overall_Angry_mean, 'b');
scatter(2*size(overall_Neutral_mean), overall_Neutral_mean, 'r')
errorbar([1, 2], [overall_Angry_mean, overall_Neutral_mean], [overall_Angry_std_error, overall_Neutral_std_error], 'k', 'LineWidth', 1.5, 'LineStyle', 'none');
plot([1, 2], [overall_Angry_mean, overall_Neutral_mean], 'k--', 'LineWidth', 1);
xlim([0.5, 2.5]);
set(gca, 'XTick', [1, 2], 'XTickLabel', {'Angry', 'Neutral'});
xlabel('Face Type');
ylabel('Mean RT');
title('Overall Mean Reaction Time for Angry and Neutral Faces');
hold off

