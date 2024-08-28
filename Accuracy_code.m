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

    %Find the total number of targets for each file
    total_targets = table_data(table_data.target == 1, :);

    %Find the correct rejections for each file
    correct_rejections = table_data(isnan(table_data.RT) & table_data.target == 0, :);

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
    
    %Filter the data table such that we get all the correct button presses
    %made
    correct_presses = IQR_filtered_data(~isnan(IQR_filtered_data.RT) & IQR_filtered_data.target == 1, :);

    %Find the false alarms (button pressed with no target) for each file
    false_alarms = IQR_filtered_data(~isnan(IQR_filtered_data.RT) & IQR_filtered_data.target == 0, :);

    %Find teh number of missed targets for each file
    missed_targets = table_data(isnan(table_data.RT) & table_data.target == 1, :);

    %Create table containing the lower and upper bounds of the RT data
    Q1_IQR_data = time_filtered_data(rt_column < lower_bound, :);
    Q4_IQR_data = time_filtered_data(rt_column > upper_bound, :);
    Q1_and_Q4_IQR = vertcat(Q1_IQR_data, Q4_IQR_data);

    %Store filtered data table along with file name
    filtered_data_tables{a, 1} = files(a).name;
    filtered_data_tables{a, 2} = correct_rejections;
    filtered_data_tables{a, 3} = time_filtered_data;
    filtered_data_tables{a, 4} = IQR_filtered_data;
    filtered_data_tables{a, 5} = correct_presses;

    %Store removed data tables along with file name
    removed_data_tables{a, 1} = files(a).name;
    removed_data_tables{a, 2} = false_alarms;
    removed_data_tables{a, 3} = missed_targets;
    removed_data_tables{a, 4} = RT_too_slow;
    removed_data_tables{a, 5} = RT_too_fast;
    removed_data_tables{a, 6} = Q1_and_Q4_IQR;
end

%Make histogram of filtered_data_tables to see the dfference between the
%number of filtered correct and correct presses
%Create counts array with number of paticipants
inclusion_counts = zeros(60, 3);

%Loop through each file name in filtered_data_tables
for b = 1:size(filtered_data_tables, 1)
    %Extract info from current file name
    file_name = filtered_data_tables{b, 1};
    file_name_parts = strsplit(file_name, '_');
    PP_number = str2double(file_name_parts{1}(3:end));

    %Iterate through each filter type (columns 2 to 5)
    for filter_type = 2:5
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
legend({'Correct Rejections', 'Time Filtered Data', 'IQR Filtered Data', 'Filtered Correct Presses'}, 'Location', 'best');


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
    for removal_type = 2:6
        %Count the number of trials removed for current removal type
        removal_counts(PP_number, removal_type-1) = size(removed_data_tables{PP_number, removal_type}, 1);
    end
end

% Plot histogram
figure;
bar(removal_counts);
xlabel('Participant Number');
ylabel('Number of Trials Removed');
title('Histogram of Trials Removed by Participant Number and Removal Type');
legend({'False alarms', 'Missed Targets', 'Too Slow RTs', 'Too Fast RTs', 'Data Outside IQR'}, 'Location', 'best');



%Calculate the accuracy per file based on the data in filtered_data_tables,
%looking at filtered correct presses and total targets
%Create a cell array to store the accuracies per file in
accuracy_scores = cell(length(filtered_data_tables), 2);

%Loop through each file stored in filtered_data_tables
for c = 1:size(filtered_data_tables, 1)
    %Retreive the number of filtered correct button presses
    hits = size(filtered_data_tables{c, 5});

    %Retreive the number of correct rejections
    corr_rejections = size(filtered_data_tables{c, 2});

    %Retreive the number of missed targets
    misses = size(removed_data_tables{c, 3});

    %Retreive the number of false alarms
    FAs = size(removed_data_tables{c, 2});

    %Calculate the accuracy per file based on Signal Detection Theory for
    %accuracy calculation (hits + correct rejections)/total trials (correct
    %rejections, false alarms, hits and missees)
    numb_hits_corr_rej = hits(1) + corr_rejections(1);
    total_trials = hits(1) + corr_rejections(1) + misses(1) + FAs(1);
    accuracy_score = numb_hits_corr_rej/total_trials; 

    %Store the accuracy score in the cell array with the file name
    accuracy_scores{c, 1} = filtered_data_tables{c, 1};
    accuracy_scores{c, 2} = accuracy_score;
end

%Save the accuracy scores as an xls file
filename = 'Accuracy_scores.xls';
writecell(accuracy_scores, filename);



%Split the data based on the different conditions

%Sort files as TMS or SHAM stimulation
%Read TMSorSHAM file
TMSorSHAM_file = dir(fullfile(folder_path, "TMSorSHAM.xlsx"));
TMSorSHAM_file_path = fullfile(folder_path, TMSorSHAM_file.name);

%Load order from TMSorSHAM xls file
TMSorSHAM_data = readtable(TMSorSHAM_file_path);

% Initialize cell arrays to store TMS and SHAM data
tms_pre_A_accuracy_data = {};
tms_post_A_accuracy_data = {};
tms_pre_N_accuracy_data = {};
tms_post_N_accuracy_data = {};
sham_pre_A_accuracy_data = {};
sham_post_A_accuracy_data = {};
sham_pre_N_accuracy_data = {};
sham_post_N_accuracy_data = {};

% Loop through each file name in filtered_data_tables
for b = 1:size(accuracy_scores, 1)
    % Extract information from the current file name
    file_name = accuracy_scores{b, 1};
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
                    tms_pre_A_accuracy_data = [tms_pre_A_accuracy_data; accuracy_scores(b, :)];
                elseif face_type == 'N'
                    tms_pre_N_accuracy_data = [tms_pre_N_accuracy_data; accuracy_scores(b, :)];
                end
            elseif session == sham_session
                if face_type =='A'
                    sham_pre_A_accuracy_data = [sham_pre_A_accuracy_data; accuracy_scores(b, :)];
                elseif face_type == 'N'
                    sham_pre_N_accuracy_data = [sham_pre_N_accuracy_data; accuracy_scores(b, :)];
                end
            end
        elseif pre_or_post == 2
            if session == tms_session
                if face_type == 'A'
                    tms_post_A_accuracy_data = [tms_post_A_accuracy_data; accuracy_scores(b, :)];
                elseif face_type == 'N'
                    tms_post_N_accuracy_data = [tms_post_N_accuracy_data; accuracy_scores(b, :)];
                end
            elseif session == sham_session
                if face_type == 'A'
                    sham_post_A_accuracy_data = [sham_post_A_accuracy_data; accuracy_scores(b, :)];
                elseif face_type =='N'
                    sham_post_N_accuracy_data = [sham_post_N_accuracy_data; accuracy_scores(b, :)];
                end
            end
        end
    end
end

%Create averages for the accuracies of the participants
% Define the condition labels
conditions = {'tms_pre_A_accuracy_data', 'tms_post_A_accuracy_data', 'tms_pre_N_accuracy_data', 'tms_post_N_accuracy_data', ...
              'sham_pre_A_accuracy_data', 'sham_post_A_accuracy_data', 'sham_pre_N_accuracy_data', 'sham_post_N_accuracy_data'};

% Define condition data
condition_data = {tms_pre_A_accuracy_data, tms_post_A_accuracy_data, tms_pre_N_accuracy_data, tms_post_N_accuracy_data, ...
                  sham_pre_A_accuracy_data, sham_post_A_accuracy_data, sham_pre_N_accuracy_data, sham_post_N_accuracy_data};

%Save each of the created arrays as an xls file, to do so loop through the
%different conditions as they were saved
for idx = 1:numel(conditions)
    %Retreive the name for the files
    filename = [conditions{idx} '.xls'];

    %write the data table into an xls file
    writecell(condition_data{idx}, filename);
end


% Initialize a structure to store the mean RTs for each participant and condition
mean_accuracies = struct();

% Initialize a cell array to store tables for each condition
condition_tables = cell(1, numel(conditions));

% Iterate through each condition
for condition_idx = 1:numel(conditions)
    % Get the current condition
    condition = conditions{condition_idx};
    
    % Initialize a structure to store mean RTs for the current condition
    mean_accuracies.(condition) = struct();
    
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
        mean_accuracy = mean(data);

        % Store the mean RT for the current participant and condition
        participant_key = ['PP', num2str(participant_number)];% Prefix with PP
        if isfield(mean_accuracies.(condition), num2str(participant_key))
            mean_accuracies.(condition).(num2str(participant_key)) = ...
                [mean_accuracies.(condition).(num2str(participant_key)), mean_accuracy];
        else
            mean_accuracies.(condition).(num2str(participant_key)) = mean_accuracy;
        end
    end

    % Get the fieldnames of the mean RTs for the current condition
    participant_numbers = fieldnames(mean_accuracies.(condition));

    % Iterate over each participant number
    for c = 1:size(participant_numbers)
        participant_keys = participant_numbers{c};
        
        % Get the means stored for the participant and condition
        participant_means = mean_accuracies.(condition).(participant_keys);
            
        % Calculate the mean of the three means
        overall_mean = mean(participant_means);
            
        % Store the overall mean for the participant and condition
        mean_accuracies.(condition).(participant_keys) = overall_mean;
    end

    % Get structure with means per participant for the current condition
    condition_accuracy_means = mean_accuracies.(condition);

    % Convert the structure to a table
    condition_accuracy_table = struct2table(condition_accuracy_means);

    % Transpose the table to flip rows and columns
    condition_accuracy_table = rows2vars(condition_accuracy_table);

    % Rename first column to participant number
    condition_accuracy_table.Properties.VariableNames{1} = 'Participant_Number';

    % Rename second column to the current condition name
    condition_accuracy_table.Properties.VariableNames{2} = condition;

    % Store the table for the current condition
    condition_tables{condition_idx} = condition_accuracy_table;
end

% Combine tables for all conditions into one table
combined_accuracy_table = condition_tables{1};
for condition_idx = 2:numel(conditions)
    combined_accuracy_table = outerjoin(combined_accuracy_table, condition_tables{condition_idx}, 'MergeKeys', true);
end

% Extract numeric part of participant number
combined_accuracy_table.Participant_Number = extractAfter(combined_accuracy_table.Participant_Number, 2);
combined_accuracy_table.Participant_Number = str2double(combined_accuracy_table.Participant_Number);

% Sort the table based on the numeric participant number
combined_accuracy_table = sortrows(combined_accuracy_table, 'Participant_Number');

% Rename the first variable of the TMSorSHAM_data table to match the variable name in combined_table
TMSorSHAM_data.Properties.VariableNames{1} = 'Participant_Number';

% Replace participant numbers with order as found in the TMSorSHAM data
% table created earlier
[~, idx] = ismember(combined_accuracy_table(:, 1), TMSorSHAM_data(:, 1));
combined_accuracy_table(:, 1) = TMSorSHAM_data(idx, 2);

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
existing_columns = combined_accuracy_table(:, 1:end);

% Add BDI scores as the first column
combined_accuracy_table = table(BDI_scores, 'VariableNames', {'BDI_Scores'});
combined_accuracy_table = table(BDI_Group, 'VariableNames', {'BDI_Group'});

% Concatenate existing columns to the combined table
combined_accuracy_table = [combined_accuracy_table, existing_columns];

% Rename the Participant Number variable name to order
combined_accuracy_table.Properties.VariableNames{2} = 'Order stim';

% Read %RMT file
percent_RMT_file = dir(fullfile(folder_path, "stim_intensity.xlsx"));
percent_RMT_file_path = fullfile(folder_path, percent_RMT_file.name);

% Load the %RMT from file
percent_RMT_data = readtable(percent_RMT_file_path);

% Extract the 3rd column with the %RMT
percent_RMT = percent_RMT_data{:, 3};

% Store existing columns of the combined table
existing_columns = combined_accuracy_table(:, 1:end);

% Add %RMT as the first column
combined_accuracy_table_BDI = table(percent_RMT, 'VariableNames', {'percent_RMT'});

% Concatenate existing columns to the combined table
combined_accuracy_table_BDI = [combined_accuracy_table_BDI, existing_columns];

% Extract the 4th column with the machine output
Machine_output = percent_RMT_data{:, 4};

% Store existing columns of the combined table
existing_columns = combined_accuracy_table(:, 1:end);

% Add Machine output as first column
combined_accuracy_table = table(Machine_output, 'VariableNames', {'Machine_output'});

% Concatenate existing columns to the combined table
combined_accuracy_table = [combined_accuracy_table, existing_columns];

% Save the table as an xls file
filename = 'combined_accuracy_table_BDI.csv'; %Name of the xls file
writetable(combined_accuracy_table_BDI, filename); %Write the table to an xls file

