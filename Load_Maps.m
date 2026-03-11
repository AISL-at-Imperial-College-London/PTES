global gpr_comp_map_PR gpr_comp_map_eta gpr_comp_map_flow
global speed_mean speed_std flow_mean flow_std PR_mean PR_std eta_mean eta_std

file_path  = 'Z_map.xlsx';
sheet_name = 'Map';

% Check if trained models already exist
model_file = 'trained_gp_models.mat';

if isfile(model_file)
    % If saved models exist, load them instead of refitting
    load(model_file, 'gpr_comp_map_flow', 'gpr_comp_map_PR', 'gpr_comp_map_eta', ...
         'speed_mean', 'speed_std', 'flow_mean', 'flow_std', ...
         'PR_mean', 'PR_std', 'eta_mean', 'eta_std');
    disp('Loaded trained GP models from file.');
else
    % Read data
    compressor_data = xlsread(file_path, sheet_name);

    % Compute scaling factors
    speed_mean = mean(compressor_data(:,1));
    speed_std  = std(compressor_data(:,1));
    Scaled_speed = (compressor_data(:,1) - speed_mean) ./ speed_std;

    flow_mean = mean(compressor_data(:,2)); 
    flow_std  = std(compressor_data(:,2));
    Scaled_flow = (compressor_data(:,2) - flow_mean) ./ flow_std;

    PR_mean = mean(compressor_data(:,3));
    PR_std  = std(compressor_data(:,3));
    Scaled_PR = (compressor_data(:,3) - PR_mean) ./ PR_std;

    eta_mean = mean(compressor_data(:,4));
    eta_std  = std(compressor_data(:,4));
    Scaled_eta = (compressor_data(:,4) - eta_mean) ./ eta_std;

    Scaled_comp_data = [Scaled_speed, Scaled_flow, Scaled_PR, Scaled_eta];

    % Train GP models
    gpr_comp_map_flow = fitrgp(Scaled_comp_data(:,1), Scaled_comp_data(:,2), 'KernelFunction', 'exponential');
    gpr_comp_map_PR   = fitrgp([Scaled_comp_data(:,1), Scaled_comp_data(:,2)], Scaled_comp_data(:,3), 'KernelFunction', 'exponential');
    gpr_comp_map_eta  = fitrgp([Scaled_comp_data(:,1), Scaled_comp_data(:,2)], Scaled_comp_data(:,4), 'KernelFunction', 'exponential');

    % Save the trained models and scaling factors
    save(model_file, 'gpr_comp_map_flow', 'gpr_comp_map_PR', 'gpr_comp_map_eta', ...
         'speed_mean', 'speed_std', 'flow_mean', 'flow_std', ...
         'PR_mean', 'PR_std', 'eta_mean', 'eta_std');
    disp('Trained and saved GP models.');
end