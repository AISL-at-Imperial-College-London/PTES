try 
    gpr_comp_map_eta;
catch
    Load_Maps;
end
% Load the final state if it exists
state_file = 'x_ss_v3_speed.mat';
modelName = 'C_PTES_sim_v2';
try
    steady_state=x_ss_v3_speed;
catch
    if isfile(state_file)
        load(state_file, 'x_ss_v3_speed');  % Load the steady-state values
        % disp('Loaded steady_state for initialization.');
        steady_state = x_ss_v3_speed;
    else
        error('State file (x_ss_v3.mat) not found! Run a steady-state simulation first.');
    end
end

load('splines.mat');
load('external.mat');

Param.arch = computer;
Param.splines= splines;
Param.water = s.Water;
Param.air = s.Air;
try
    Param.calc_HTC = calc_HTC;
catch
    Param.calc_HTC = 0;
end
Param.REC.N = 60;
Param.HHX.N = 50;
Param.CHX.N = 30;
Param.COL.N = 6;
hxNames = {'REC','HHX','CHX','COL'};

% hxNames = {'REC','HHX','COL'};

for i=1:numel(hxNames)
    HX= hxNames{i};
    newSize = Param.(HX).N;
    newvals=[];
    for k = 1:numel(steady_state.signals)
        blockName=steady_state.signals(k).blockName;
        if contains(blockName,HX)
            dim=length(steady_state.signals(k).values);
            if dim == newSize
                disp('no change')
            else
                disp('new size wrong')
                newvals=[repelem(mean(steady_state.signals(k).values(1:dim/2)),newSize/2),...
                    repelem(mean(steady_state.signals(k).values(dim/2+1:dim)),newSize/2)];
                steady_state.signals(k).values=newvals;
            end
        end
    end
end





sfuncBlocks = find_system(modelName);
for k = 1:numel(sfuncBlocks)
    if contains(sfuncBlocks{k},'_Func','IgnoreCase',true)
        set_param(sfuncBlocks{k},'Parameters','Param')
    end
end

set_param(modelName, 'LoadInitialState','on','InitialState','steady_state')



% set_param(modelName,...
%     'SaveOutput','on','OutputSaveName','yout',...
%     'SaveState','on','StateSaveName','xout',...
%     'SaveTime','on','TimeSaveName','tout',...
%     'SaveFormat','StructureWithTime');
set_param(modelName,...
    'SaveOutput','on','OutputSaveName','yout',...
    'StateSaveName','xout',...
    'SaveTime','on','TimeSaveName','tout',...
    'SaveFormat','StructureWithTime');


set_param('C_PTES_sim_v2','AlgebraicLoopMsg','none');

try
    Kp_tau;
    Ki_tau;
    stepping_time;
    DeltaW;
catch
    Kp_tau=289.6844;
    Ki_tau=23.1297;
    tchange = 100000;
    % stepping_time=input('Step time:');
    stepping_time = 100000;
    % DeltaW=input('Power change: ');
    DeltaW=0;
    % omega_step=input('Omega_step (376 nominal): ');
    omega_step = 376;
end