%% Contour map with efficiency island
figure(1);
% Define flows and speeds ranges
flows  = (0:0.1:70)';
speeds = (1:10:380)';

% Initialize matrices for PR and efficiency
p_ratio_ss = zeros(length(flows), length(speeds));
eff        = zeros(length(flows), length(speeds));
hold on

% Loop over the range of speeds and flows to predict PR and efficiency
for j = 1:length(speeds)
    % out.omega_comp = speeds(j);

    for i = 1:length(flows)
        omega_scaled = (speeds(j) / 2 / pi * 60 - speed_mean) / speed_std; % Scale the speed
        flow_scaled  = (flows(i) - flow_mean) / flow_std;                    % Scale the flow

        % Predict PR and efficiency using GPR models
        p_ratio_ss(i, j) = predict(gpr_comp_map_PR, [omega_scaled, flow_scaled]) * PR_std + PR_mean * 1.3;
        eff(i, j)        = predict(gpr_comp_map_eta, [omega_scaled, flow_scaled]) * eta_std + eta_mean;
        % eff(i, j)        = max(eff(i, j), 0.15);
        % eff(i, j)        = min(eff(i, j), 0.9);
    end

    % Find and plot the surge line
    [y_max, z] = max(p_ratio_ss(:, j));
    surge(j, :) = [y_max, flows(z)];
    plot(flows, p_ratio_ss(:, j), '-.','Color', [0 0 0]);
end
s=[];
% Meshgrid for flow and corresponding pressure ratios
[x, y] = meshgrid(flows, linspace(min(p_ratio_ss(:)), max(p_ratio_ss(:)), length(speeds)));

% Plot the efficiency islands using contour
contour(x, y, eff',15,'ShowText', 'off','LineWidth', 1); % Efficiency island with labels
contour(x, y, eff',5,'ShowText', 'on','LineWidth', 1);
% Plot other components (surge line, mass flow vs PR curve)
axis([0 70 1 3.3])
s(:,1)=[1.89536 1.96873 2.04146 2.11012 2.19768 2.28181]*1.3;
s(:,2)=[22.2062 26.6923 31.2905 36.2252 40.3749 44.861];
plot(s(:, 2), s(:, 1), 'r', 'LineWidth', 2);   % Surge line
plot(out_comp.mass, out_comp.P_ratio, 'r--', 'LineWidth', 2);         % Mass flow vs PR curve

plot(out_comp.mass(end), out_comp.P_ratio(end),'bo','LineWidth', 2);
plot(out_comp.mass(1), out_comp.P_ratio(1),'bo','LineWidth', 2);
text(out_comp.mass(end), out_comp.P_ratio(end), 'finish');
text(out_comp.mass(1), out_comp.P_ratio(1), 'start');

title('Compressor map and operating routine');
xlabel('Mass flow (kg/s)')
ylabel('Pressure ratio')
grid on;
hold off;

