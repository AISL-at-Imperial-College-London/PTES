%% Contour map with efficiency island
figure(1);
% subplot(1,2,1);
% Define flows and speeds ranges
flows  = (0:0.1:70)';
speeds = (200:10:380)';

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
surge_line=[];
% Meshgrid for flow and corresponding pressure ratios
[x, y] = meshgrid(flows, linspace(min(p_ratio_ss(:)), max(p_ratio_ss(:)), length(speeds)));

% Plot the efficiency islands using contour
% contour(x, y, eff',15,'ShowText', 'off','LineWidth', 1); % Efficiency island with labels
% contour(x, y, eff',5,'ShowText', 'on','LineWidth', 1);
% Plot other components (surge line, mass flow vs PR curve)
axis([0 70 1 3.3])
surge_line(:,1)=[1.89536 1.96873 2.04146 2.11012 2.19768 2.28181]*1.3;
surge_line(:,2)=[22.2062 26.6923 31.2905 36.2252 40.3749 44.861];
plot(surge_line(:, 2), surge_line(:, 1), 'r', 'LineWidth', 2);   % Surge line
plot(out.mass, out.P_ratio, 'r--', 'LineWidth', 2);         % Mass flow vs PR curve

plot(out.mass(end), out.P_ratio(end),'bo','LineWidth', 2);
plot(out.mass(1), out.P_ratio(1),'bo','LineWidth', 2);
text(out.mass(end), out.P_ratio(end), 'finish');
text(out.mass(1), out.P_ratio(1), 'start');

title('Compressor map and operating routine');
xlabel('Mass flow (kg/s)')
ylabel('Pressure ratio')
grid on;
hold off;

% %% T-S diagram
% subplot(1,2,2);
% % Define temperature range in Kelvin
% Tmin = 220;                           % Minimum temperature 
% Tmax = 304.12;                        % Maximum temperature
% T_range = linspace(Tmin, Tmax, 200);  % Create a temperature range
% 
% % Preallocate arrays for entropy and temperature
% s_liquid = zeros(size(T_range));     % Entropy of liquid phase
% s_vapor = zeros(size(T_range));      % Entropy of vapor phase
% 
% % Loop over the temperature range to get saturation properties
% for m = 1:length(T_range)
%     T = T_range(m);
% 
%     % Get the entropy of the liquid and vapor at saturation
%     s_liquid(m) = PropsSI('S','T',T,'Q',0,'CO2'); % Saturated liquid (Q=0)
%     s_vapor(m)  = PropsSI('S','T',T,'Q',1,'CO2'); % Saturated vapor (Q=1)
% end
% 
% % Concatenate entropy and temperature values to form a complete curve
% s_complete = [s_liquid, fliplr(s_vapor)]; % Fliplr flips the array from left to right
% T_complete = [T_range, fliplr(T_range)];
% 
% % Plot the complete saturation curve on T-s diagram
% plot(s_complete, T_complete, 'b-', 'LineWidth', 2,'DisplayName','Saturation Line');
% hold on
% 
% % T-s diagram of PTES
% T_array = [
% out.Comp_T_out(end)
% out.HHX_CO2_out(end)
% out.REC_hot_out(end)
% out.Exp_out(end)
% out.Cooler_out(end)
% out.CHX_CO2_out(end)
% out.REC_cold_out(end)]';
% P_array = [repelem(out.P_high(end)*1e3,3),repelem(out.P_low(end)*1e3,4)];
% P1       = out.P_low(end)*1e3;
% P2       = out.P_high(end)*1e3;
% S_array = zeros(1,length(T_array));
% 
% for i = 1:length(T_array)
%     S_array(i) = PropsSI('S', 'T', T_array(i), 'P', P_array(i), 'CO2');
% end
% % Draw the high pressure isobar
% T_min_h   = min(T_array(1:3));
% T_max_h   = max(T_array(1:3));
% T_range_h = linspace(T_min_h, T_max_h, 100);
% S_high    = arrayfun(@(T) PropsSI('S', 'T', T, 'P', P2, 'CO2'), T_range_h);  % Calculate entropy values along the high-pressure isobar
% 
% % Draw the low pressure isobar
% T_min_l   = min(T_array(4:end));
% T_max_l   = max(T_array(4:end));
% T_range_l = linspace(T_min_l, T_max_l, 100);
% S_low     = arrayfun(@(T) PropsSI('S', 'T', T, 'P', P1, 'CO2'), T_range_l);  % Calculate entropy values along the low-pressure isobar 
% 
% % Plot the T-s diagram, including the process line and isobars
% hold on;
% 
% n = 15;
% i = 1.4633748931767434;
% b = (P2*1.2-P1*0.85)/(n^i);
% P_values = P1*0.85+(linspace(1,15,15)).^i*b;
% T_values = linspace(250,950,50);
% for i = P_values
%     s_values = PropsSI('S', 'T', T_values,'P',i,'CO2');
%     plot(s_values,T_values,LineWidth=1.1, color=[1 0 0 0.4],HandleVisibility='off')
% end
% 
% % Plot the process line of compression and expansion
% plot([S_array,S_array(1)], [T_array,T_array(1)], 'ro', 'LineWidth', 2, 'DisplayName', 'Process data');
% plot([S_array(3:4)], [T_array(3:4)], 'r-', 'LineWidth', 2, 'HandleVisibility', 'off');
% plot([S_array(end),S_array(1)], [T_array(end),T_array(1)], 'r-', 'LineWidth', 2, 'HandleVisibility', 'off');
% 
% % Plot the low-pressure isobar
% plot(S_low, T_range_l, 'r-', 'LineWidth', 2, 'DisplayName', ['Isobar at P_{low}']);
% 
% % Plot the high-pressure isobar
% plot(S_high, T_range_h, 'r-', 'LineWidth', 2, 'DisplayName', ['Isobar at P_{high}']);
% 
% % Add labels, title, and legend
% xlabel('Entropy (J/kg.K)');
% ylabel('Temperature (K)');
% title('T-s Diagram with Isobars');
% legend('Location', 'best');
% grid on;
% hold off;