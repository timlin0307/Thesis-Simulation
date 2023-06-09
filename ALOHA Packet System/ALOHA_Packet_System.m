%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of the simulation
% clear all; % Initialise toutes les variables
close all; % Ferme toutes les fenetres ouvertes
clc; % Clear command window

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 1 : ALOHA SLOT PACKET SYSTEM
CR = 0:1:100; % Capture Ratio
beta = 1./(10.^(CR./10)); % Power Ratio
gamma_1 = 1./(1-beta); % Total Transmission (Synchronous)
sigma_max = beta + (1-beta).*exp(-gamma_1); % Actual Throughput (Maximum)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 2 : SATELLITE DELAY
sigma_sat = 0:0.01:1; % Actual Throughput
beta_sat = 0; % Power Ratio
m_sat = [1 2]; % [Synchronous Asynchronous]
gamma_sat = zeros(length(m_sat), length(sigma_sat)); % Total Transmission

% T_sat = 1400/(50*10^3); % Packet Time Length (50 kB channel using 1400 bit packets)
% C_sat = 0.27; % Propagation Delay
% D_sat = (C_sat + 6.*T_sat).*(gamma_sat./sigma_sat) - 5.*T_sat;
D_sat = zeros(length(m_sat), length(sigma_sat)); % Total Delay (Synchronous & Asynchronous)

for m = 1:length(m_sat)
    for s_sat = 1:length(sigma_sat)
        gamma_sat(m, s_sat) = fsolve(@(gamma)(beta_sat/m)*(1-exp(-gamma*m))+(1-beta_sat)*gamma*exp(-gamma*m)-sigma_sat(s_sat), 0);
    end
    D_sat(m, :) = 0.438.*(gamma_sat(m, :)./sigma_sat) - 0.14; % Total Delay
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 3 : ALOHA SLOT SYSTEM
CR_gnd_array = [0 1.5 3 6 Inf]; % Capture Ratio
sigma_gnd = 0:0.01:1; % Actual Throughput
m_gnd = 1; % Synchronous
gamma_gnd = zeros(1, length(sigma_gnd)); % Total Transmission

% T_gnd = 1400/(50*10^3); % Packet Time Length (50 kB channel using 1400 bit packets)
% C_gnd = 0.27; % Propagation Delay
% D_gnd = (2.*C_sat + 7.*T_sat).*(gamma_sat./sigma_sat) - (C_sat + 6.*T_sat);
D_gnd = zeros(length(CR_gnd_array), length(sigma_sat)); % Total Delay (Synchronous & Asynchronous)

for cr_gnd = 1:length(CR_gnd_array)
    beta_gnd = 1/(10^(CR_gnd_array(cr_gnd)/10)); % Power Ratio
    for s_gnd = 1:length(sigma_sat)
        gamma_gnd(s_gnd) = fsolve(@(gamma)(beta_gnd/m_gnd)*(1-exp(-gamma*m_gnd))+(1-beta_gnd)*gamma*exp(-gamma*m_gnd)-sigma_gnd(s_gnd), 0);
    end
    D_gnd(cr_gnd, :) = 0.196.*(gamma_gnd./sigma_gnd) - 0.168; % Total Delay
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIGURE 4 : ALOHA SLOT SYSTEM
CR_array = [0 1.5 3 6 Inf]; % Capture Ratio
gamma_x = 0:0.1:10; % Total Transmission
sigma_y = zeros(length(CR_array), length(gamma_x)); % Actual Throughput

for cr = 1:length(CR_array)
    beta_cr = 1/(10^(CR_array(cr)/10)); % Power Ratio
    for g = 1:length(gamma_x)
        sigma_y(cr, g) = beta_cr*(1-exp(-gamma_x(g))) + (1-beta_cr)*gamma_x(g)*exp(-gamma_x(g));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT THE FIGURE
% figure(1);
% loglog(CR, gamma_1, 'k-', CR, sigma_max, 'g-', CR, beta, 'm-'); % Result of analysis
% ylim([0.05 5])
% xlabel('Capture Ratio - DB');
% ylabel('γ, σ, β');
% legend('Total Transmission (γ)', 'Actual Throughput (σ)', 'Power Ratio (β)');
% title('Max Channel Utilization vs. FM CR for evenly populated area around STA');
% % hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% saveas(figure(1), 'figure1.jpg');

% figure(2);
% loglog(sigma_sat, D_sat(1, :), 'b-', sigma_sat, D_sat(2, :), 'r-'); % Result of analysis
% xlim([0.01 1])
% ylim([0.1 3])
% xlabel('Channel Utilization (σ)');
% ylabel('Delay - Sec.');
% legend('Asynchronous', 'Synchronous Slots');
% title('Asynchronous and Slots Systems for 50 KB Using 1400 Bit Packets');
% % hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% saveas(figure(2), 'figure2.jpg');

% figure(3);
% loglog(sigma_gnd, D_gnd(1, :), 'k-', sigma_gnd, D_gnd(2, :), 'g-', sigma_gnd, D_gnd(3, :), 'm-', sigma_gnd, D_gnd(4, :), 'r-', sigma_gnd, D_gnd(5, :), 'b-'); % Result of analysis
% xlim([0.01 1])
% ylim([0.01 1])
% xlabel('Actual Throughput (σ)');
% ylabel('Delay - Sec.');
% legend('CR = 0', 'CR = 1.5 DB', 'CR = 3 DB', 'CR = 6 DB', 'CR = ∞');
% title('Ground Use Delay for 1400 Bit / 50 KB');
% % hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% saveas(figure(3), 'figure3.jpg');

% figure(4);
% loglog(gamma_x, sigma_y(1, :), 'k-', gamma_x, sigma_y(2, :), 'g-', gamma_x, sigma_y(3, :), 'm-', gamma_x, sigma_y(4, :), 'r-', gamma_x, sigma_y(5, :), 'b-'); % Result of analysis
% xlim([0 10])
% ylim([0.1 2])
% xlabel('Total Transmission (γ)');
% ylabel('Actual Throughput (σ)');
% legend('CR = 0', 'CR = 1.5 DB', 'CR = 3 DB', 'CR = 6 DB', 'CR = ∞');
% title('σ vs. γ For Capture Ratios from 0-∞');
% % hold on; % Can keep former plotting and add new plotting
% grid on; % Display the grid
% saveas(figure(4), 'figure4.jpg');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%