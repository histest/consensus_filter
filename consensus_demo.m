clear
clc
close all

%% target
% parameter initialization
T = 0.1;
nSteps = 200;
nStates = 2;
x_tgt_init = [3; 2];
x_tgt = zeros(nSteps, nStates);
sigma_v = 0.1;
F = [1, T; 0, 1];
Gamma = [T^2 / 2; T];
% target update
for i = 1 : nSteps
    if i == 1
        x_tgt(i, :) = (F * x_tgt_init + Gamma * normrnd(0, sigma_v) )';
    else
        x_tgt(i, :) = (F * (x_tgt(i - 1, :) )' + Gamma * normrnd(0, sigma_v) )';
    end
end
% testing
figure
subplot(1, 2, 1)
plot(1:nSteps, x_tgt(:, 1))
xlabel('time')
ylabel('target position')
subplot(1, 2, 2)
plot(1:nSteps, x_tgt(:, 2))
xlabel('time')
ylabel('target velocity')

%% measurement
% parameter initialization
H = [1, 0];
R = [1, 1, 1];
nSensors = 3;
Z = zeros(nSteps, nSensors);
% measurement update
for i = 1 : nSteps
    for j = 1 : nSensors
        Z(i, j) = H * (x_tgt(i, :) )' + normrnd(0, sqrt(R(j) ) );
    end
end
% testing
figure
for i = 1 : nSensors
    subplot(1, 3, i)
    plot(1:nSteps, x_tgt(:, 1), 1:nSteps, Z(:, i), '*');
    str = sprintf('position measurement from sensor %d', i);
    xlabel('time')
    ylabel(str)
end

%% filter
% parameter initialization
x_fltr_init = [2; 1];
P_fltr_init = diag( [0.1, 0.5] );
x_fltr = cell(nSteps, nSensors);
P_fltr = cell(nSteps, nSensors);
cm = [0 0 1; 1 0 0; 0 1 0];
Q = Gamma * (Gamma)' * (sigma_v)^2;
x_chnl_prev = cell(nSensors, nSensors);
P_chnl_prev = cell(nSensors, nSensors);
x_chnl_curr = cell(nSensors, nSensors);
P_chnl_curr = cell(nSensors, nSensors);
% main loop
for i = 1 : nSteps
    % local Kalman filter
    for j = 1 : nSensors 
        if i == 1
            x_KF_upd = x_fltr_init;
            P_KF_upd = P_fltr_init;
        else
            x_KF_upd = x_fltr{i - 1, j};
            P_KF_upd = P_fltr{i - 1, j};
        end
        [x_KF_pred, P_KF_pred] = KF_pred(x_KF_upd, P_KF_upd, F, Q);
        z = Z(i, j);
        [x_KF_upd, P_KF_upd] = KF_upd(z, x_KF_pred, P_KF_pred, H, R(j) );
        % update
        x_fltr{i, j} = x_KF_upd;
        P_fltr{i, j} = P_KF_upd; 
    end
    % distributed Kalman filter, DKF (consensus step 1)
    for j = 1 : nSensors
        for k = 1 : nSensors
            if cm(k, j) == 1
                z = Z(i, k);
                x_DKF = x_fltr{i, j};
                P_DKF = P_fltr{i, j};
                [x_DKF, P_DKF] = KF_upd(z, x_DKF, P_DKF, H, R(k) );
                % update
                x_fltr{i, j} = x_DKF;
                P_fltr{i, j} = P_DKF;
            end
        end
    end
    % channel filter (consensus step 2)
    nConsensus = 10;
    for j = 1 : nConsensus
        % communication
        for k = 1 : nSensors
            for l = 1 : nSensors
                if cm(l, k) == 1
                    x_chnl_prev{l, k} = x_chnl_curr{l, k};
                    P_chnl_prev{l, k} = P_chnl_curr{l, k};
                    x_chnl_curr{l, k} = x_fltr{i, l};
                    P_chnl_curr{l, k} = P_fltr{i, l};
                end
            end 
        end
        % backup local estimate
        x_fltr_backup = x_fltr;
        P_fltr_backup = P_fltr;
        % filter
        for k = 1 : nSensors
            x_1 = x_fltr_backup{i, k};
            P_1 = P_fltr_backup{i, k};
            for l = 1 : nSensors
                if cm(l, k) == 1
                    x_2 = x_chnl_curr{l, k};
                    P_2 = P_chnl_curr{l, k};
                    if j == 1
                        if i == 1
                            assert(isempty(x_chnl_prev{l, k}) )
                            x_p = F * x_fltr_init;
                            P_p = F * P_fltr_init * F' + Q;
                        else
                            x_p = F * x_chnl_prev{l, k};
                            P_p = F * P_chnl_prev{l, k} * F' + Q;
                        end
                    else
                        x_p = x_chnl_prev{l, k};
                        P_p = P_chnl_prev{l, k};
                    end
                    [x_CF, P_CF] = chnl_fltr(x_p, P_p, x_1, P_1, x_2, P_2);
                    % update
                    x_fltr{i, k} = x_CF;
                    P_fltr{i, k} = P_CF;
                end
            end
        end       
    end  
end
% testing
figure
for i = 1 : nSensors
    t = cell2mat(x_fltr(:, i) );
    t_p = t(1:2:end);
    t_v = t(2:2:end);
    subplot(2, nSensors, i)
    plot(1:nSteps, x_tgt(:, 1), 1:nSteps, t_p);
    xlabel('time')
    ylabel('position')
    str = sprintf('local estimation from sensor %d', i);
    legend('target', str)
    hold on
    subplot(2, nSensors, i + nSensors)
    plot(1:nSteps, x_tgt(:, 2), 1:nSteps, t_v);
    xlabel('time')
    ylabel('velocity')
    str = sprintf('local estimation from sensor %d', i);
    legend('target', str)
end
% bg = biograph(cm);
% view(bg)