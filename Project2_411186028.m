close all;

Fs=1000;
fD=5;
L=3;

N = 1500;
N_lms=1500;
N_RLS = 200;
M = 8;
ite500=N+M;
ite1500 = N_lms+M;
ite200=N_RLS+M;

n = 1:3;  %列向量 [1 2 3]
chan = comm.RayleighChannel(...
    'SampleRate', Fs, ...
    'MaximumDopplerShift', fD, ...
    'PathDelays', (0:L-1)/Fs, ...
    'AveragePathGains', zeros(1,L));

x = 2 * randi([0 1], N, 1) - 1;


delay = 5; 
noise_var = 0.001;
runs = 500;        %Monte Carlo runs
mu = 0.005;        %LMS initial parameter
delta = 250;       %RLS initial parameter
lambda = 0.98;

%% 模擬出三條多徑通道
chan1 = comm.RayleighChannel( ...
    'SampleRate', Fs, ...
    'MaximumDopplerShift', fD, ...
    'PathDelays', (0:L-1)/Fs, ...
    'AveragePathGains', zeros(1, L), ...
    'PathGainsOutputPort', true);
x = ones(N,1);
[y, H] = chan1(x);
n = 1:size(H,1);
figure;
stem(n, abs(H(:,1)), 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineStyle', 'none', 'DisplayName', 'Path 1'); hold on;
stem(n, abs(H(:,2)), 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineStyle', 'none', 'DisplayName', 'Path 2');
stem(n, abs(H(:,3)), 'o', 'MarkerSize', 3, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'LineStyle', 'none', 'DisplayName', 'Path 3');

xlabel('Time Index n');
ylabel('|h_k[n]|');
title('Sample Paths of Rayleigh Fading Coefficients |h_k[n]|');
legend('Location', 'northeast');
grid on;
%% 雜訊 v[n]
v = sqrt(noise_var) * randn(N,1);
figure;
plot(v);
xlabel('Time Index n');
ylabel('Noise v[n]');
title('Sample Path of AWGN Noise v[n]');
grid on;
%% ==== Delay Sweep 測試 ====
delay_range = 3:9;
runs_per_delay = 10;
final_mse = zeros(size(delay_range));
final_std = zeros(size(delay_range)); 

for i = 1:length(delay_range)
    d = delay_range(i);
    mse_set = zeros(runs_per_delay, 1);
    
    for r = 1:runs_per_delay
        mse_tmp = lms(runs, mu, N+M, M, chan, noise_var, d);
        mse_set(r) = mean(mse_tmp(end-100:end));
    end

    final_mse(i) = mean(mse_set);  % 對多次的結果再取平均
    final_std(i) = std(mse_set);
end

figure;
errorbar(delay_range, final_mse, final_std, '-o', 'LineWidth', 2);
xlabel('Delay value');
ylabel('Steady-state MSE');
title('Delay Sweep: Steady-state MSE vs. Delay');
grid on;
%% %different mu
%mse_lms1 = lms(runs, 0.075, ite1500, M, chan, noise_var, delay);
mse_lms2 = lms(runs, 0.0075, ite1500, M, chan, noise_var, delay);
mse_lms3 = lms(runs, 0.025, ite1500, M, chan, noise_var, delay);
mse_lms4 = lms(runs, 0.01, ite1500, M, chan, noise_var, delay);

%semilogy(mse_lms1(M:ite1500))
%hold on;
semilogy(mse_lms2(M:ite1500))
hold on;
semilogy(mse_lms3(M:ite1500))
hold on;
semilogy(mse_lms4(M:ite1500))
hold on;

xlim([0 N_lms])
%legend(['\mu = ',num2str(0.075)],['\mu = ',num2str(0.0075)],['\mu = ',num2str(0.025)],['\mu = ',num2str(0.01)]);
legend(['\mu = ',num2str(0.0075)],['\mu = ',num2str(0.025)],['\mu = ',num2str(0.01)]);
xlabel('Number of adaptation cycles with LMS, n')
ylabel('MSE')
title('Effect of Step Size μ on LMS Performance');
hold off;
grid on;
%% %different lambda of RLS
mse_rlsL1 = rls(runs, delta, 1.0, ite200, M, chan, noise_var, delay);
mse_rlsL2 = rls(runs, delta, 0.9, ite200, M, chan, noise_var, delay);
mse_rlsL3 = rls(runs, delta, 0.8, ite200, M, chan, noise_var, delay);
mse_rlsL4 = rls(runs, delta, 0.7, ite200, M, chan, noise_var, delay);

figure;
semilogy(mse_rlsL1(M:ite200))
hold on;
semilogy(mse_rlsL2(M:ite200))
hold on;
semilogy(mse_rlsL3(M:ite200))
hold on;
semilogy(mse_rlsL4(M:ite200))
xlim([0 N_RLS])
legend(['λ = ',num2str(1)],['λ = ',num2str(0.9)],['λ = ',num2str(0.8)],['λ = ',num2str(0.7)]);
xlabel('Number of adaptation cycles with RLS , n')
ylabel('MSE')
title('Effect of Forgetting Factor λ on RLS Performance');
hold off;
grid on;
%% %LMS vs RLS
mse_lms = lms(runs, mu, ite500, M, chan, noise_var, delay);
mse_rls = rls(runs, delta, lambda, ite500, M, chan, noise_var, delay);

figure;
semilogy(mse_lms(M:ite500))
hold on;
semilogy(mse_rls(M:ite500))
xlim([0 N])
legend('LMS','RLS');
xlabel('Number of adaptation cycles, n')
ylabel('MSE')
title('Comparison of LMS and RLS Convergence ')
hold off;
grid on;
%% % LMS tap-weight trajectories
[mse, W_record] = lms_recorded(runs, mu, N, M, chan, noise_var, delay);

figure;
disp(max(abs(W_record(:))))
plot(W_record')
xlabel('Number of adaptation cycles, n')
ylabel('Tap Weights')
title('LMS Tap-weight Trajectories')

legendStrings = cell(1, M);
for i = 1:M
    legendStrings{i} = ['w_', num2str(i-1)];
end
legend(legendStrings, 'Location', 'best');
hold off;
grid on;
%% % calculate steady-state misadjustment
steady_start = 300;
steady_end = ite500;
misadjust_lms = mean(mse_lms(steady_start:steady_end));
misadjust_rls = mean(mse_rls(steady_start:steady_end));

figure;
bar([misadjust_lms, misadjust_rls]);
set(gca, 'XTickLabel', {'LMS', 'RLS'});

ylabel('Steady-state MSE');
title('Steady-state Misadjustment Comparison');
grid on;

%% BER VS SNR
SNR_dB = 0:2:20;  % SNR 以 dB 表示
BER_LMS = zeros(size(SNR_dB));
BER_RLS = zeros(size(SNR_dB));

for idx = 1:length(SNR_dB)
    snr_db = SNR_dB(idx);
    snr_linear = 10^(snr_db/10);
    noise_var_local = 1 / snr_linear;  % 根據 SNR 決定雜訊強度

    bit_errors_LMS = 0;
    bit_errors_RLS = 0;
    total_bits = 0;

    for r = 1:runs
        chan_snr = comm.RayleighChannel( ...
            'SampleRate', Fs, ...
            'MaximumDopplerShift', fD, ...
            'PathDelays', (0:L-1)/Fs, ...
            'AveragePathGains', zeros(1,L));
        x = 2 * randi([0 1], N, 1) - 1;
        reset(chan_snr);
        y = chan_snr(x);
        u = real(y) + sqrt(noise_var_local) * randn(N,1);

        w_lms = zeros(M,1);          % LMS
        for n = max(M, delay+1):N
            u_n = u(n:-1:n-M+1);
            y_hat = w_lms' * u_n;
            x_hat = sign(y_hat);  % 判決為 ±1
            d = x(n-delay);
            if x_hat ~= d
                bit_errors_LMS = bit_errors_LMS + 1;
            end
            e = d - y_hat;
            w_lms = w_lms + mu * u_n * e;
        end

        w_rls = zeros(M,1);           % RLS
        P = delta * eye(M);
        for n = max(M, delay+1):N
            u_n = u(n:-1:n-M+1);
            y_hat = w_rls' * u_n;
            x_hat = sign(y_hat);
            d = x(n-delay);
            if x_hat ~= d
                bit_errors_RLS = bit_errors_RLS + 1;
            end
            e = d - y_hat;
            k = P*u_n / (lambda + u_n'*P*u_n);
            w_rls = w_rls + k * e;
            P = (P - k*u_n'*P) / lambda;
        end

        total_bits = total_bits + (N - max(M, delay+1) + 1);  
    end

    BER_LMS(idx) = bit_errors_LMS / total_bits;
    BER_RLS(idx) = bit_errors_RLS / total_bits;
end
semilogy(SNR_dB, BER_LMS, '-o', 'DisplayName', 'LMS');
hold on;
semilogy(SNR_dB, BER_RLS, '-x', 'DisplayName', 'RLS');
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs. SNR for LMS and RLS');
legend;
grid on;

%% %LMS FOR tap-weight
function  [mse, W_record] = lms_recorded(runs, mu, N, M, chan, noise_var, delay)
    mse = zeros(N,1); 
    W_record = zeros(M, N);
    for i=1:runs
	    w=zeros(M,1);
	    x = 2 *randi([0 1], N + M-1, 1)-1;
        reset(chan); 
        y = chan(x);  
        y = real(y); 
        y = y(1:N) + sqrt(noise_var) * randn(N,1);
	    for n = max(M, delay+1):N
		    u_n=y(n:-1:n-M+1);
		    e=x(n-delay)-w'*u_n;
		    w=w+mu*u_n*e;
            W_record(:,n) = W_record(:,n) + w;
		    mse(n)=mse(n)+ abs(e)^2;
	    end
    end
    mse=mse/runs;
    W_record = W_record / runs;
end

%% %LMS 
function mse = lms(runs, mu, N, M, chan, noise_var, delay)
    mse = zeros(N,1); %產生一個長度 N 的「全是零的直向向量」。
    mse(1:max(M, delay+1)-1) = NaN;
    for i=1:runs
	    w=zeros(M,1);
	    x = 2 *randi([0 1], N, 1)-1;
        reset(chan);
        y = chan(x);
        y = real(y) + sqrt(noise_var) * randn(N, 1);
        for n = max(M, delay+1):N
		    u_n=y(n:-1:n-M+1);
		    e=x(n-delay)-w'*u_n;
		    w=w+mu*u_n*e;
		    mse(n)=mse(n)+abs(e)^2;
	    end
    end
    mse=mse/runs;
end

%% %RLS
function mse = rls(runs, delta, lambda, N, M, chan, noise_var, delay)
    mse = zeros(N,1);
    for i = 1:runs
        w = zeros(M,1);                        
        P = delta * eye(M);                    
        x = 2 * randi([0 1], N, 1) - 1;       

        reset(chan); 
        y = chan(x); 
        y = real(y) + sqrt(noise_var) * randn(N, 1);
        for n = max(M, delay+1):N
            u_n = y(n:-1:n-M+1);              
            e = x(n - delay) - w' * u_n;      
            k = (P * u_n) / (lambda + u_n' * P * u_n);  
            w = w + k * e;                    
            P = (P - k * u_n' * P) / lambda;  
            mse(n) = mse(n) + abs(e)^2;       
        end
    end
    mse = mse / runs;  
end






