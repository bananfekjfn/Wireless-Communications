
clear; clc; close all;

%% 參數設定 
fc = 2e9;
v1 =  0.0006;
v2 = 0.006;
c= 3e8;
fm1 = v1*fc/c;
fm2 = v2*fc/c;

fs  = 1024;
N   = 100;    % 可嘗試加大 N 以獲得更平滑的分佈，但運算時間也會增加

% 產生對應 fm1 的通道 
rayleigh_time1 = generateRayleighFading(N, fm1, fs);
rayleigh_time2 = generateRayleighFading(N, fm1, fs);
rayleigh_phase1 = generateRayleighphase(N, fm1, fs);
rayleigh_phase2 = generateRayleighphase(N, fm1, fs);

% 相位(度數)
rayleigh_phase_respone1 = angle((rayleigh_phase1 + rayleigh_phase2)/2) * 180/pi;

% 包絡 (轉 dB)
rayleigh_fading1 = sqrt((rayleigh_time1).^2 + (rayleigh_time2).^2);
rayleigh_fading1 = 10*log10(rayleigh_fading1 / mean(rayleigh_fading1));

% --- 產生對應 fm2 的通道 ---
rayleigh_time3 = generateRayleighFading(N, fm2, fs);
rayleigh_time4 = generateRayleighFading(N, fm2, fs);
rayleigh_phase3 = generateRayleighphase(N, fm2, fs);
rayleigh_phase4 = generateRayleighphase(N, fm2, fs);

rayleigh_phase_respone2 = angle((rayleigh_phase3 + rayleigh_phase4)/2) * 180/pi;

rayleigh_fading2 = sqrt((rayleigh_time3).^2 + (rayleigh_time4).^2);
rayleigh_fading2 = 10*log10(rayleigh_fading2 / mean(rayleigh_fading2));

%% (1) Rayleigh Magnitude Response
figure;
plot(rayleigh_fading1, 'Color','blue','LineWidth',1.5); hold on;
plot(rayleigh_fading2, 'Color',[1, 0, 0],'LineWidth',1.5);
xlim([0,length(rayleigh_fading2)]);
xlabel('Samples');
ylabel('Amplitude (dB)');
title('Rayleigh Magnitude Response');
legend(['fm = ', num2str(fm1),'Hz'], ['fm = ', num2str(fm2),'Hz']);
hold off;

%% (2) 繪製理論 U 型都卜勒頻譜 
figure;
Ns_plot = 512;
f_axis  = linspace(-0.05, 0.05, Ns_plot);

H1 = zeros(size(f_axis));
H2 = zeros(size(f_axis));
% 計算理論 U 型 (只有 |f| <= fd 才有值)
for k = 1:Ns_plot
    if abs(f_axis(k)) <= fm1
        H1(k) = 1/( pi*fm1 * sqrt( 1 - (f_axis(k)/fm1)^2 ));
    else
        H1(k) = 0;
    end
end


for k = 1:Ns_plot
    if abs(f_axis(k)) <= fm2
        H2(k) = 1 / (pi * fm2 * sqrt(1 - (f_axis(k)/fm2)^2));
    else
        H2(k) = 0;
    end
end
H1 = H1 / max(H1);%正規化

H2 = H2 / max(H2);
% 繪圖
plot(f_axis, H1, 'b-', 'LineWidth', 1.5, 'DisplayName', ['fm = ', num2str(fm1), ' Hz']); hold on;
plot(f_axis, H2, 'r', 'LineWidth', 1.5, 'DisplayName', ['fm = ', num2str(fm2), ' Hz']);
grid on;
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title(['U-Shaped Doppler Spectrum Comparison']);
xlim([-0.05, 0.05]);
legend(['fm = ', num2str(fm1), ' Hz'], ['fm = ', num2str(fm2), ' Hz']);

hold off;


%% (3) Ploting PDF with histogram
figure; 
%histcounts 將一組資料依照分箱（bins）進行統計，並計算每個區間內的數量或機率密度
[arg,bin]=histcounts(rayleigh_fading1,20,'Normalization', 'probability');
arg(end+1) =0;   %% 為了補齊最後一個 bin 的值，使繪圖完整
plot(bin,arg,'Color','blue','LineWidth',1.5);
hold on;
[arg,bin]=histcounts(rayleigh_fading2,20,'Normalization', 'probability');
arg(end+1) =0;
plot(bin,arg,'Color', [1, 0, 0],'LineWidth',1.5);

hold off;

xlabel('Data(dB)');
ylabel('Density');
title('Rayleigh Fading PDF  Approximation');
legend(['fm = ', num2str(fm1),'Hz'], ['fm = ', num2str(fm2),'Hz']);
xlim([-50,30]);   % 限制 X 軸範圍
hold off;


%% (4)Ploting Autocorrelation
figure;
% 產生「複數通道」
h_complex1 = generateRayleighFading(N, fm1, fs) + 1i*generateRayleighphase(N, fm1, fs);
h_complex2 = generateRayleighFading(N, fm2, fs) + 1i*generateRayleighphase(N, fm2, fs);

% 計算樣本自相關（normalized）
[acf_full1, lags1] = xcorr(h_complex1,'coeff');
[acf_full2, lags2] = xcorr(h_complex2,'coeff');

% 只取正半邊
idx_pos1 = find(lags1 >= 0);
acf_pos1 = acf_full1(idx_pos1);
tau1 = lags1(idx_pos1) / fs;  % 時間延遲 (秒)

idx_pos2 = find(lags2 >= 0);
acf_pos2 = acf_full2(idx_pos2);
tau2 = lags2(idx_pos2) / fs;  % 時間延遲 (秒)

% 計算理論 Jakes 自相關: R_theory = J_0(2 π fm τ)
R_theory1 = besselj(0, 2*pi*fm1*tau1); % 零階貝索函數
R_theory2 = besselj(0, 2*pi*fm2*tau2);

plot(tau1, real(acf_pos1), 'b', 'LineWidth', 1.5); hold on;
plot(tau1, R_theory1, 'magenta', 'LineWidth', 1.5);
plot(tau2, real(acf_pos2), 'r', 'LineWidth', 1.5);
plot(tau2, R_theory2, 'g', 'LineWidth', 1.5);
grid on;
xlabel('Time delay \tau (s)');
ylabel('Autocorrelation(0.004)');
legend(['fm = ', num2str(fm1), ' Hz'], ['J_0(2\pi f_m \tau), f_m = ', num2str(fm1), ' Hz'], ['fm = ', num2str(fm2), ' Hz'], ['J_0(2\pi f_m \tau), f_m = ', num2str(fm2), ' Hz']);
title('Sample Autocorrelation vs. Jakes Theory');
hold off;

function rayleigh_time = generateRayleighFading(N, fm, fs)  
   %Step 1：設定頻率網格
    df = 2*fm/(N-1);
    fd = -fm:df:fm;   % 頻率範圍，從 -fm 到 fm
   %Step 2：FFT 點數
    M = round(fs/df);   %在此頻率解析度下，對應的頻域取樣長度
    %Step 3：設計 Jakes U 型功率譜
    S_Ez = 1.5./(pi*fm*sqrt(1-(fd./fm).^2));
    S_Ez(1)   = 10*S_Ez(2);
    S_Ez(end) = 10*S_Ez(end-1);
    sqrt_S_Ez = sqrt(S_Ez);  % 濾波器 H(f) = sqrt(S(f))
    % Step 4： 複數高斯雜訊, I/Q 分支, Hermitian 對稱
    signal = (randn(1, (N/2)) + 1i*randn(1, (N/2)));  %signal 大小 = N/2，
    sample = [flip(conj(signal)), signal];  % 將一半的訊號取共軛並翻轉
    %Step 5：頻域濾波 + 零填充(用 sqrt_S_Ez 進行頻譜塑形 (Doppler 濾波))
    U_sample_padded = [zeros(1,round((M-N)/2)), sample.*sqrt_S_Ez, zeros(1,round((M-N)/2))];  
    %Step 6：組合頻域資料以對應 ifft 格式
    x = length(U_sample_padded); % x是整個頻域向量長度
    ifft_sample = [0, U_sample_padded(x/2+1:x), 0, U_sample_padded(1:(x/2))];
    %Step 7：IFFT 轉換回時域
    rayleigh_time = real(ifft(ifft_sample)); %取實部 real()確保輸出為實數訊號
end



function rayleigh_phase = generateRayleighphase(N, fm, fs)
    df = 2*fm/(N-1);
    fd = -fm:df:fm;
    M = round(fs/df);

    S_Ez = 1.5./(pi*fm*sqrt(1-(fd./fm).^2));
    S_Ez(1)   = 10*S_Ez(2);
    S_Ez(end) = 10*S_Ez(end-1);
    sqrt_S_Ez = sqrt(S_Ez);
    
    signal = (randn(1, N) + 1i*randn(1, N));   % 直接用長度 N 的複數雜訊
    U_sample_padded = [zeros(1,round((M-N)/2)), signal.*sqrt_S_Ez, zeros(1,round((M-N)/2))];
    x = length(U_sample_padded);
    ifft_sample = [0, U_sample_padded(x/2+1:x), 0, U_sample_padded(1:(x/2))];
    rayleigh_phase = ifft(ifft_sample);   % 輸出為複數 Q 分量
end




function rayleigh = outputRayleigh(N, fm, fs)
    % 產生兩組 Rayleigh 衰落通道的 I 分量（實數訊號）
    rayleigh_time1 = generateRayleighFading(N, fm, fs);
    rayleigh_time2 = generateRayleighFading(N, fm, fs);
     % 合併兩組 I 分量形成 Rayleigh 包絡（r(t) = sqrt(I^2 + Q^2)）
    rayleigh_fading1 = sqrt((rayleigh_time1).^2 + (rayleigh_time2).^2);
     % 正規化平均功率後轉換為 dB 單位
    rayleigh_fading1 = 10*log10(rayleigh_fading1 / mean(rayleigh_fading1));
    rayleigh = rayleigh_fading1;
end
