clear; clc; close all;

Nt = 2;  % Number of transmit antennas
Nr = 2;  % Number of receive antennas
mod_order = 4;  % QPSK
SNR_dB = 0:1:30;
num_bits = 1e5;

ber_stbc = zeros(size(SNR_dB));
ber_zf   = zeros(size(SNR_dB));
ber_mmse = zeros(size(SNR_dB));

for snr_idx = 1:length(SNR_dB)
    snr_db = SNR_dB(snr_idx);
    snr_lin = 10^(snr_db/10);
    noise_var = 1/snr_lin;

    % Generate bits and modulate
    bits = randi([0 1], num_bits, 1);
    symbols = qpsk_mod(bits);

    % Mode 1: Alamouti STBC (2x1)
    [symbol_pairs, tx_block] = alamouti_encode(symbols);
    ber_stbc(snr_idx) = mimo_stbc_simulate(symbol_pairs, tx_block, noise_var);

    % Mode 2: Spatial Multiplexing
    [ber_zf(snr_idx), ber_mmse(snr_idx)] = mimo_spatial_multiplexing(symbols, snr_lin);
end

%  BER vs. SNR for 2x1 Alamouti (QPSK over Rayleigh Fading)
valid_idx = find(ber_stbc > 0);  % 找出非零的 index
SNR_dB_valid = SNR_dB(valid_idx);
ber_valid = ber_stbc(valid_idx);

% 計算理論 BER（Alamouti QPSK over Rayleigh fading）

SNR_dB_theory = SNR_dB;  
SNR_lin = 10.^(SNR_dB_theory / 10);
Pb_theory = 0.5 * (1 - sqrt(SNR_lin ./ (1 + SNR_lin)));


figure;
semilogy(SNR_dB_valid, ber_valid, 'o-', 'DisplayName', 'Simulated BER');
hold on;
semilogy(SNR_dB_theory, Pb_theory, '--k', 'DisplayName', 'Theoretical BER');
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs. SNR for 2x1 Alamouti QPSK over Rayleigh Fading');
legend('Location', 'southwest');
grid on;



% BER vs. SNR for 2x2 Spatial Multiplexing
figure;
semilogy(SNR_dB_valid, ber_valid, 'o-', 'DisplayName', 'Simulated BER');
hold on;
semilogy(SNR_dB, ber_zf, '-s', SNR_dB, ber_mmse, '-x');
legend('Alamouti STBC','ZF', 'MMSE');
xlabel('SNR (dB)');
ylabel('BER');
title('BER vs. SNR');
grid on;

%% Throughput
rate_stbc = 0.5;   %Alamouti STBC
tp_stbc = rate_stbc * (1 - ber_stbc);

rate_sm = 2;  %Spatial Multiplexing
tp_zf   = rate_sm * (1 - ber_zf);
tp_mmse = rate_sm * (1 - ber_mmse);



%  Throughput vs. SNR for 2x2 Spatial Multiplexing
figure;
plot(SNR_dB, tp_stbc, 'o-',SNR_dB, tp_zf, '-s', SNR_dB, tp_mmse, '-x');
legend('Alamouti STBC','ZF', 'MMSE');
xlabel('SNR (dB)');
ylabel('Throughput (bits/symbol)');
title('Throughput vs. SNR for 2x2 MIMO');
grid on;



function symbols = qpsk_mod(bits)
    bits = reshape(bits, [], 2);
    symbol_map = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2);
    idx = bi2de(bits, 'left-msb') + 1;
    symbols = symbol_map(idx).';
end


function [symbol_pairs, tx_block] = alamouti_encode(symbols)
    symbols = reshape(symbols, 2, []);
    symbol_pairs = symbols;
    tx_block = [symbols(1,:); symbols(2,:); ...
               -conj(symbols(2,:)); conj(symbols(1,:))];
end


function ber = mimo_stbc_simulate(symbol_pairs, tx_block, noise_var)
    num_blocks = size(tx_block, 2);
    errors = 0;
    total_bits = 0;

    for k = 1:num_blocks
        h = (randn(2,1) + 1j*randn(2,1)) / sqrt(2);  % 2x1 channel

        s1 = symbol_pairs(1,k);
        s2 = symbol_pairs(2,k);

        x1 = s1;
        x2 = s2;
        x3 = -conj(s2);
        x4 =  conj(s1);

        % Transmit using Alamouti
        n1 = sqrt(noise_var/2)*(randn + 1j*randn);
        n2 = sqrt(noise_var/2)*(randn + 1j*randn);

        y1 = h(1)*x1 + h(2)*x2 + n1;
        y2 = h(1)*x3 + h(2)*x4 + n2;

        r1 = conj(h(1))*y1 + h(2)*conj(y2);  % 對應 s1
        r2 = conj(h(2))*y1 - h(1)*conj(y2);  % 對應 s2
       
        s_hat = [r1; r2] / (norm(h)^2);

        bits_rx = qpsk_demod(s_hat);
        bits_tx = qpsk_demod([s1; s2]);
        errors = errors + sum(bits_rx ~= bits_tx);
        total_bits = total_bits + length(bits_tx);
    end

    ber = errors / total_bits;
end




function [ber_zf, ber_mmse] = mimo_spatial_multiplexing(symbols, snr)
    symbols = reshape(symbols, 2, []);
    num_sym = size(symbols, 2);
    errors_zf = 0;
    errors_mmse = 0;
    noise_var = 1/snr;

    for k = 1:num_sym
        H = (randn(2,2) + 1j*randn(2,2)) / sqrt(2); % 建立 2x2 Rayleigh 通道
        x = symbols(:,k);
        n = sqrt(noise_var/2)*(randn(2,1)+1j*randn(2,1)); % 加入複數高斯雜訊
        y = H * x + n;  % 接收訊號

        x_zf = pinv(H) * y;
        x_mmse = H' * inv(H*H' + noise_var*eye(2)) * y;
         % 與原始訊號比較，計算錯誤位元數
        errors_zf = errors_zf + sum(qpsk_demod(x_zf) ~= qpsk_demod(x));
        errors_mmse = errors_mmse + sum(qpsk_demod(x_mmse) ~= qpsk_demod(x));
    end

    ber_zf = errors_zf / (2 * num_sym);
    ber_mmse = errors_mmse / (2 * num_sym);
end



function bits = qpsk_demod(symbols)
    ref = [1+1j, -1+1j, -1-1j, 1-1j] / sqrt(2);
    bits = zeros(2 * length(symbols), 1);
    for i = 1:length(symbols)
        [~, idx] = min(abs(symbols(i) - ref));
        bits(2*i-1:2*i) = de2bi(idx - 1, 2, 'left-msb').';
    end
end

