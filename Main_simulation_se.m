%% Full-digital (ZF)
clc, clear

% 1: Capacity, 2: BER, 3: per cap
mode = 3;

% 전송 파라미터 설정
fft_len = 64;               % OFDM 부반송파의 수
mod_type = 1;               % 변조 차수 ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 2;                % 수신기의 수 (수신기의 안테나는 1개)
tx_ant = 32;                % 기지국의 안테나 수
snr = -20:5:15;             % 전송 채널 SNR 범위
path = 3;
scatter = 10;
iter = 200;                 % 전송 반복 횟수

% 코드북 파라미터 설정
b_bit = log2(tx_ant);
h_bit = 6;
cs_dim = path;
Ntx = tx_ant;
load_codebook
h_cb = p_cscb;
% h_cb = c_cscb;
sdw_sig = 4;

% 설정 파라미터 기반의 프로그램 실행
rx_ant = 4;                 % 수신기의 개별 안테나 수
plot_format = 'm-.s';
ZF_quantization             % 실행 프로그램
hold on


%% Quantized Hybrid Beamforming (QHB): Beam Selection Based on DFT Matrix
clc, clear

% 1: Capacity, 2: BER
mode = 3;

% 전송 파라미터 설정
fft_len = 64;               % OFDM 부반송파의 수
mod_type = 1;               % 변조 차수 ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 2;                % 수신기의 수 (수신기의 안테나는 1개)
tx_ant = 32;                % 기지국의 안테나 수
snr = -20:5:15;             % 전송 채널 SNR 범위
path = 3;
scatter = 10;
iter = 200;                 % 전송 반복 횟수

% 코드북 파라미터 설정
b_bit = log2(tx_ant);
h_bit = 6;
cs_dim = path;
Ntx = tx_ant;
load_codebook
h_cb = p_cscb;
% h_cb = c_cscb;
sdw_sig = 4;

sel_node = rx_node;

% 설정 파라미터 기반의 프로그램 실행
num_rf = sel_node * 2;       % 기지국의 RF chain의 수 *(필수 조건: rx_node <= num_rf <= tx_ant)*
rx_ant = 4;                 % 수신기의 개별 안테나 수
plot_format = 'b-s';
Hybrid_quantization_dft     % 실행 프로그램
hold on


%% Quantized Hybrid Beamforming (QHB): SEP
clc, clear

% 1: Capacity, 2: BER
mode = 3;

% 전송 파라미터 설정
fft_len = 64;               % OFDM 부반송파의 수
mod_type = 1;               % 변조 차수 ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 2;                % 수신기의 수 (수신기의 안테나는 1개)
tx_ant = 32;                % 기지국의 안테나 수
snr = -20:5:15;             % 전송 채널 SNR 범위
path = 3;
scatter = 10;
iter = 200;                 % 전송 반복 횟수

% 코드북 파라미터 설정
b_bit = log2(tx_ant);
h_bit = 6;
cs_dim = path;
Ntx = tx_ant;
load_codebook
h_cb = p_cscb;
% h_cb = c_cscb;
sdw_sig = 4;

% 설정 파라미터 기반의 프로그램 실행
num_rf = rx_node * 2;       % 기지국의 RF chain의 수 *(필수 조건: rx_node <= num_rf <= tx_ant)*
rx_ant = 4;                 % 수신기의 개별 안테나 수
plot_format = 'r-.s';
Hybrid_quantization_sep     % 실행 프로그램
hold on


%% Beam steering method
clc, clear

% 1: Capacity, 2: BER
mode = 1;

% 전송 파라미터 설정
fft_len = 64;               % OFDM 부반송파의 수
mod_type = 2;               % 변조 차수 ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 4;                % 수신기의 수 (수신기의 안테나는 1개)
tx_ant = 32;                % 기지국의 안테나 수
snr = -5:5:20;              % 전송 채널 SNR 범위
path = 3;
scatter = 10;
iter = 400;                 % 전송 반복 횟수

% 설정 파라미터 기반의 프로그램 실행
rx_ant = 4;                 % 수신기의 개별 안테나 수
plot_format = 'r-o';
Hybrid_beamforming          % 실행 프로그램
hold on
