%% Full-digital (ZF)
clc, clear

% 1: Capacity, 2: BER, 3: per cap
mode = 3;

% ���� �Ķ���� ����
fft_len = 64;               % OFDM �ιݼ����� ��
mod_type = 1;               % ���� ���� ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 2;                % ���ű��� �� (���ű��� ���׳��� 1��)
tx_ant = 32;                % �������� ���׳� ��
snr = -20:5:15;             % ���� ä�� SNR ����
path = 3;
scatter = 10;
iter = 200;                 % ���� �ݺ� Ƚ��

% �ڵ�� �Ķ���� ����
b_bit = log2(tx_ant);
h_bit = 6;
cs_dim = path;
Ntx = tx_ant;
load_codebook
h_cb = p_cscb;
% h_cb = c_cscb;
sdw_sig = 4;

% ���� �Ķ���� ����� ���α׷� ����
rx_ant = 4;                 % ���ű��� ���� ���׳� ��
plot_format = 'm-.s';
ZF_quantization             % ���� ���α׷�
hold on


%% Quantized Hybrid Beamforming (QHB): Beam Selection Based on DFT Matrix
clc, clear

% 1: Capacity, 2: BER
mode = 3;

% ���� �Ķ���� ����
fft_len = 64;               % OFDM �ιݼ����� ��
mod_type = 1;               % ���� ���� ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 2;                % ���ű��� �� (���ű��� ���׳��� 1��)
tx_ant = 32;                % �������� ���׳� ��
snr = -20:5:15;             % ���� ä�� SNR ����
path = 3;
scatter = 10;
iter = 200;                 % ���� �ݺ� Ƚ��

% �ڵ�� �Ķ���� ����
b_bit = log2(tx_ant);
h_bit = 6;
cs_dim = path;
Ntx = tx_ant;
load_codebook
h_cb = p_cscb;
% h_cb = c_cscb;
sdw_sig = 4;

sel_node = rx_node;

% ���� �Ķ���� ����� ���α׷� ����
num_rf = sel_node * 2;       % �������� RF chain�� �� *(�ʼ� ����: rx_node <= num_rf <= tx_ant)*
rx_ant = 4;                 % ���ű��� ���� ���׳� ��
plot_format = 'b-s';
Hybrid_quantization_dft     % ���� ���α׷�
hold on


%% Quantized Hybrid Beamforming (QHB): SEP
clc, clear

% 1: Capacity, 2: BER
mode = 3;

% ���� �Ķ���� ����
fft_len = 64;               % OFDM �ιݼ����� ��
mod_type = 1;               % ���� ���� ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 2;                % ���ű��� �� (���ű��� ���׳��� 1��)
tx_ant = 32;                % �������� ���׳� ��
snr = -20:5:15;             % ���� ä�� SNR ����
path = 3;
scatter = 10;
iter = 200;                 % ���� �ݺ� Ƚ��

% �ڵ�� �Ķ���� ����
b_bit = log2(tx_ant);
h_bit = 6;
cs_dim = path;
Ntx = tx_ant;
load_codebook
h_cb = p_cscb;
% h_cb = c_cscb;
sdw_sig = 4;

% ���� �Ķ���� ����� ���α׷� ����
num_rf = rx_node * 2;       % �������� RF chain�� �� *(�ʼ� ����: rx_node <= num_rf <= tx_ant)*
rx_ant = 4;                 % ���ű��� ���� ���׳� ��
plot_format = 'r-.s';
Hybrid_quantization_sep     % ���� ���α׷�
hold on


%% Beam steering method
clc, clear

% 1: Capacity, 2: BER
mode = 1;

% ���� �Ķ���� ����
fft_len = 64;               % OFDM �ιݼ����� ��
mod_type = 2;               % ���� ���� ex) 1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
rx_node = 4;                % ���ű��� �� (���ű��� ���׳��� 1��)
tx_ant = 32;                % �������� ���׳� ��
snr = -5:5:20;              % ���� ä�� SNR ����
path = 3;
scatter = 10;
iter = 400;                 % ���� �ݺ� Ƚ��

% ���� �Ķ���� ����� ���α׷� ����
rx_ant = 4;                 % ���ű��� ���� ���׳� ��
plot_format = 'r-o';
Hybrid_beamforming          % ���� ���α׷�
hold on
