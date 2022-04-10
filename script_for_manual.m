%% 기본 초기화
clc, clear


% 모델 초기화 (객체의 이름은 자류롭게 설정 ex. model_A, node_1, node_2)
model = SCM();

% 거리 비율 설정
model.distance_rate = 1;

% 다중 경로 수 설정
model.n_path = 7;

% Pathloss exponent 설정 (거리에 따른 전력 감쇄)
model.exp_beta = 3;

% LOS 환경 사용 여부 결정
% 1로 설정하면 절대좌표가 LOS 환경일 시 LOS 반영
% 0으로 설정하면 절대좌표와 관계없이 무조건 non_LOS 환경 반영
model.los = 1;

% 신호의 중심 주파수 설정(default는 이동통신 주파수 대역인 800MHz)
model.fc = 800e6;


%% SISO 신호 전송


% 신호 생성 
mod_type = 2;
fft_len = 64;
bit = randi([0 1], 1, fft_len * mod_type);
sym = base_mod(bit, mod_type);


% 모델 생성(node A와의 채널과 조금 더 거리가 먼 node B와의 채널을 생성)
% A 채널(distance_rate = 1 - default)
node_A = SCM();

% B 채널(distance_rate = 2)
node_B = SCM();
node_B.distance_rate = 2;

% 모델들을 벡터를 통해 관리
model = [node_A node_B];

% 채널 계수 생성
H_A = model(1).FD_channel(fft_len);
H_B = model(2).FD_channel(fft_len);


% 신호 전송 및 수신(FD_fading은 모든 모델에서 동일하게 작동)
y_A = awgn_noise(  model(1).FD_fading(sym, H_A),  20);
y_B = awgn_noise(  model(1).FD_fading(sym, H_B),  20);


%% MIMO 채널 모델
clc, clear


% 모델 생성
Nt = 4; Nr = 4;
model = SCM();
model.ant(Nr, Nt);

% 송신단과 수신단의 안테나 간격 설정 (4 * 파장)
% MIMO의 공간 다중화로 안정적인 성능을 위해서는
% 안테나 간격이 최소 반파장 이상은 벌어져야 함
model.tx_ant(3) = 4;
model.rx_ant(3) = 4;


% 신호 생성 
mod_type = 2;
fft_len = 64;
bit = randi([0 1], Nt, fft_len * mod_type);
sym = base_mod(bit, mod_type);

% 채널 계수 생성
H = model.FD_channel(fft_len);

% 신호 전송 및 수신
y = awgn_noise(  model.FD_fading(sym, H),  20);


%% Fast Fading Channel
clc, clear


% 모델 생성
Nt = 4; Nr = 4;
model = SCM();
model.ant(Nr, Nt);
model.los = 1;

% 반송파 및 샘플링 주파수 설정(반송파 주파수 = 6GHz, 샘플링 주파수 = 20MHz)
model.fc = 60e9;
model.fs = 20e6;

% 송신단과 수신단의 안테나 간격 설정 (4 * 파장)
model.tx_ant(3) = 4;
model.rx_ant(3) = 4;


% 신호 생성 
mod_type = 2;
fft_len = 1000;
bit = randi([0 1], Nt, fft_len * mod_type);
sym = base_mod(bit, mod_type);

% 채널 계수 생성(상대속도 300km/h)
H = model.FD_channel(fft_len, 300);

% 신호 전송 및 수신
y = awgn_noise(  model.FD_fading(sym, H),  20);


