function W = steer_precoding(fc, ant, angle, mode)
% 송수신 각도를 통해 Beamforming matrix를 계산하는 함수
% ant_num: 안테나의 배열 [ 수평방향의 안테나 수, 수직방향의 안테나 수, 수평방향의 안테나 간격/lamda, 수직방향의 안테나 간격/lamda ]
% rx_angle: 채널의 송수신 각도 [theta; phi]
% Nrf: RF chain (전송 스트림)의 수
% mode: 빔포머 구조(1: All connected, 2: Sub connected)

% 변수 초기화
if nargin < 4, mode = 1; end
if ant(1)*ant(2) == 1, mode = 2; end
[~, Nrf] = size(angle);

% 기본 변수 설정
%==========================================================
% 배열 변수
c = 3e8;                        % 빛의 속도
lamda = c/fc;                   % 신호의 파장
k = 2*pi / lamda;               % 파수
dy = ant(3) * lamda;            % 송신 안테나 간 거리 (가로)
dz = ant(4) * lamda;            % 송신 안테나 간 거리 (세로)
N = ant(1) * ant(2);

% 안테나 위치 행렬
temp1 = repmat(0:ant(1)-1, ant(2), 1);
temp2 = repmat(0:ant(2)-1, 1, ant(1));
ant_mat = [ zeros(N,1) temp1(:)*dy (temp2.')*dz];

% 좌표계 변환 함수
trans_f = @(t_theta, t_phi) [ sin(t_theta).*cos(t_phi); sin(t_theta).*sin(t_phi); cos(t_theta) ];


% 빔 matrix 계산 
if mode == 1
    % All connected array를 사용하는 경우
    
    W = zeros(N, Nrf);
    for i = 1:Nrf
        theta = angle(1,i);    phi = angle(2,i);
        W(:,i) = exp( -1j*k* ant_mat * trans_f(theta, phi) );
    end
    
elseif mode == 2
    % Sub connected array를 사용하는 경우
    
    W = zeros(N*Nrf, Nrf);
    for i = 1:Nrf
        theta = angle(1,i);    phi = angle(2,i);
        W(1+(i-1)*N : i*N, i) = exp( -1j*k* ant_mat * trans_f(theta, phi) );
    end
end


end