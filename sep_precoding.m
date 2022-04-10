function W = sep_precoding(t_H, N_rf, mode)
% Sspideh 논문의 MU-MISO를 구현한 빔포머(디지털 - ZF)
% t_H(필수): 2차원 채널 행렬(수신안테나, 송신안테나 또는 3차원 채널 행렬(부반송파, 수신안테나, 송신안테나)
% N_rf(필수): RF chain의 수 (N_rf >= N_s);
% mode(선택): 1은 fully-connected(default), 2는 partially-connected를 의미함
% W: 아날로그 프리코더


% 기본 변수 설정 (dim == 3: selective fading, mode = 2: flat fading)
dim = size(t_H);
if length(dim) == 3
    [fft_len, ~, N_tx] = size(t_H);
    
    % 채널 행렬의 부반송파 평균 계산
    A = 0;
    for i = 1:fft_len
        H(:,:) = t_H(i,:,:);
        A = A + (H' * H); 
    end
    A = A / fft_len;
    
elseif dim == 2
    [~, N_tx] = size(t_H);
    H = t_H;
    A = H' * H;
end

% 아날로그 구조 설정
if nargin < 3, mode = 1; end
if mode == 2
    B = zeros(N_tx, N_rf);
    s_num = N_tx / N_rf;
    for j = 1:N_rf, B((j-1)*s_num+1 : j*s_num, j) = ones(s_num,1); end
    N = s_num;
    
elseif mode == 1
    B = ones(N_tx, N_rf);
    N = N_tx;
end


% 최적화 변수 초기화(각 최적화에서 이전과의 차이 및 이전 행렬)
W = ones(N_tx, N_rf);
a = A * W;
diff_W = 100;
prev_W = W;

% 아날로그 프리코딩 행렬의 수렴을 확인
while diff_W > 4
    
    % 아날로그 프리코딩 행렬의 열 최적화
    for j = 1:N_rf
        
        % f_opt 계산
        a_s = 0;
        for l = 1:N_rf
            a_s = a_s + ( a(:,l) * a(:,l)' );
            if l == j, continue; end
        end
        
        % f_opt 정규화
        fopt = inv( a_s + eye(N_tx) ) * a(:,j);
        fopt = fopt / sqrt( trace( fopt * fopt' ) );
        
        % W와 a의 열벡터 업데이트
        W(:,j) = exp( 1j * angle(B(:,j).*fopt) ) / N;
        a(:,j) = A * W(:,j);
    end
    
    % 아날로그 프리코딩 행렬 차이 계산
    diff_W = sum( sum( abs(prev_W - W) ) );
    prev_W = W;
end


end