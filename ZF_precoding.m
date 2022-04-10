function [ sym_hat, K, W ] = ZF_precoding(sym, H)
% ZF precoding을 수행하는 함수
% H: 통신 채널 행렬
% sym_hat: precoding이 적용된 송신 신호
% K: 각 부반송파마다 계산된 precoding 정규화 벡터
% rep: 아날로그 빔포밍을 위한 브렌치의 수

% 채널 행렬의 사이즈를 구한다.
[fft_len, Mr, N_tx] = size(H);
[N_s, fft_len] = size(sym);

% 반복 행렬 선언
K = zeros(1,fft_len);
W = zeros(fft_len,N_tx,N_s);
sym_hat = zeros(N_tx,fft_len);

% 모든 서브케리어에 대해 precoding 수행
for k = 1 : fft_len
    
    % 서브케리어 하나의 채널 행렬을 저장
    t_H(:,:) = H(k,:,:);
    
    % Precoding matrix 계산
    G = t_H' * inv( t_H * t_H' );
    
    % Precoding 수행
    sym_hat(:,k) = G * sym(:,k);
    W(k,:,:) = G;
    K(k) = sqrt( trace(G * G') );
end




