function W = dft_precoding(t_H, N_rf)
% DFT 행렬을 통해 다중 사용자 채널을 위한 프리코딩 행렬을 계산
% t_H(필수): 2차원 채널 행렬(수신안테나, 송신안테나 또는 3차원 채널 행렬(부반송파, 수신안테나, 송신안테나)
% N_rf(필수): RF chain의 수 (N_rf >= N_s);
% W: 아날로그 프리코더

% 기본 변수 설정 (dim == 3: selective fading, mode = 2: flat fading)
dim = size(t_H);
if length(dim) == 3
    [fft_len, N_d, N_tx] = size(t_H);
elseif dim == 2
    [N_d, N_tx] = size(t_H);
    fft_len = 1;
end

% DFT 행렬 생성
f = -0.5 : 1/N_tx : 0.5 - 1/N_tx;
cb = exp( -2i*pi*( (0:N_tx-1).' * f ) );

% DFT 행렬로부터 벡터 생성
acc = 0;
for i = 1:fft_len
   if length(dim) == 3, H(1:N_d,:) = t_H(i,:,:);
   elseif length(dim) == 2, H = t_H;
   end
   
   % DFT 행렬의 적용 결과를 누적
   acc = acc + abs(H * cb).^2;
end

[~, idx] = max(acc,[],2);
W = cb(:,idx);

end