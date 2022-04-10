function [h_hat, h_basis, h_coef] = cs_est(rx_H, b_cb, num, basis_)
% OMP 기반의 Compressed sensing(CS)를 이용하여 신호를 분석하는 함수
% rx_H: beam sweeping을 통해 수신한 파일럿 채널
% b_cb: beam sweeping에 사용된 beam set
% basis: OMP에 사용되는 기저 벡터
% h_hat: 추정된 채널 행렬
% h_basis: 추정에 사용된 기저 벡터
% h_coef: 기저벡터를 통해 추정된 계수

% OMP를 위한 반복 횟수
if nargin < 3, num = 2; end

% 변수 초기화
[fft_len, N_rx, ~] = size(rx_H);
[N_tx, b_num] = size(b_cb);

% mode: 1 - 임의 기저벡터 생성, 2 - 주어진 기저벡터 사용
if nargin < 4
    mode = 1;
    M = b_num;
    h_basis = exp(2j*pi * (0:N_tx-1).' * (-0.5 : 1/M : 0.5 - 1/M) );
else
    mode = 2;
    [~, M] = size(basis_);
    h_basis = basis_;
end

% 빔 코드북과 기저벡터를 이용하여 탐색을 위한 행렬을 계산
phi = (b_cb.' * h_basis);


% 각 부반송파 별로 채널을 추정
h_hat = zeros(fft_len, N_rx, N_tx);
h_coef = zeros(fft_len, M);
for k = 1:fft_len
    
    % 수신 행렬을 직렬화
    H_(1:N_rx, 1:b_num) = rx_H(k,:,:);
    r = H_(:);
    
    % 3개의 기저벡터를 이용하여 채널을 추정
    for j = 1:num
        % 내적을 통한 성분 크기 비교 및 감산
        res = phi' * r;
        [~, idx] = max(abs(res));
        h_coef(k,idx) = res(idx) / sum( abs(phi(:,idx) ).^2 );
        r = r - (h_coef(k,idx) * phi(:,idx));
    end
    
    % 추정된 매개변수를 기반으로 채널을 복원
    h_hat(k,:,:) = h_basis * h_coef(k,:).';

end


