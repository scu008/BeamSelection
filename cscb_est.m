function [h_hat, b_idx, h_basis, h_coef] = cscb_est(rx_H, b_cb, h_cb, basis_, No)
% Compressed Sensing 기반의 채널 추정 후 양자화를 수행하는 함수
% rx_H: beam sweeping을 통해 수신한 파일럿 채널
% b_cb: beam sweeping에 사용된 beam set
% h_cb: effective channel feedback에 사용되는 코드북
% basis_: 기저벡터를 선택하기 위한 조건을 명시
% No: AWGN 잡음의 전력

% 변수 초기화
[fft_len, N_rx, ~] = size(rx_H);
[N_tx, b_num] = size(b_cb);
[num, ~] = size(h_cb);
[~, m_idx(:)] = sort( sum( sum(abs(rx_H).^2, 1), 2), 'descend' );

% 기저벡터 선택 방식을 구현
if nargin < 4
    M = b_num;
    h_basis = exp(2j*pi * (0:N_tx-1).' * (-0.5 : 1/M : 0.5 - 1/M) );
    
    % 빔 코드북과 기저벡터를 이용하여 탐색을 위한 행렬을 계산
    phi = (b_cb.' * h_basis);
    
elseif (basis_ == 'p') && (nargin < 5)
    M = num;
    h_basis = b_cb(:,m_idx(1:M)).^-1;
    
    % 빔 코드북과 기저벡터를 이용하여 탐색을 위한 행렬을 계산
    phi = (b_cb.' * h_basis);
    
elseif (basis_ == 'n')
    M = b_num;
    h_basis = exp(2j*pi * (0:N_tx-1).' * (-0.5 : 1/M : 0.5 - 1/M) );
    phi = (b_cb.' * h_basis + No*eye(M));
    
    % MMSE 솔루션을 사용하여 기저벡터 추출
    res = zeros(fft_len, N_rx * b_num);
    for k = 1:fft_len
        H_(1:N_rx, 1:b_num) = rx_H(k,:,:);
        res(k,:) = (phi' * H_(:)).';
    end
    
    % 기저벡터를 선별
    [~, m_idx(:)] = sort( sum(abs(res).^2, 1), 'descend' );
%     phi = phi(:, m_idx(1:num));
%     h_basis = h_basis(:, m_idx(1:num));
    
    phi(:, m_idx(num+1:end)) = zeros(size( phi(:, m_idx(num+1:end)) ));
    h_basis(:, m_idx(num+1:end)) = zeros(size( h_basis(:, m_idx(num+1:end)) ));
    
else
    [~, M] = size(basis_);
    h_basis = basis_;
    
    % 빔 코드북과 기저벡터를 이용하여 탐색을 위한 행렬을 계산
    phi = (b_cb.' * h_basis);
end
b_idx = m_idx(1);


% 각 부반송파 별로 채널을 추정
h_hat = zeros(fft_len, N_rx, N_tx);
h_coef = zeros(fft_len, M);
for k = 1:fft_len
    
    % 수신 행렬을 직렬화
    H_(1:N_rx, 1:b_num) = rx_H(k,:,:);
    r = H_(:);
    
    % 기저벡터를 이용하여 채널을 추정
    idx = zeros(1,num);
    for j = 1:num
        % 내적을 통한 성분 크기 비교 및 감산
        res = phi' * r;
        [~, idx(j)] = max(abs(res));
        h_coef(k,idx(j)) = res(idx(j)) / sum( abs(phi(:,idx(j)) ).^2 );
        r = r - (h_coef(k,idx(j)) * phi(:,idx(j)));
    end
    
    % 추정된 매개변수 및 코드북을 기반으로 채널을 복원
    u_coef = h_coef(k,idx).' / sqrt( sum(abs(h_coef(k,idx)).^2) );
    [~, h_idx] = min( sum( abs(u_coef - h_cb).^2 ) );
    h_hat(k,:,:) = h_basis(:,idx) * h_cb(:,h_idx) * sqrt( sum(abs(h_coef(k,idx)).^2) );

end







