function [h_hat, m_idx] = bsweep_est(rx_H, b_cb, h_cb_)
% Beam sweeping based channel feedback and estimation
% rx_H: beam sweeping을 통해 수신한 파일럿 채널
% b_cb: beam sweeping에 사용된 beam set
% h_cb: effective channel feedback에 사용되는 코드북

% mode: 1 - 선형 양자화, 2, 비선형 양자화
if ismatrix(h_cb_)
    mode = 1;
    [~, h_num] = size(h_cb_);
else
    mode = 2;
    [~, ~, h_num] = size(h_cb_);
end

% 변수 초기화
[fft_len, N_rx, b_num] = size(rx_H);
[N_tx, ~] = size(b_cb);
h_hat = zeros(fft_len, N_rx, N_tx);


% 각 부반송파 별로 채널을 추정
[~, m_idx] = max( sum( sum(abs(rx_H).^2, 1), 2) );
for k = 1:fft_len
    
    % 각 수신 안테나별로 채널을 추정
    H_(1:N_rx, 1:b_num) = rx_H(k,:,:);
    for r = 1:N_rx
        
        % 채널 벡터 정규화
        beta = sum( abs(H_(r,:)).^2 );
        H_(r,:) = H_(r,:) / sqrt(beta);
        
        % 선형 비선형 분류
        if mode == 1, h_cb = h_cb_;
        elseif mode == 2
            [~, b_idx] = max( abs(H_).^2 );
            h_cb(1:N_tx, 1:h_num) = h_cb_(b_idx,:,:);
            h_cb(1:N_tx, 1:h_num) = h_cb_(m_idx,:,:);
        end
        
        % 추정 행렬 정규화
        tmp = b_cb.' * h_cb;
        tmp = tmp ./ sqrt( sum( abs(tmp).^2, 1 ) );
        
        % 단위 채널 벡터 추정
        dist = sum( abs( H_(r,:).' - tmp ).^2, 1 );
        [~, idx] = min(dist);
        h_u = h_cb(:,idx).';
        
        % 채널 추정
        beta_h = beta / sum( abs(h_u * b_cb).^2, 2 );
        h_hat(k,r,:) = h_u * sqrt(beta_h);
    end
end




