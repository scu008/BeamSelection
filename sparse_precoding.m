function [W, sum_coef] = sparse_precoding(h_coef, Nrf, b_cb)
% sparse vector를 입력으로 받아 하이브리드 빔포밍을 위한 아날로그 빔포밍 행렬을 계산
% h_coef: sparse vector의 크기 vector
% Nrf: RF chain (전송 스트림)의 수
% 빔포밍 행렬 구성을 위한 빔 후보 행렬

% 초기 변수 설정
[N_d, N_tx] = size(h_coef);
rx_idx = cell(1,4);

% Sparse vector 요소들의 크기에 따른 정렬을 수행
comp_coef = sum(h_coef,1);
[~, idx] = sort(comp_coef, 'descend');

% 각 수신기의 경로와 그게 따른 합성 이득 구분 저장
for i = 1:N_d
    tmp_idx = find(h_coef(i,:));
    rx_idx{i} = [ tmp_idx; comp_coef(tmp_idx) ];
end

% 초기 빔 선택 후보군 설정
b_idx = idx(1:Nrf);

% rank 및 중복 선택 점검
esc_flg = zeros(1,N_d);
rep_idx = [];
for i = 1:Nrf
    
    % 각 수신기마다 해당하는 빔이 있는지 점검
    for j = 1:N_d
        rx_tmp = rx_idx{j};
        if sum(rx_tmp(1,:) == b_idx(i))
            esc_flg(j) = esc_flg(j) + 1; 
            
            % 중복이 있을 경우 기록
            if esc_flg(j) > 1, rep_idx = [i rep_idx]; end
        end
    end
end

% 해당하는 빔이 없는 수신기를 탐색
n_idx = find(esc_flg == 0);
for i = n_idx
    rx_tmp = rx_idx{i};
    [~, m_idx] = max(rx_tmp(2,:));
    b_idx(rep_idx(1)) = rx_tmp(1,m_idx);
    rep_idx(1) = [];
end


% 빔포밍 행렬 결정
W = b_cb(:,b_idx);
sum_coef = sum( comp_coef(b_idx) );








