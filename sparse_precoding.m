function [W, sum_coef] = sparse_precoding(h_coef, Nrf, b_cb)
% sparse vector�� �Է����� �޾� ���̺긮�� �������� ���� �Ƴ��α� ������ ����� ���
% h_coef: sparse vector�� ũ�� vector
% Nrf: RF chain (���� ��Ʈ��)�� ��
% ������ ��� ������ ���� �� �ĺ� ���

% �ʱ� ���� ����
[N_d, N_tx] = size(h_coef);
rx_idx = cell(1,4);

% Sparse vector ��ҵ��� ũ�⿡ ���� ������ ����
comp_coef = sum(h_coef,1);
[~, idx] = sort(comp_coef, 'descend');

% �� ���ű��� ��ο� �װ� ���� �ռ� �̵� ���� ����
for i = 1:N_d
    tmp_idx = find(h_coef(i,:));
    rx_idx{i} = [ tmp_idx; comp_coef(tmp_idx) ];
end

% �ʱ� �� ���� �ĺ��� ����
b_idx = idx(1:Nrf);

% rank �� �ߺ� ���� ����
esc_flg = zeros(1,N_d);
rep_idx = [];
for i = 1:Nrf
    
    % �� ���ű⸶�� �ش��ϴ� ���� �ִ��� ����
    for j = 1:N_d
        rx_tmp = rx_idx{j};
        if sum(rx_tmp(1,:) == b_idx(i))
            esc_flg(j) = esc_flg(j) + 1; 
            
            % �ߺ��� ���� ��� ���
            if esc_flg(j) > 1, rep_idx = [i rep_idx]; end
        end
    end
end

% �ش��ϴ� ���� ���� ���ű⸦ Ž��
n_idx = find(esc_flg == 0);
for i = n_idx
    rx_tmp = rx_idx{i};
    [~, m_idx] = max(rx_tmp(2,:));
    b_idx(rep_idx(1)) = rx_tmp(1,m_idx);
    rep_idx(1) = [];
end


% ������ ��� ����
W = b_cb(:,b_idx);
sum_coef = sum( comp_coef(b_idx) );








