function W = sep_precoding(t_H, N_rf, mode)
% Sspideh ���� MU-MISO�� ������ ������(������ - ZF)
% t_H(�ʼ�): 2���� ä�� ���(���ž��׳�, �۽ž��׳� �Ǵ� 3���� ä�� ���(�ιݼ���, ���ž��׳�, �۽ž��׳�)
% N_rf(�ʼ�): RF chain�� �� (N_rf >= N_s);
% mode(����): 1�� fully-connected(default), 2�� partially-connected�� �ǹ���
% W: �Ƴ��α� �����ڴ�


% �⺻ ���� ���� (dim == 3: selective fading, mode = 2: flat fading)
dim = size(t_H);
if length(dim) == 3
    [fft_len, ~, N_tx] = size(t_H);
    
    % ä�� ����� �ιݼ��� ��� ���
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

% �Ƴ��α� ���� ����
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


% ����ȭ ���� �ʱ�ȭ(�� ����ȭ���� �������� ���� �� ���� ���)
W = ones(N_tx, N_rf);
a = A * W;
diff_W = 100;
prev_W = W;

% �Ƴ��α� �����ڵ� ����� ������ Ȯ��
while diff_W > 4
    
    % �Ƴ��α� �����ڵ� ����� �� ����ȭ
    for j = 1:N_rf
        
        % f_opt ���
        a_s = 0;
        for l = 1:N_rf
            a_s = a_s + ( a(:,l) * a(:,l)' );
            if l == j, continue; end
        end
        
        % f_opt ����ȭ
        fopt = inv( a_s + eye(N_tx) ) * a(:,j);
        fopt = fopt / sqrt( trace( fopt * fopt' ) );
        
        % W�� a�� ������ ������Ʈ
        W(:,j) = exp( 1j * angle(B(:,j).*fopt) ) / N;
        a(:,j) = A * W(:,j);
    end
    
    % �Ƴ��α� �����ڵ� ��� ���� ���
    diff_W = sum( sum( abs(prev_W - W) ) );
    prev_W = W;
end


end