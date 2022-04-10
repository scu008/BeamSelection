function [ sym_hat, K, W ] = ZF_precoding(sym, H)
% ZF precoding�� �����ϴ� �Լ�
% H: ��� ä�� ���
% sym_hat: precoding�� ����� �۽� ��ȣ
% K: �� �ιݼ��ĸ��� ���� precoding ����ȭ ����
% rep: �Ƴ��α� �������� ���� �귻ġ�� ��

% ä�� ����� ����� ���Ѵ�.
[fft_len, Mr, N_tx] = size(H);
[N_s, fft_len] = size(sym);

% �ݺ� ��� ����
K = zeros(1,fft_len);
W = zeros(fft_len,N_tx,N_s);
sym_hat = zeros(N_tx,fft_len);

% ��� �����ɸ�� ���� precoding ����
for k = 1 : fft_len
    
    % �����ɸ��� �ϳ��� ä�� ����� ����
    t_H(:,:) = H(k,:,:);
    
    % Precoding matrix ���
    G = t_H' * inv( t_H * t_H' );
    
    % Precoding ����
    sym_hat(:,k) = G * sym(:,k);
    W(k,:,:) = G;
    K(k) = sqrt( trace(G * G') );
end




