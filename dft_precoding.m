function W = dft_precoding(t_H, N_rf)
% DFT ����� ���� ���� ����� ä���� ���� �����ڵ� ����� ���
% t_H(�ʼ�): 2���� ä�� ���(���ž��׳�, �۽ž��׳� �Ǵ� 3���� ä�� ���(�ιݼ���, ���ž��׳�, �۽ž��׳�)
% N_rf(�ʼ�): RF chain�� �� (N_rf >= N_s);
% W: �Ƴ��α� �����ڴ�

% �⺻ ���� ���� (dim == 3: selective fading, mode = 2: flat fading)
dim = size(t_H);
if length(dim) == 3
    [fft_len, N_d, N_tx] = size(t_H);
elseif dim == 2
    [N_d, N_tx] = size(t_H);
    fft_len = 1;
end

% DFT ��� ����
f = -0.5 : 1/N_tx : 0.5 - 1/N_tx;
cb = exp( -2i*pi*( (0:N_tx-1).' * f ) );

% DFT ��ķκ��� ���� ����
acc = 0;
for i = 1:fft_len
   if length(dim) == 3, H(1:N_d,:) = t_H(i,:,:);
   elseif length(dim) == 2, H = t_H;
   end
   
   % DFT ����� ���� ����� ����
   acc = acc + abs(H * cb).^2;
end

[~, idx] = max(acc,[],2);
W = cb(:,idx);

end