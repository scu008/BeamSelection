function W = steer_precoding(fc, ant, angle, mode)
% �ۼ��� ������ ���� Beamforming matrix�� ����ϴ� �Լ�
% ant_num: ���׳��� �迭 [ ��������� ���׳� ��, ���������� ���׳� ��, ��������� ���׳� ����/lamda, ���������� ���׳� ����/lamda ]
% rx_angle: ä���� �ۼ��� ���� [theta; phi]
% Nrf: RF chain (���� ��Ʈ��)�� ��
% mode: ������ ����(1: All connected, 2: Sub connected)

% ���� �ʱ�ȭ
if nargin < 4, mode = 1; end
if ant(1)*ant(2) == 1, mode = 2; end
[~, Nrf] = size(angle);

% �⺻ ���� ����
%==========================================================
% �迭 ����
c = 3e8;                        % ���� �ӵ�
lamda = c/fc;                   % ��ȣ�� ����
k = 2*pi / lamda;               % �ļ�
dy = ant(3) * lamda;            % �۽� ���׳� �� �Ÿ� (����)
dz = ant(4) * lamda;            % �۽� ���׳� �� �Ÿ� (����)
N = ant(1) * ant(2);

% ���׳� ��ġ ���
temp1 = repmat(0:ant(1)-1, ant(2), 1);
temp2 = repmat(0:ant(2)-1, 1, ant(1));
ant_mat = [ zeros(N,1) temp1(:)*dy (temp2.')*dz];

% ��ǥ�� ��ȯ �Լ�
trans_f = @(t_theta, t_phi) [ sin(t_theta).*cos(t_phi); sin(t_theta).*sin(t_phi); cos(t_theta) ];


% �� matrix ��� 
if mode == 1
    % All connected array�� ����ϴ� ���
    
    W = zeros(N, Nrf);
    for i = 1:Nrf
        theta = angle(1,i);    phi = angle(2,i);
        W(:,i) = exp( -1j*k* ant_mat * trans_f(theta, phi) );
    end
    
elseif mode == 2
    % Sub connected array�� ����ϴ� ���
    
    W = zeros(N*Nrf, Nrf);
    for i = 1:Nrf
        theta = angle(1,i);    phi = angle(2,i);
        W(1+(i-1)*N : i*N, i) = exp( -1j*k* ant_mat * trans_f(theta, phi) );
    end
end


end