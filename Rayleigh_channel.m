% Function for generating Rayleigh channel coeffecients made in 2017 %

function [ch_coef, rx_power] = Rayleigh_channel(channel_shape, distance_rate, sig_power, mode, profile)


% �Ű������� ���� �ʱ�ȭ
if nargin < 5, profile = 0; end;
if nargin < 4, mode = 1; end;
if nargin < 3, sig_power = 1; end;
if nargin < 2, distance_rate = 1; end;


% ���� ���� Ȯ��
Multipath = channel_shape(1,1);

[~, c] = size(channel_shape);
if c ~= 1
    row = channel_shape(1,2);
else
    row = 1;
end


% �������� ����
if profile == 0
    % �⺻ ä�� ��������
    CH_Profile = exp( -(1:Multipath) / 5 );
else
    % �Է� ä�� ��������
    CH_Profile = profile;
    [~, Multipath] = size(profile);
end


% �������� ����ȭ �� �ݺ�
CH_Profile = ones(row,1) * CH_Profile / sum(CH_Profile);


% �������ϰ� ���� ����� ��
rnd = ( randn(row,Multipath) + 1j*randn(row,Multipath) ) * sqrt(0.5);


% Racian �� ����
if mode == 2
        
        % LOS power rate
        K = zeros(row,Multipath);
        K(:,1) = 10^(15/10);
        
        % power rate ����
        rnd = sqrt( K./(K+1) ) + rnd .* sqrt( 1./(K+1) );
end
ch_coef = rnd .* sqrt(CH_Profile);


% ����� �Ÿ��� ���� ���� ���� ����=====================================

% �Ϲ��� 3 ~ 4, �̻��� 2
exp_beta = 3;

% �Ϲ��� 5 ~ 12, LTE: 10
sigma = 0;
shadowing = randn(1) * sigma;

% + log10(k) : �ۼ��� ���׳� gain �� ���׳��� ���� �ս� = 1 �� ����ȭ
path_loss = -10 * exp_beta * log10(distance_rate);

% path loss ����
loss = 10^( (path_loss + shadowing) / 10);

% �Ÿ��� ���� ���� ���� ����
ch_coef = ch_coef * sqrt(loss) * sqrt(sig_power);

% ����ȭ �� ���� ����
rx_power = loss * sig_power;


    
end



