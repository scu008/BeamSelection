clear, clc
% 양자화 성능을 테스트
iter = 20000;

% 양자화 bit 수(아날로그 빔, 채널 정보)
b_bit = 3;
h_bit = 4;
cs_dim = 2;

% 채널 모델 변수 설정
fft_len = 64;
Ntx = 2^b_bit;
Nrx = 1;
path = 2;
ray = 10;
model = SCM();
model.n_mray = ray;
model.n_path = path;
model.ant(Nrx, Ntx);
snr = 0:5:30;

% DFT 기반의 코드북 생성
b_cb = dft_cb(Ntx);
b_cb = b_cb ./ sqrt(sum(abs(b_cb).^2,1));
r_cb = dft_cb(Nrx);

% 채널 양자화를 위한 코드북 로드
load_codebook

result = zeros(3,length(snr));
for j = 1:iter
%     tic
    
    for s = 1:length(snr)
        
        % 채널 계수 생성 및 FFT
        H_ = model.FD_channel(1);
        H(1:path,1:Nrx,1:Ntx) = H_(:,1,:,:);
        H_f = fft(H, fft_len, 1);
        
        % Channel after beam sweeping
        rx_tmp = zeros(fft_len, Nrx, 2^b_bit);
        for k = 1:fft_len
            H_f_(1:Nrx,:) = H_f(k,:,:);
            [rx_tmp(k,:,:), No] = awgn_noise(H_f_ * b_cb, snr(s));
        end
        
        % Rx beam sweeping
        [~, Mr] = size(r_cb);
        if Nrx > 1
            acc = 0;
            rx_tmp_ = zeros(fft_len, Mr, 2^b_bit);
            for k = 1:fft_len
                H_f__(1:Nrx,:) = rx_tmp(k,:,:);
                rx_tmp_(k,:,:) = r_cb.' * H_f__;
                acc = acc + sum( abs(r_cb.' * H_f__).^2, 2);
            end
            [~, idx] = max(acc);
            rx_H(:,1,:) = rx_tmp_(:,idx,:);
        else
            rx_H = rx_tmp;
        end
        
        
        % 1: CS를 이용한 채널 추정
        H_z = cs_est(rx_H, b_cb, cs_dim);
        
        % 2: 선형 코드북을 통한 채널 정보 선택
        H_c = cscb_est(rx_H, b_cb, c_cscb, 'n', No);
        
        % 3: 아날로그 빔 선택 및 비선형 코드북을 통한 채널 정보 선택
        H_p = cscb_est(rx_H, b_cb, p_cscb, 'n', No);
        
        
        % 두 코드북으로 추정된 채널의 loss를 계산
        z_loss = zeros(fft_len,1);
        c_loss = zeros(fft_len,1);
        p_loss = zeros(fft_len,1);
        for i = 1:fft_len
            
            if Nrx > 1
                H_tmp_(1:Nrx,:) = H_f(i,:,:);
                H_tmp(1,:) = r_cb(:,idx).' * H_tmp_;
            else
                H_tmp(1:Nrx,:) = H_f(i,:,:);
            end
            
            z_tmp(1,:) = H_z(i,:,:);
            c_tmp(1,:) = H_c(i,:,:);
            p_tmp(1,:) = H_p(i,:,:);
            deno = trace(H_tmp * H_tmp');
            tmp_z = abs( H_tmp - z_tmp ).^2;
            tmp_c = abs( H_tmp - c_tmp ).^2;
            tmp_p = abs( H_tmp - p_tmp ).^2;
            z_loss(i) = sum( sum( tmp_z ) ) / deno;
            c_loss(i) = sum( sum( tmp_c ) ) / deno;
            p_loss(i) = sum( sum( tmp_p ) ) / deno;
        end
        
        
        result(1,s) = result(1,s) + sum(z_loss);
        result(2,s) = result(2,s) + sum(c_loss);
        result(3,s) = result(3,s) + sum(p_loss);
        
    end
    
%     toc
    
end

result = result / iter / fft_len;
result = 10*log10(result);
plot(snr.', result.');
legend('CS-OMP', 'QCS-L', 'QCS-NL')
grid on
hold on



