% Hybrid beamforming test (all-array)

% 전송 데이터 스펙 설정
model = SCM();
model.n_mray = scatter;
model.n_path = path;
cp_len = fft_len / 4;
data_len = fft_len * mod_type;

% 송신단 및 수신단에서의 총 데이터 스트림 설정
N_s = rx_node;
N_d = rx_node;
N_rf = rx_node;

% 송수신 각 배열 안테나의 안테나 수 (수신, 송신)
model.ant(rx_ant, tx_ant);
N_tx = model.Ntx;
N_rx = model.Nrx;

% ber 결과 출력 벡터들
result_ber = zeros(1, length(snr) );
result_cap = zeros(1, length(snr) );

% 채널 관련 파라미터 초기화
H = zeros(path, fft_len+cp_len, N_rx * N_d, N_tx);
sel_angle = zeros(4, N_d);
t_He = zeros(path, N_d, N_rf);
He = zeros(fft_len, N_d, N_s);


tic
for i = 1:iter
    for j = 1:length(snr)
        
        %============= 채널 계산 ============================================
        
        % 채널 계수 생성
        for d = 1:N_d
            [temp, rx_angle] = model.FD_channel(fft_len + cp_len);
            H(:,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
            
            % 각도 탐색
            [~, idx] = max( abs( temp(:,1,1,1) ) );
            sel_angle(:,d) = rx_angle(:,idx);
        end
        
        % Beamforming 계수 계산 (송신, 수신 계수)
        Wt = steer_precoding(model.fc, model.tx_ant, sel_angle(1:2,:));
        Wr = steer_precoding(model.fc, model.rx_ant, sel_angle(3:4,:), 2);
        
        % 시간 영역 Effective 채널 계수 계산
        for k = 1:path
            tmp(:,:) = H(k,1,:,:);
            t_He(k,:,:)= Wr.' * tmp * Wt;
        end
        
        
        %============== 신호 송신 ===========================================
        
        % 전송 심볼 생성
        bit = randi([0 1], N_s, data_len);
        Dsym = base_mod(bit, mod_type);
        
        % precoding 수행
        He = fft(t_He, fft_len, 1);
        [Dsym, ~, Wd] = ZF_precoding(Dsym, He);
        
        % 정규화 상수 계산
        factor = zeros(1,fft_len);
        for k = 1:fft_len
            t_Wd(:,:) = Wd(k,:,:);
            temp = Wt * t_Wd;
            factor(k) = sqrt( trace( temp * temp' ) );
        end
        
        % OFDM 심볼 생성
        Dsym = Dsym ./ factor;
        Isym = ifft(Dsym, fft_len, 2) * sqrt(fft_len);
        tx_ofdm = [ Isym(:, fft_len - cp_len + 1 : end) Isym ];
        
        
        % 송신 빔계수 적용
        tx_ofdm = Wt * tx_ofdm;
        
        
        %============== 신호 수신 ===========================================
        % ========================= BER 계산 ==============================
        % 채널 통과
        [rx_ofdm, No] = awgn_noise( model.FD_fading( tx_ofdm, H ), snr(j) );
        
        % 수신 빔포밍 계수 적용
        rx_ofdm = Wr.' * rx_ofdm;
        
        % OFDM 복조
        rx_Isym = rx_ofdm(:, cp_len + 1 : fft_len + cp_len);
        rx_Dsym = fft(rx_Isym, fft_len, 2) / sqrt(fft_len);
        rx_Dsym = rx_Dsym .* factor;
        
        % base demodulation
        rx_bit = base_demod(rx_Dsym, mod_type);
        result_ber(j) = result_ber(j) + sum( sum( bit ~= rx_bit, 2 ), 1 ) / (data_len * N_s);
        
        
        % ======================== Capacity 계산 ==========================
                
        % 각 부반송파에 대한 자기 채널 및 간섭 채널을 계산
        H_A = zeros(N_d,fft_len);
        H_I = zeros(N_d,fft_len);
        tmp_ = fft(H, fft_len, 1);
        for k = 1:fft_len
            t_Wd(:,:) = Wd(k,:,:);
            tmp__(:,1:N_tx) = tmp_(k,1,:,:);
            
            % 자기 채널
            H_A(:,k) = Wr.' * tmp__ * Wt * t_Wd / factor(k) .* eye(N_s) * ones(N_s,1);
            
            % 간섭 채널
            H_I(:,k) = Wr.' * tmp__ * Wt * t_Wd / factor(k) .* (ones(N_s) - eye(N_s)) * ones(N_s,1);
            n = Wr.' * awgn_noise(zeros(N_rx*N_d,1), snr(j));
            H_I(:,k) = H_I(:,k) + n;
        end
        
        % 최종 sum rate 계산
        result_cap(j) = result_cap(j) + mean( sum( log2( ones(N_d,1) + abs(H_A).^2 ./ abs(H_I).^2 ) ) );
        
    end
end
toc

% 결과 평균화
result_ber = result_ber / iter;
result_cap = result_cap / iter;


if mode == 1
    % Capacity 결과 출력
    plot(snr, result_cap, plot_format);
    title('Sum Rate Performance')
    legend('Sum Capacity')
    ylabel('Average Spectral Efficiency (bps/Hz)')
    xlabel('SNR (dB)')
    grid on
    
elseif mode == 2
    % BER 결과 출력
    semilogy(snr, result_ber, plot_format);
    title('BER Performance')
    legend('Detection result')
    ylabel('BER')
    xlabel('SNR (dB)')
    grid on
    axis([snr(1) max(snr) 10^-5 1])
end

