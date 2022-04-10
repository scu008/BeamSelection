% Beam sweeping 과 quantized feedback 기반의 hybrid beamforming 기법

% 전송 데이터 스펙 설정
model = SCM();
model.sdw_std = sdw_sig;
model.n_mray = scatter;
model.n_path = path;
cp_len = fft_len / 4;
data_len = fft_len * mod_type;

% 송신단 및 수신단에서의 총 데이터 스트림 설정
N_s = sel_node;
N_d = rx_node;
N_rf = num_rf;

% 송수신 각 배열 안테나의 안테나 수
model.ant(rx_ant, tx_ant);
N_tx = model.Ntx;
N_rx = model.Nrx;

% 채널 관련 파라미터 초기화
H = zeros(path, 2*(fft_len+cp_len), N_rx * N_s, N_tx);
t_He = zeros(path, N_s, N_rf);
He = zeros(fft_len, N_s, N_rf);

% Beam sweep을 위한 codebook 생성
b_cb = dft_cb(N_tx);
b_cb = b_cb ./ sqrt(sum(abs(b_cb).^2,1));
r_cb = dft_cb(N_rx);


% 반복 측정 시작
temp_num = iter;
result_ber = zeros(1, length(snr) );
result_cap = zeros(1, length(snr) );
for i = 1:iter
    
%         tic
    for j = 1:length(snr)
        
        %============= 채널 계산 ============================================
        
        % 채널 계수 생성
        max_ang = zeros(1,N_d);
        est_ang = zeros(1,N_d);
        for d = 1:N_d
            [temp, rx_angle] = model.FD_channel((fft_len + cp_len)*2);
            H(1:path,:,1+(d-1)*N_rx:d*N_rx,:) = temp;
        end
        
        % 시간 채널 주파수 변환
        tmp_H(1:path, 1:N_rx*N_d, :) = H(:,1,:,:);
        H_f = fft(tmp_H, fft_len, 1);
        
        % Beam sweeping
        rx_tmp = zeros(fft_len, N_d*N_rx, 2^b_bit);
        for k = 1:fft_len
            H_f_(1:N_rx*N_d,:) = H_f(k,:,:);
            [rx_tmp(k,:,:), No] = awgn_noise( H_f_ * b_cb, snr(j) );
        end
        
        % Rx beam sweeping
        Wr = zeros(N_rx*N_d,N_d);
        if N_rx > 1
            for d = 1:N_d
                acc = 0;
                rx_tmp_ = zeros(fft_len, N_rx, N_tx);
                for k = 1:fft_len
                    H_f__(1:N_rx,:) = rx_tmp(k,1+(d-1)*N_rx:d*N_rx,:);
                    rx_tmp_(k,:,:) = r_cb.' * H_f__;
                    acc = acc + sum( abs(r_cb.' * H_f__).^2, 2);
                end
                [~, idx] = max(acc);
                Wr(1+(d-1)*N_rx:d*N_rx,d) = r_cb(:,idx);
                rx_H(:,d,:) = rx_tmp_(:,idx,:);
            end
        else
            rx_H = rx_tmp;
            Wr = eye(N_d);
        end
        
        % Channel feedback and estimation
        H_hat = zeros(fft_len,N_d,N_tx);
        acc_coef = zeros(N_d, 2^b_bit);
        for d = 1:N_d
            [H_hat(:, d, :), ~, ~, h_coef] = cscb_est(rx_H(:,d,:), b_cb, h_cb, 'n', No);
            acc_coef(d,:) = sum(abs(h_coef).^2, 1);
        end
        
        % Beamforming 계수 계산 (송신, 수신 계수)
        if N_s < N_d
            comb_idx = nchoosek(1:N_d, N_s).';
            [~, comb_num] = size(comb_idx);
            sum_coef = zeros(1,comb_num);
            
            % 수신기 조합에 따른 이득을 계산
            for k = 1:comb_num
                [~, sum_coef(k)] = sparse_precoding(acc_coef(comb_idx(:,k),:), N_rf, dft_cb(N_tx));
            end
            
            % 최적 조합을 결정
            [~, opt_comb] = max(sum_coef);
            opt_idx = comb_idx(:,opt_comb);
            Wt = sparse_precoding(acc_coef(opt_idx,:), N_rf, dft_cb(N_tx));
        else
            Wt = sparse_precoding(acc_coef, N_rf, dft_cb(N_tx));
            opt_idx = (1:N_d).';
        end
        
        % 채널에서 선택된 수신기의 인덱스를 설정
        sel_idx = ( repmat(opt_idx-1,1,N_rx)*N_rx + repmat(1:N_rx,N_s,1) ).';
        sel_idx = sel_idx(:);
        
        %============== 신호 송신 ===========================================
        
        % 전송 심볼 생성
        bit = randi([0 1], N_s, data_len);
        Dsym = base_mod(bit, mod_type);
        p_sym = base_mod(bit, mod_type);
        
        % precoding 수행
        for k = 1:fft_len
            tmp_H_f(:,:) = H_hat(k,opt_idx,:);
            He(k,:,:) = tmp_H_f * Wt;
        end
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
        tx_ofdm = [Wt * tx_ofdm Wt * tx_ofdm];
        
        
        %============== 신호 수신 ===========================================
        % ========================= BER 계산 ==============================
        % 채널 통과
        [rx_ofdm, No] = awgn_noise( model.FD_fading( tx_ofdm, H(:,:,sel_idx,:) ), snr(j) );
        
        % 수신 빔포밍 계수 적용
        rx_ofdm = Wr(sel_idx,opt_idx).' * rx_ofdm;
        
        % OFDM 복조
        rx_Isym = rx_ofdm(:, cp_len + 1 : fft_len + cp_len);
        rx_Dsym = fft(rx_Isym, fft_len, 2) / sqrt(fft_len);
        rx_Dsym = rx_Dsym .* factor;
        
        % 파일럿 복조
        i2 = fft_len + cp_len;
        p_Isym = rx_ofdm(:, i2 + cp_len + 1 : i2 + fft_len + cp_len);
        p_Dsym = fft(p_Isym, fft_len, 2) / sqrt(fft_len);
        p_Dsym = p_Dsym .* factor;
        p_h = p_Dsym ./ p_sym;
        rx_Dsym = rx_Dsym ./ p_h;
        
        
        % base demodulation
        rx_bit = base_demod(rx_Dsym, mod_type);
        result_ber(j) = result_ber(j) + sum( sum( bit ~= rx_bit, 2 ), 1 ) / (data_len * N_s);
        
        
        % ======================== Capacity 계산 ==========================
        
        % 각 부반송파에 대한 자기 채널 및 간섭 채널을 계산
        H_A = zeros(N_s,fft_len);
        H_I = zeros(N_s,fft_len);
        tmp_ = fft(H(:,:,sel_idx,:), fft_len, 1);
        for k = 1:fft_len
            t_Wd(:,:) = Wd(k,:,:);
            tmp__(:,1:N_tx) = tmp_(k,1,:,:);
            
            % 자기 채널
            H_A(:,k) = Wr(sel_idx,opt_idx).' * tmp__ * Wt * t_Wd / factor(k) .* eye(N_s) * ones(N_s,1);
            
            % 간섭 채널
            H_I(:,k) = Wr(sel_idx,opt_idx).' * tmp__ * Wt * t_Wd / factor(k) .* (ones(N_s) - eye(N_s)) * ones(N_s,1);
            n = Wr(sel_idx,opt_idx).' * awgn_noise(zeros(N_rx*N_s,1), snr(j));
            H_I(:,k) = H_I(:,k) + n;
        end
        
        % 최종 sum rate 계산
        res_temp = mean( sum( log2( ones(N_s,1) + abs(H_A).^2 ./ abs(H_I).^2 ) ) );
        if isnan(res_temp)
            temp_num = temp_num - 1;
            continue;
        else, result_cap(j) = result_cap(j) + res_temp; end
        
    end
%         toc
    
end


% 결과 평균화
result_ber = result_ber / iter;
result_cap = result_cap / temp_num;

if mode == 1
    % Capacity 결과 출력
    plot(snr, result_cap, plot_format);
    title('Sum Rate Performance')
    legend('Sum Capacity')
    ylabel('Sum Rate (bps/Hz)')
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
    
elseif mode == 3
    % Capacity 결과 출력
    plot(snr, result_cap / N_s, plot_format);
    title('Sum Rate Performance')
    legend('Sum Capacity')
    ylabel('Average Sum Rate (bps/Hz/receiver)')
    xlabel('SNR (dB)')
    grid on
end
