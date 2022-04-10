% 3GPP TR 38.900 based three-dimensional Spatial Channel Model(SCM) model
% Version 0.8 made in 2020-05-10

classdef SCM < handle
    
    % ä�θ� ȯ�� ����
    properties
        % Small scale ����
        fc              % ��ȣ�� �߽� ���ļ� (Center Frequency), default: 800MHz
        lamda           % ��ȣ�� ����
        fs              % ä���� sampling frequency
        Ts              % ä���� sampling period
        tx_ant          % �۽Ŵ��� ���׳� �迭 ���� [ row  col  row_dis  col_dis ]
        rx_ant          % ���Ŵ��� ���׳� �迭 ���� [ row  col  row_dis  col_dis ]
        tx_d            % �۽Ŵ��� ���׳� ��ġ ���
        rx_d            % ���Ŵ��� ���׳� ��ġ ���
        Ntx             % �۽Ŵ��� ���׳� ��
        Nrx             % ���Ŵ��� ���׳� ��
        n_path          % ä���� path(cluster)�� ��
        n_mray          % ä���� path�� ray�� �ִ� ����
        n_ray           % ä���� path�� ray�� ��
        asd             % ASD ��
        zsd             % ZSD ��
        asa             % ASA ��
        zsa             % ZSA ��
        xpr_mu          % ���� ���ϰ� ���� ���� ��� 10^(x/10), x-normal dist
        xpr_std         % ���� ���ϰ� ���� ���� ǥ������ 10^(x/10), x-normal dist
        pdp             % Power Delay Profile (PDP)
        
        % Large scale ����
        Gt              % �۽Ŵ��� ���׳� Gain [dB]
        Gr              % ���Ŵ��� ���׳� Gain [dB]
        L               % �ý��� loss factor [dB]
        distance_rate   % �Ÿ� 1�� ä�ο� ���� ������� �Ÿ� ����
        exp_beta        % Path loss exponent
        sdw_std         % Shadowing ����� ǥ������, default: 0
        los             % LOS(Line of Sight) ȯ���� �ݿ����� �����ϴ� ���� 1: LOS, 0: non-LOS
        K               % LOS�� non-LOS ��ȣ ���� ������ ���� [dB]
        No              % AWGN�� PSD, default: -174 [dB/Hz] (-204 [dbm/Hz])
        ZoD_L           % LOS ������ ZoD
        AoD_L           % LOS ������ AoD
        ZoA_L           % LOS ������ ZoA
        AoA_L           % LOS ������ AoA
        p_src           % �۽Ŵ��� 3���� ��ġ
        p_dst           % ���Ŵ��� 3���� ��ġ
        abr_src         % �۽Ŵ��� 3���� ���׳� ���� ���� (alpha, beta, gamma)
        abr_dst         % ���Ŵ��� 3���� ���׳� ���� ���� (alpha, beta, gamma)
        los_flag        % LOS ȯ�� check�� ���� ���� ����
        
        % ��� ���� �Լ� �� ����
        tx_theta        % �۽Ŵ��� ���� �������
        tx_phi          % �۽Ŵ��� ���� �������
        rx_theta        % ���Ŵ��� ���� �������
        rx_phi          % ���Ŵ��� ���� �������
        
        % ��Ÿ �Լ� �� �� ����
        cvt_S2R         % ������ǥ�� ������ǥ�� ��ȯ�ϴ� �Լ�
        R_mat           % GCS�� LCS�� ��ȯ�� ���� ��� ��� �Լ�
        model_var       % �� �Լ� ������ ���� ���� ����
    end
    
    
    % ä�θ� �Լ�
    methods
        % ������ ===========================================================
        function obj = SCM()
            
            % 3���� ��ġ �� �ۼ��� ���� ���� �ʱ�ȭ
            obj.p_src = [0 0 0];
            obj.p_dst = [1 0 0];
            obj.abr_src = [0 0 0];
            obj.abr_dst = [pi 0 0];
            
            % Small scale �ʱⰪ ����
            obj.fc = 800e6;
            obj.lamda = [];
            obj.fs = 20e6;
            obj.Ts = [];
            obj.tx_ant = [1 1 0.5 0.5];
            obj.rx_ant = [1 1 0.5 0.5];
            obj.Ntx = [];
            obj.Nrx = [];
            obj.n_path = 7;
            obj.n_mray = 15;
            obj.n_ray = [];
            obj.asd = 3;
            obj.zsd = 3;
            obj.asa = 3;
            obj.zsa = 3;
            obj.xpr_mu = 8;
            obj.xpr_std = 3;
            obj.pdp = [];
            
            % Large scale �ʱⰪ ����
            obj.Gt = 0;
            obj.Gr = 0;
            obj.L = 0;
            obj.distance_rate = 1;
            obj.exp_beta = 3;
            obj.sdw_std = 0;
            obj.los = 0;
            obj.los_flag = 0;
            obj.K = 15;
            obj.No = -174;
            obj.ZoD_L = pi/2;
            obj.AoD_L = 0;
            obj.ZoA_L = pi/2;
            obj.AoA_L = 0;
            
            % ��� ���� �ʱⰪ ����
            obj.tx_theta = @(phi, theta) 1;
            obj.tx_phi = @(phi, theta) 0;
            obj.rx_theta = @(phi, theta) 1;
            obj.rx_phi = @(phi, theta) 0;
            
            % �Լ� �ʱ�ȭ
            obj.cvt_S2R = @(theta, phi) [sin(theta).*cos(phi); sin(theta).*sin(phi); cos(theta)];
            obj.R_mat = @(alpha, beta, gamma) ...
                [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1] * ...
                [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)] * ...
                [1 0 0; 0 cos(gamma) -sin(gamma); 0 sin(gamma) cos(gamma) ];
        end
        
        
        % �ۼ��� ��ġ�� �ʱ�ȭ �ϴ� �Լ�
        function [res_ang, los_flag] = position(obj, p_src, p_dst, abr_src, abr_dst)
            
            % �ۼ��� ���� ���� �ʱ�ȭ(����ǥ�� ����)
            if nargin < 5
                % �۽ű� ���׳� ���� ���� �ʱ�ȭ
                R = obj.R_mat(0, 0, 0);
                src_tmp = p_dst - p_src;
                src_tmp = src_tmp / norm(src_tmp);
                obj.abr_src = [angle([1 1j 0] * R.' * src_tmp.') -( pi/2 - acos([0 0 1] * R.' * src_tmp.') ) 0];
                
                % ���ű� ���׳� ���� ���� �ʱ�ȭ
                dst_tmp = p_src - p_dst;
                dst_tmp = dst_tmp / norm(dst_tmp);
                obj.abr_dst = [angle([1 1j 0] * R.' * dst_tmp.') -( pi/2 - acos([0 0 1] * R.' * dst_tmp.') ) 0];
            else
                obj.abr_src = abr_src;
                obj.abr_dst = abr_dst;
            end
            
            % 3���� ��ġ ���� �ʱ�ȭ(������ǥ�� ����)
            obj.p_src = p_src;
            obj.p_dst = p_dst;
            
            % �۽� LOS ���� ���
            src_tmp = obj.p_dst - obj.p_src;
            src_tmp = src_tmp / norm(src_tmp);
            R = obj.R_mat(obj.abr_src(1), obj.abr_src(2), obj.abr_src(3));
            obj.ZoD_L = acos([0 0 1] * R.' * src_tmp.');
            obj.AoD_L = angle([1 1j 0] * R.' * src_tmp.');
            
            % ���� LOS ���� ���
            dst_tmp = obj.p_src - obj.p_dst;
            dst_tmp = dst_tmp / norm(dst_tmp);
            R = obj.R_mat(obj.abr_dst(1), obj.abr_dst(2), obj.abr_dst(3));
            obj.ZoA_L = acos([0 0 1] * R.' * dst_tmp.');
            obj.AoA_L = angle([1 1j 0] * R.' * dst_tmp.');
            
            % LOS�� ������ ȯ������ ����
            flag = zeros(1,4);
            if (obj.ZoD_L >= 0) && ( abs(obj.ZoD_L) <= pi ), flag(1) = 1; end
            if abs(obj.AoD_L) <= pi/2, flag(2) = 1; end
            if (obj.ZoA_L >= 0) && ( abs(obj.ZoA_L) <= pi ), flag(3) = 1; end
            if abs(obj.AoA_L) <= pi/2, flag(4) = 1; end
            if sum(flag) == 4, obj.los_flag = 1; end
            
            % ��� ��ȯ
            los_flag = obj.los_flag;
            res_ang = [obj.ZoD_L obj.AoD_L obj.ZoA_L obj.AoA_L];
        end
        
        
        % ���簢 �迭 ���׳��� ������� ���׳� ��ġ ����� �ʱ�ȭ�ϴ� �Լ� ====
        function [] = init_d(obj)
            
            % �۽Ŵ��� ���׳� ��ġ ���
            obj.Ntx = obj.tx_ant(1) * obj.tx_ant(2);
            tdy = obj.tx_ant(3) * obj.lamda;    
            tdz = obj.tx_ant(4) * obj.lamda;
            temp1 = repmat(0:obj.tx_ant(1)-1, obj.tx_ant(2), 1);
            temp2 = repmat(0:obj.tx_ant(2)-1, 1, obj.tx_ant(1));
            obj.tx_d = [ zeros(obj.Ntx,1) temp1(:)*tdy (temp2.')*tdz];
            
            % ���Ŵ��� ���׳� ��ġ ���
            obj.Nrx = obj.rx_ant(1) * obj.rx_ant(2);
            rdy = obj.rx_ant(3) * obj.lamda;
            rdz = obj.rx_ant(4) * obj.lamda;
            temp3 = repmat(0:obj.rx_ant(1)-1, obj.rx_ant(2), 1);
            temp4 = repmat(0:obj.rx_ant(2)-1, 1, obj.rx_ant(1));
            obj.rx_d = [ zeros(obj.Nrx,1) temp3(:)*rdy (temp4.')*rdz];
        end
        
        
        % �۽Ŵܰ� ���Ŵ��� ���׳� ���� �ʱ�ȭ �ϴ� �Լ� =====================
        function [Rx_ant, Tx_ant] = ant(obj, N_rx, N_tx)
            % �Է����� ���� �ۼ��� ���� ���׳� ���� 1���� �迭 ���׳��� ����
            Rx_ant = [N_rx 1 0.5 0.5];
            Tx_ant = [N_tx 1 0.5 0.5];
            
            % ȯ�溯���� ���� ����
            obj.rx_ant = Rx_ant;
            obj.tx_ant = Tx_ant;
            obj.Nrx = Rx_ant(1) * Rx_ant(2);
            obj.Ntx = Tx_ant(1) * Tx_ant(2);
        end
        
        
        % Cluster �� ray�� �Ҵ�Ǵ� ��������� ����ϴ� �Լ� =================
        function [] = def_pow(obj)
            
            % obj.pdp ���� ���� Ȯ��
            if isempty(obj.pdp) == 1
                
                % ���������� ������� �� cluster �� ��� ���� �Ҵ�
                pw = exp( -(1:obj.n_path) / 5 );
                obj.pdp = ( pw / sum(pw) ).';
                
            else
                obj.n_path = length(obj.pdp);
                obj.pdp = ( obj.pdp / sum(obj.pdp) ).';
            end
            
            % obj.n_ray ���� ���� Ȯ��
            if isempty(obj.n_ray) == 1
                obj.n_ray = ones(1,obj.n_path) * obj.n_mray;
            end
        end
        
        
        % Cluster �� ray�� ZoD, AoD, ZoA, AoA�� ����ϴ� �Լ� ==============
        function [res_angle, angle] = gen_angle(obj)
            % c_ang: �ۼ��� ����, c_ang = [ AoD(1:n_path); ZoD(1:n_path); AoA(1:n_path); ZoA(1:n_path); ]
            % res_ang: ��ü �ۼ��� ����
    
            % �� cluster�� ZoD, AoD, ZoA, AoA �߽� �� ����
            angle(1,:) = rand(1, obj.n_path)*pi;         % ZoD
            angle(2,:) = -pi/2 + rand(1, obj.n_path)*pi; % AoD
            angle(3,:) = rand(1, obj.n_path)*pi;         % ZoA
            angle(4,:) = -pi/2 + rand(1, obj.n_path)*pi; % AoA
            
            % �� ray�� ZoD, AoD, ZoA, AoA ����
            res_angle = cell(1, obj.n_path);
            for i = 1 : obj.n_path
                
                % ray�� ���� 1�� ��쿡�� �߽� ������ �״�� �̿�
                if obj.n_ray == 1, tmp_angle = angle;
                else
                    tmp_angle = randn(4, obj.n_ray(i));
                    tmp_angle(1,:) = tmp_angle(1,:) * (obj.zsd * pi/180) + angle(1,i);
                    tmp_angle(2,:) = tmp_angle(2,:) * (obj.asd * pi/180) + angle(2,i);
                    tmp_angle(3,:) = tmp_angle(3,:) * (obj.zsa * pi/180) + angle(3,i);
                    tmp_angle(4,:) = tmp_angle(4,:) * (obj.asa * pi/180) + angle(4,i);
                end
                
                res_angle{i} = tmp_angle;
            end
        end
        
        
        % PAS�� ���� ray�� ��������� ����ϴ� �Լ� =========================
        function pw = pas(obj, res_angle, ray_num)
            % res_angle: ���� cluster�� ���ϴ� ray�� ����
            % ray_num: ��������� ����ؾ��ϴ� ray �� ���� ray�� ��ȣ
            
            % ���� �ʱ�ȭ
            [~, lay_len] = size(res_angle);
            
            % default ������� �Ҵ�
            pw = 1 / sqrt(lay_len);
        end
        
        
        % Ray �� ä�� ����� ����ϴ� �Լ� ===================================
        function [subpath_coeff] = ray_cal(obj, sample_len, ZoD, AoD, ZoA, AoA, xpr, vel)
            
            % ���� ������� ���� ���
            trx_coef = [ obj.rx_theta(ZoA); obj.rx_phi(AoA) ].';
            if xpr == 0, trx_coef = trx_coef * ( exp(2j*pi*rand(1)) .* [1 0; 0 -1] );
            else, trx_coef = trx_coef * ( exp(2j*pi*rand(2)) .* [1 1/sqrt(xpr); 1/sqrt(xpr) 1] ); end
            trx_coef = trx_coef * [ obj.tx_theta(ZoD); obj.tx_phi(AoD)];
            
            % �ۼ��� ���׳� ���� ���� ���
            rx_r = obj.cvt_S2R(ZoA, AoA);
            sub_rx = exp(2j*pi * obj.rx_d * rx_r / obj.lamda);
            tx_r = obj.cvt_S2R(ZoD, AoD);
            sub_tx = exp(2j*pi * obj.tx_d * tx_r / obj.lamda);
            trx_tmp(1,:,:) = trx_coef * sub_rx * sub_tx.';
            
            % ���÷� ���� ���
            dop_tmp = zeros(sample_len,1,2);
            if vel == 0, dop_tmp(:,1,1) = ones(sample_len,1,1);
            else
                t_sample = 0 : obj.Ts : obj.Ts * (sample_len-1);
                dop_tmp(:,1,1) = exp(2j*pi * vel * rx_r / obj.lamda * t_sample);
            end
            
            % subpath �� ����
            subpath_coeff = repmat(trx_tmp,sample_len,1,1) .* repmat(dop_tmp(:,1,1), 1, obj.Nrx, obj.Ntx);
        end
        
        
        % Cluter �� ä�� ����� ����ϴ� �Լ� ==================================
        function [r_coeff, c_ang, res_ang] = FD_channel(obj, sample_len, i_vel)
            % sample_len: �ð� ���� ä�� ���� (�۽� ��ȣ�� ���� ���̿� ����)
            % i_vel: �� ���ÿ� ���� �ӵ� ����(3����) e.g. [160 0 0]: x �������� 160km/h
            % c_ang: �ۼ��� ���� c_ang = [ AoD(1:n_path); ZoD(1:n_path); AoA(1:n_path); ZoA(1:n_path); ]
            % ang: c_ang�� ���� subcluster ����
            
            % ���� �ʱ�ȭ(�ӵ�, ����, ���� ����, ���׳� ���, ù ��° ����� �ε���)
            if nargin < 3, vel = 0; else, vel = i_vel * 5/18; end
            if length(vel) < 3, vel = [vel 0 0]; end
            obj.lamda = (3e8) / obj.fc;
            obj.Ts = 1/obj.fs;
            obj.init_d();
            f_idx = 0;
            
            % �� path�� ��� ���� �� ������ ���
            obj.def_pow();
            [res_ang, c_ang] = obj.gen_angle();
            
            % ������� ���� ���� ���
            xpr = 10.^( ( randn(obj.n_path, obj.n_mray) * obj.xpr_std + obj.xpr_mu ) / 10 );
            
            % �� clusster�� ä�� ��� ���
            coeff = zeros(obj.n_path+1, sample_len, obj.Nrx, obj.Ntx);
            for i = 1:obj.n_path
                
                % 0 ��������� �Ҵ�� ��θ� ��� ���� �� ù ��� index ����
                if obj.pdp(i) == 0, continue; 
                elseif f_idx == 0, f_idx = i;   end
                
                % �� ray�� ä�� ��� ���
                tmp_coeff = zeros(sample_len, obj.Nrx, obj.Ntx);
                for j = 1:obj.n_ray(i)
                    ang = res_ang{i};
                    sub_tmp = obj.ray_cal(sample_len, ang(1,j), ang(2,j), ang(3,j), ang(4,j), xpr(i,j), vel);
                    sub_tmp = sub_tmp * obj.pas(ang,j);
                    tmp_coeff = tmp_coeff + sub_tmp;
                end
                
                % �� cluster�� ��� �Ҵ�
                coeff(i,:,:,:) = tmp_coeff * sqrt(obj.pdp(i));
            end
            
            % LOS ��� ����
            if (obj.los & obj.los_flag) == 1
                Kr = 10^(obj.K/10);
                coeff = sqrt( 1 / (Kr + 1) ) * coeff;
                nlos_tmp = zeros(sample_len, obj.Nrx, obj.Ntx);
                nlos_tmp(:,:,:) = coeff(f_idx,:,:,:);
                coeff_los = obj.ray_cal(sample_len, obj.ZoD_L, obj.AoD_L, obj.ZoA_L, obj.AoA_L, 0, vel);
                coeff(f_idx,:,:,:) = nlos_tmp  +  sqrt( Kr / (Kr + 1) ) * coeff_los;
            end
            
            % ���Ÿ��� �ݿ�
            p_loss = -10 * obj.exp_beta * log10(obj.distance_rate);
            shadowing = randn(1) * obj.sdw_std;
            loss = 10^( (p_loss + shadowing) / 10 );
            coeff = coeff * sqrt(loss);
            
            % ��°� ����
            r_coeff = coeff(1:obj.n_path,:,:,:);
        end
        
        
        % ä�� ����� �̿��Ͽ� ���� ��ȣ�� ����ϴ� �Լ� =====================
        function [rx_sig] = FD_fading(obj, sig, coeff)
            
            % �Ű����� �ʱ�ȭ
            [tap_len, sym_len, N_rx, N_tx] = size(coeff);
            
            % �۽� ��ȣ�� ä�� ��� ����
            temp = zeros(tap_len, sym_len, N_rx);
            for i = 1:N_rx
                for j = 1:N_tx
                    % ���� ���׳����� ä�ΰ� �۽� ��ȣ ����
                    temp(:,:,i) = temp(:,:,i) + coeff(:,:,i,j) .* ( ones(tap_len,1) * sig(j,:) ) ;
                end
            end
            
            % ���� ��� ����
            rx_sig = zeros(N_rx, sym_len + tap_len -1);
            for m = 1:N_rx
                % ���� ��� �迭 ����
                flat_m = zeros(tap_len, sym_len + tap_len -1);
                
                for i = 1:tap_len
                    flat_m(i, i:i + sym_len - 1) = temp(i,:,m);
                end
                
                % ���� ��� ��ø
                rx_sig(m,:) = sum( flat_m, 1 );
            end
        end
        
        % �־��� �����ð� �� ��������� ���� �ý��ۿ� �´� PDP�� ����ϴ� �Լ�
        function [] = pdp_interp(obj, delay, power)
            % Delay [sec]
            % power [Watt]
            
            % ���̰� �ٸ� ��� ����
            if length(delay) ~= length(power)
                disp( 'The sizes of delay and power doesn not match ')
                return
            end
            
            % PDP ���� ����
            ts = 1/obj.fs;
            
            % index ���� ��ġ �� ���� ������� ���� �л����� ���� ���� ���� '2'
            len = ceil( max(delay) / ts ) + 2;
            obj.pdp = zeros(1, len);
            
            % �ø��Ͽ� ��ġ ���
            coef = floor( delay / ts);
            relative_d = abs(delay / ts - coef);
            p0 = (1-relative_d) .* power;
            p1 = relative_d .* power;
            
            % Power ���� ����Ͽ� ����
            for i = 1:length(delay)
                idx = coef(i)+1;
                obj.pdp( idx ) = obj.pdp( idx ) + p0(i);
                obj.pdp( idx+1) = obj.pdp( idx+1 ) + p1(i) ;
            end
            
            % ����ȭ ����
            obj.pdp = obj.pdp / sum(obj.pdp);
            
        end
        
        % �Ÿ� �� �۽� ������ ������� ��ȣ�� ���� SNR�� ����ϴ� �Լ� ========
        function [snr] = path_loss(obj, tx_psd, bandwidth, distance)
            % tx_psd: ���� ��ȣ�� ��� PSD(Power Spectral Density) [dB/Hz]
            % bandwith: ���� ��ȣ�� �뿪��
            
            % ���� �ʱ�ȭ
            if nargin < 4
                distance  =  sqrt( sum( (obj.p_src - obj.p_dst).^2 ) );
            end
            
            % �۽�, ���� ���� �� path loss ���
            tx_power = 10^(tx_psd/10) * bandwidth;
            N_power = 10^(obj.No/10) * bandwidth;
            p_loss = -20 .* log10( 4*pi / obj.lamda ) - 10 .* obj.exp_beta .* log10( distance ) - obj.L + obj.Gt + obj.Gr;
            
            % Shadowing ������ �ݿ��Ͽ� SNR�� ���
            shadowing = randn(1) * obj.sdw_std;
            loss = 10^( (p_loss + shadowing) / 10 );
            rx_power = tx_power * loss;
            snr = 10*log10(rx_power / N_power);
        end
        
        
    end
end