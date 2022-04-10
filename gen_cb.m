function [cb, sample, i_sample] = gen_cb(a, n_bit, ang, num)
% 주어진 A 행렬을 이용하여 Grassman packing을 통해 codebook을 계산하는 함수
% 채널의 계수가 Rayleigh 분포를 따른다고 가정
% a: 주어진 특정 기저 벡터
% n_bit: 생성할 코드북의 크기를 결정하는 bit 수
% ang: 양자화에 고려되는 타겟 벡터와의 이격 각도
% num: 양자화를 위해 고려되는 샘플의 수(코드북 크기 * num);

% 각 행벡터의 크기를 1로 정규화 하는 함수
nom = @(x) x ./ sqrt(sum(abs(x).^2, 2));

% 변수 초기화
if nargin < 4, num = 1500000; end
if nargin < 3, ang = 2*pi; end
iter  = 1000;
[dim, ~] = size(a);
if sum(imag(a)) == 0, mode = 'r'; else mode = 'c'; end
M = 2^n_bit;    % 양자화 상태의 수
Q = num;    % 코드북이 대표하는 샘플들의 수
if mode == 'r', i_sample = nom( randn(Q,dim) );
elseif mode == 'c', i_sample = nom( 0.5 * (randn(Q,dim) + 1j*randn(Q,dim)) );
end

% 주어진 기저벡터와 유사도가 높은 벡터를 먼저 선별 후 코드북 후보 초기화
sample = i_sample;
idx = find( ang < acos( real( sample * nom(a.')' ) ) );
sample(idx,:) = [];
[Q, ~] = size(sample);
cb = sample(1:M,:);

% 최적화 루틴
dis = zeros(M,Q);
group = zeros(1,Q);
flag = 0;
for i = 1:iter
    % 가까운 샘플들을 그룹화
    for j = 1:M
        dis(j,:) = sum( abs(ones(Q,1)*cb(j,:) - sample).^2, 2 ).';
    end
    [~, n_group] = min(dis);
    
    % 코드북의 벡터들을 그룹의 평균 위치로 이동
    for j = 1:M
        idx = find(group == j);
        cb(j,:) = nom( mean( sample(idx,:) ) );
    end
    
    % 최적화 중지 조건을 검토
    if sum(group == n_group), flag = flag + 1;end
    if flag > 4, break; end
    
    % 변수 재할당
    group = n_group;
end

% 결과 행렬 전치
i_sample = i_sample.';
sample = sample.';
cb = cb.';

end