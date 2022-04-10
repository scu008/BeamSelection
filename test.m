clear, clc, close all
% 특정 벡터와 유사한 공간에 대한 packing 기능을 테스트

% 타겟 벡터를 이용하여 코드북 생성
a = randi([1 10], 3, 1);
[cb, s_sample, sample] = gen_cb(a, 2, pi/2, 1000);

% 타겟 벡터의 주변 샘플 및 코드북 벡터 plot
plot3(sample(1,:), sample(2,:), sample(3,:), '.g');
hold on
plot3(s_sample(1,:), s_sample(2,:), s_sample(3,:), '.b');
plot3(cb(1,:), cb(2,:), cb(3,:), '.r', 'MarkerSize', 20);

% 원점으로 부터 타겟 방향으로의 벡터
a = [zeros(3,1) a];
plot3(a(1,:), a(2,:), a(3,:), '*-m','MarkerSize', 20);