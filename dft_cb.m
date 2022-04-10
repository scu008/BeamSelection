function [cb] = dft_cb(N_tx, num)
% DFT 기반의 프리코딩 코드북 matrix를 생성
% N_tx: 송신을 위한 안테나 수
% num: Full 코드북에서 등 간격으로 선택되는 submatrix
% (num 미입력 시 full matrix return)

% 변수 초기화
if nargin < 2, num = N_tx; end

% exp 매개변수 생성
f = -0.5 : 1/num : 0.5 - 1/num;
cb = exp( -2i*pi*( (0:N_tx-1).' * f ) );