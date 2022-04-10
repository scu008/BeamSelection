% 채널 피드백을 위한 코드북 생성 (선형 코드북: c_cscb)
c_cscb = [];
name = "c_cscb_tx" + num2str(Ntx) + "_h" + num2str(h_bit) + "_cs" + num2str(cs_dim);
if exist(name + ".mat") == 2, load(name + ".mat")
else
    c_cscb(:,:) = gen_cscb(1j*ones(cs_dim,1), h_bit);
    save(name, 'c_cscb')
end

% 채널 피드백을 위한 코드북 생성 (비선형 코드북: p_cscb)
p_cscb = [];
name = "p_cscb_tx" + num2str(Ntx) + "_h" + num2str(h_bit) + "_cs" + num2str(cs_dim);
if exist(name + ".mat") == 2, load(name + ".mat")
else
    p_cscb(:,:) = gen_cscb(1j*ones(cs_dim,1), h_bit, 0);
    save(name, 'p_cscb')
end