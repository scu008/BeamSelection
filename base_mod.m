function [mod_data] = base_mod(data,mod_scheme)

if mod_scheme == 1       % BPSK modulation
    mod_data = data*2-1;

elseif mod_scheme == 2   % QPSK modulation
    % 가로가 홀수라면 맨뒤에 0을 채워줌
    [MP_row,MP_col]=size(data);  
    if mod(MP_col,2) ~= 0
        data(:,MP_col+1)=0;
        MP_col=MP_col+1;
    end

    odd_data=data(:,1:2:MP_col-1); % 홀수만 분리
    odd_data=(odd_data*2-1)*j;     
    even_data=data(:,2:2:MP_col);  % 짝수만 분리
    even_data=even_data*2-1;

    mod_data=(odd_data+even_data).*0.7071;
    
elseif mod_scheme == 4   % 16QAM modulation
    % 가로가 4의 배수가 아니라면 뒤에 0을 채워 4의 배수로 만들어줌
    [MP_row,MP_col]=size(data);
    while mod(MP_col,4) ~= 0
        data(:,MP_col+1) = 0;
        MP_col = MP_col+1;
    end

    first_data = data(:,1:4:MP_col-3);   % 1st
    first = first_data*4-2;     
    second_data = data(:,2:4:MP_col-2);  % 2nd
    second = (second_data~=first_data)*2-1;
    third_data = data(:,3:4:MP_col-1);   % 3rd
    third = (third_data*4-2)*j;     
    fourth_data = data(:,4:4:MP_col);    % 4th
    fourth = ((fourth_data~=third_data)*2-1)*j;
    
    mod_data=(first+second+third+fourth).*0.3162;
  
elseif mod_scheme == 6   % 64QAM modulation
    % 가로가 6의 배수가 아니라면 뒤에 0을 채워 6의 배수로 만들어줌
    [MP_row,MP_col]=size(data);
    while mod(MP_col,6) ~= 0
        data(:,MP_col+1) = 0;
        MP_col = MP_col+1;
    end

    first_data=data(:,1:6:MP_col-5);   % 1st
    first=first_data*8-4;     
    second_data=data(:,2:6:MP_col-4);  % 2nd
    second=(second_data~=first_data)*4-2;
    third_data=data(:,3:6:MP_col-3);   % 3rd
    third=(xor(xor(first_data,second_data),third_data))*2-1;
    fourth_data=data(:,4:6:MP_col-2);  % 4th
    fourth=(fourth_data*8-4)*j;
    fifth_data=data(:,5:6:MP_col-1);   % 5th
    fifth=((fourth_data~=fifth_data)*4-2)*j;     
    sixth_data=data(:,6:6:MP_col);     % 6th
    sixth=((xor(xor(fourth_data,fifth_data),sixth_data))*2-1)*j;

    mod_data=(first+second+third+fourth+fifth+sixth).*0.1543;
    
elseif mod_scheme == 8   % 256QAM modulation
    % 가로가 8의 배수가 아니라면 뒤에 0을 채워 8의 배수로 만들어줌
    [MP_row,MP_col]=size(data);
    while mod(MP_col,8) ~= 0
        data(:,MP_col+1) = 0;
        MP_col = MP_col+1;
    end

    first_data=data(:,1:8:MP_col-7);   % 1st
    first=first_data*16-8;     
    second_data=data(:,2:8:MP_col-6);  % 2nd
    second=(second_data~=first_data)*8-4;
    third_data=data(:,3:8:MP_col-5);   % 3rd
    third=(xor(xor(first_data,second_data),third_data))*4-2;
    fourth_data=data(:,4:8:MP_col-4);  % 4th
    fourth=(xor(xor(xor(first_data,second_data),third_data),fourth_data))*2-1;
    fifth_data=data(:,5:8:MP_col-3);   % 5th
    fifth=(fifth_data*16-8)*j;     
    sixth_data=data(:,6:8:MP_col-2);   % 6th
    sixth=((sixth_data~=fifth_data)*8-4)*j;
    seventh_data=data(:,7:8:MP_col-1);   % 7th
    seventh=((xor(xor(fifth_data,sixth_data),seventh_data))*4-2)*j;
    eighth_data=data(:,8:8:MP_col);     % 8th
    eighth=((xor(xor(xor(fifth_data,sixth_data),seventh_data),eighth_data))*2-1)*j;
    
    mod_data=(first+second+third+fourth+fifth+sixth+seventh+eighth).*0.0767;
        
else
    disp('매개변수를 확인해 주세요');
    mod_data=0;
    
end