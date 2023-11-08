%% testing the gradient code

% grab the calibration
calibFile_coeffs = 'D:\Utah Field Data\BubbleData\9-25-2022\axis_cal_dltCoefs.csv';
all_coeffs = load(calibFile_coeffs);

% grab some positions 
stmFile = 'D:\Utah Field Data\BubbleData\9-25-2022\bub_7\rays_out_cpp.h5';
frames = 1:50;
[X,T,E] = readSTM(stmFile,frames);
pos = X(500,:);

%% gradients wrt position 

% calculate analytically
C = DLTlin_pos(pos,all_coeffs);

% estimate numerically
C_num = DLTlin_pos_num(pos,all_coeffs);

% plot
figure
subplot(1,2,1)
imagesc(C), colorbar

subplot(1,2,2)
imagesc(C_num), colorbar

%% gradients wrt calibration coefficients
% pick a camera
coeffs = all_coeffs(:,4);

% calculate analytically
[vx,vy] = DLTlin_coeff_vec(pos,coeffs);

% estimate numerically
[vx_num,vy_num] = DLTlin_coeff_vec_num(pos,coeffs);

% plot
figure
subplot(2,1,1)
plot(vx,'o'), hold on
plot(vx_num,'x')

subplot(2,1,2)
plot(vy,'o'), hold on
plot(vy_num,'x')

figure
plot((abs(vx-vx_num)),'x'), hold on
plot((abs(vy-vy_num)),'o')

