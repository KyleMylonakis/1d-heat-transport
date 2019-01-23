load('soe_weights_and_nodes.mat');
gamma_real = real(gamma);
gamma_imag = imag(gamma);
R_real = real(R);
R_imag = imag(R);



gamma_real_fileID = fopen('gamma_real.bin', 'w');
fwrite(gamma_real_fileID,gamma_real,'double');
fclose(gamma_real_fileID);

gamma_imag_fileID = fopen('gamma_imag.bin', 'w');
fwrite(gamma_imag_fileID,gamma_imag,'double');
fclose(gamma_imag_fileID);

R_real_fileID = fopen('R_real.bin', 'w');
fwrite(R_real_fileID,R_real,'double');
fclose(R_real_fileID);

R_imag_fileID = fopen('R_imag.bin', 'w');
fwrite(R_imag_fileID,R_imag,'double');
fclose(R_imag_fileID);