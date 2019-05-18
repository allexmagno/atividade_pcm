% Allex Magno Andrade
% Exercício de Sistemas de Comunicação
% Modulação PCM

% Script desenvolvido utilizando octave carregando os pacotes signal e communications
pkg load signal
pkg load communications

%%
clear all
close all
clc

% Lendo um som
figure(1)
[y Fs] = wavread('you-dont-know-how.wav');
subplot(311);
x_n = y(:, 1:end/2);
plot(x_n)
title('Som original')
%sound(x_n,Fs);
Y_f = fftshift(fft((x_n)/length(x_n)));
f = -Fs/2:1/(length(y)/Fs):Fs/2-(1/(length(y)/Fs));
subplot(312)
plot(f, Y_f)
plot(f, transpose(abs(Y_f)));
title('Componentes Frequenciais')
xlim([-10e3 10e3])  
Vpp = 2*ceil(max([max(x_n) abs(min(x_n))]));

% qtd de bits do quantizador
%k = 3;
%k = 5; 
%k = 8;
k = 13;

% Passo de quantizacao
passo = Vpp/2^k; 

%%
x_desl = x_n + ((Vpp/2) - (passo/2)); % Deslocando o sinal
x_desl2 = x_desl/passo; % Dividindo pelo passo para ter os niveis

x_qtz = round(x_desl2);

%aux_1 =  x_qtz == 2^k; % Caso daria 8
%aux_2 =  x_qtz == -1; % caso surgisse -1
%x_qtz = x_qtz - aux_1 + aux_2; % Para ficar dentro de 0 e 2^k
x_bin = de2bi(x_qtz); % vai gerar uma matriz
[x, y] = size(x_bin);
% Fim da quantizacao

%% Processo de compressao
i = 1;
x_c = [];
while i <= x
  cont = 2;
  flag = 1;
  while cont <= 8
    if x_bin(i,cont:cont) == 1
      if flag == 1
        j = cont;
        flag = 0;
      end
    end
    cont = cont + 1;
  end
  
  if flag == 0
    if j == 2
      s = [1 1 1];
      comp = [x_bin(i,1:1) s x_bin(i,j+1:6)];
    elseif j == 3
      s = [1 1 0];
      comp = [x_bin(i,1:1) s x_bin(i,j+1:7)]; 
    elseif j == 4
      s = [1 0 1];
      comp = [x_bin(i,1:1) s x_bin(i,j+1:8)];
    elseif j == 5
      s = [1 0 0];
      comp = [x_bin(i,1:1) s x_bin(i,j+1:9)];
    elseif j == 6
      s = [0 1 1];
      comp = [x_bin(i,1:1) s x_bin(i,j+1:10)];
    elseif j == 7
      s = [0 1 0];
      comp = [x_bin(i,1:1) s x_bin(i,j+1:11)];
    elseif j == 8
      s = [0 0 1];
      comp = [x_bin(i,1:1) s x_bin(i,j+1:12)];
    end
    else
      s = [0 0 0];
      comp = [x_bin(i,1:1) s x_bin(i,9:12)];
  end
  i = i + 1;
  x_c = [x_c; comp];

end

%x_bin_linha = reshape(transpose(x_bin), 1, x*y); % Para colocar tudo em uma unica linha 




% Observando o sinal transmitido
subplot(313)
filtro_NRZ = ones(1, 100); % filtro
%tx = filter(filtro_NRZ, 1, x_bin_linha);
%tx = filter(filtro_NRZ, 1, upsample(x_bin_linha, 100));
%tx = filter(filtro_NRZ, 1, upsample(x_comp, 100));
%plot(tx);
%ylim([-0.5 2])
%xlim([2.3e6 2.35e6])
%title('Transmissão do sinal')

%% Recuperando o Sinal
% rrr = reshape(xxx, coluna, linha)
% reshape(rrr', linha, coluna)
%x_compressado = reshape(x_comp, 8, length(x_comp)/8);
%x_c = reshape(x_compressado', length(x_comp)/8, 8);

%% expansão
i = 1;
x_exp = [];
[x, y] = size(x_c);
while i <= x
  s = bi2de([x_c(i,2:4)]);
  
  if s == 0
    b = [x_c(i,1:1) 0 0 0 0 0 0 0 x_c(i,5:8) 1];   

  elseif s == 1
    b = [x_c(i,1:1) 0 0 0 0 0 0 1 x_c(i,5:8) 1];
    
  elseif s == 2
    b = [x_c(i,1:1) 0 0 0 0 0 1 x_c(i,5:8) 1 0];
    
  elseif s == 3
    b = [x_c(i,1:1) 0 0 0 0 1 x_c(i,5:8) 1 0 0];
    
  elseif s == 4
    b = [x_c(i,1:1) 0 0 0 1 x_c(i,5:8) 1 0 0 0];
    
  elseif s == 5
    b = [x_c(i,1:1) 0 0 1 x_c(i,5:8) 1 0 0 0 0];
    
  elseif s == 6
    b = [x_c(i,1:1) 0 1 x_c(i,5:8) 1 0 0 0 0 0];
    
  elseif s == 7
    b = [x_c(i,1:1) 1 x_c(i,5:8) 1 0 0 0 0 0 0];
  end
  i = i + 1;
  x_exp = [x_exp; b];
end

%x_exp_dec = reshape(x_exp, k, length(x_exp)/k);
%xdec_qtz = reshape(x_exp_dec', length(x_exp)/k, k); % sinal ainda quantizado
%xdec_qtz1 = transpose(xdec_qtz); % o sinal original em uma coluna
xdec = bi2de(x_exp); % ou xdec = bi2de(transpose(xdec_qtz1)) 
xdec_passo = xdec * passo;
xdec_rec = xdec_passo - ((Vpp/2) - (passo/2));
f_PB = fir1(50, 10e3/Fs);
yy = conv(f_PB, xdec_rec);
%plot(yy);
%sound(yy, Fs)
%sound(xdec_rec, Fs)
%% analizando o espctro
X_f = fftshift(fft(x_n)/length(x_n));
Xdec_f = fftshift(fft(yy)/length(yy));
f = -Fs/2:1/(length(Xdec_f)/Fs):Fs/2-(1/length(Xdec_f)/Fs);%-1;
figure(2)
subplot(211)
plot(yy)
subplot(212)
plot(f, (abs(Xdec_f)))
xlim([-10e3 10e3])