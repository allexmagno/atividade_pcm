clear all;
close all;
clc;
%%
%pkg load signal
%pkg load communications
Vpp = 2; % escala de tensao
k = 13; % quantidade de bits do quantizador

[x_n,Fs_in] = audioread('you-dont-know-how.wav');
x_n = x_n(:,1); % pega somente a primeira coluna

audio_original = audioplayer(x_n, Fs_in);
play(audio_original)

end_t = (length(x_n)/Fs_in-(1/Fs_in));
t = [0:1/Fs_in:end_t];
f = [-Fs_in/2:1/end_t:Fs_in/2];
x_n_fft = fftshift(fft(x_n)/length(x_n));
figure(1)
subplot(121)
plot(t, x_n)
title('Áudio original no tempo')
xlabel('Tempo [s]')
ylabel('Amplitude [V]')
subplot(122)
plot(f, abs(x_n_fft))
title('Espectro do áudio original')
xlabel('Frequência [Hz]')
ylabel('Amplitude [V]')
pause(3) % intervalo para reproducao do audio

passo_q = Vpp/(2^k); % passo de quantizacao
x_quant_2 = x_n/passo_q; % desloca para o nivel de tensão
x_quant = round(x_quant_2); % arredonda para valor inteiro

negativos = x_quant < 0; % guarda valores negativos
[linhas, colunas] = size(negativos);
x_quant = abs(x_quant); % deixa positivo para executar o de2bi
x_bin_1 = de2bi(x_quant, 13, 'left-msb');

% Observando o sinal transmitido
figure(2)
[x, y] = size(x_bin_1);
x_bin_linha = reshape(transpose(x_bin_1), 1, x*y); % Para colocar tudo em uma unica linha 
subplot(311)
filtro_NRZ = ones(1, 100); % filtro
%tx = filter(filtro_NRZ, 1, x_bin_linha);
tx = filter(filtro_NRZ, 1, upsample(x_bin_linha, 100));
plot(tx);
ylim([-0.5 2])
xlim([2.3e6 2.31e6])
title('Transmissão do sinal original 13 bits')


% Loop que recupera os valores que eram negativos
i = 1;
while i<=linhas
    if negativos(i) == 1
        x_bin_1(i, 1) = 1;
    end
    i = i + 1;
end
x_bin = x_bin_1;

%% Compressao

[x, y] = size(x_bin);
x_comp = zeros(x, 8);
i = 1;
while i<=x
    j = 2;
    while j <= 8
       if x_bin(i, j) == 1
           break;
       end
       j = j + 1;
    end
   if j == 2
       b = j+1;
       c = [x_bin(i, 1) 1 1 1 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 3
       b = j+1;
       c = [x_bin(i, 1) 1 1 0 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 4
       b = j+1;
       c = [x_bin(i, 1) 1 0 1 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 5
       b = j+1;
       c = [x_bin(i, 1) 1 0 0 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 6
       b = j+1;
       c = [x_bin(i, 1) 0 1 1 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 7
       b = j+1;
       c = [x_bin(i, 1) 0 1 0 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 8
       b = j+1;
       c = [x_bin(i, 1) 0 0 1 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   else
       b = j+1;
       c = [x_bin(i, 1) 0 0 0 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   end
   i = i + 1;
end

% Gera vetor de 0 e 1 que serao transmitidos
[linha_comp, coluna_comp] = size(x_comp);
x_send = reshape(transpose(x_comp), 1, linha_comp*coluna_comp);

%Visualizando o Sinal
subplot(312)
tx_send = filter(filtro_NRZ, 1, upsample(x_send, 100));
plot(tx_send);
ylim([-0.5 2])
xlim([2.3e6 2.31e6])
title('Transmissão do sinal Comprimido 8 bits')


%% Descompressao
x_rec = transpose(reshape(x_send, coluna_comp, linha_comp));
x_expand = zeros(linha_comp, 13);
i = 1;
while i<=linha_comp
    val = bi2de(x_rec(i,2:4), 'left-msb');
    if val == 0
        c = [x_rec(i,1) 0 0 0 0 0 0 0 x_rec(i,5:8) 1];
        x_expand(i,:) = c;
    elseif val == 1
        c = [x_rec(i,1) 0 0 0 0 0 0 1 x_rec(i,5:8) 1];
        x_expand(i,:) = c;
    elseif val == 2
        c = [x_rec(i,1) 0 0 0 0 0 1 x_rec(i,5:8) 1 0];
        x_expand(i,:) = c;
    elseif val == 3
        c = [x_rec(i,1) 0 0 0 0 1 x_rec(i,5:8) 1 0 0];
        x_expand(i,:) = c;
    elseif val == 4
        c = [x_rec(i,1) 0 0 0 1 x_rec(i,5:8) 1 0 0 0];
        x_expand(i,:) = c;
    elseif val == 5
        c = [x_rec(i,1) 0 0 1 x_rec(i,5:8) 1 0 0 0 0];
        x_expand(i,:) = c;
    elseif val == 6
        c = [x_rec(i,1) 0 1 x_rec(i,5:8) 1 0 0 0 0 0];
        x_expand(i,:) = c;
    else
        c = [x_rec(i,1) 1 x_rec(i,5:8) 1 0 0 0 0 0 0];
        x_expand(i,:) = c;
    end
    i = i + 1;
end

subplot(313)
[x, y] = size(x_expand);
rx_bin_linha = reshape(transpose(x_expand), 1, x*y); % Para colocar tudo em uma unica linha 
rx_exp = filter(filtro_NRZ, 1, upsample(rx_bin_linha, 100));
plot(rx_exp);
ylim([-0.5 2])
xlim([2.3e6 2.31e6])
title('Transmissão do sinal expandido 13 bits')

%% Recupera valores "originais"
x_dec_1 = bi2de(x_expand(:,2:13), 'left-msb'); % bi2de apenas unsigned

% Loop que recupera os valores negativos
i = 1;
while i<=linha_comp
    if x_expand(i, 1) == 1
        x_dec_1(i) = -x_dec_1(i);
    end
    i = i + 1;
end

x_dec = transpose(x_dec_1);
x_n_rec = x_dec * passo_q; % retorna aos valores originais

f_cut = 3.4e3;
filtro_PB = fir1(100, (f_cut*2)/Fs_in);
y_n = filter(filtro_PB, 1, x_n_rec);
x_n_rec_fft = fftshift(fft(x_n_rec)/length(x_n_rec));
y_fft = fftshift(fft(y_n)/length(y_n));


figure(2)
subtitle(['Quantização de ',num2str(k),' bits'])
subplot(221)
plot(t, x_n_rec)
title('Áudio no tempo')
xlabel('Tempo [s]')
ylabel('Amplitude [V]')
subplot(222)
plot(f, abs(x_n_rec_fft))
title('Espectro do áudio')
xlabel('Frequência [Hz]')
ylabel('Amplitude [V]')
subplot(223)
plot(t, y_n)
title('Áudio com filtro 8 KHz')
xlabel('Tempo [s]')
ylabel('Amplitude [V]')
subplot(224)
plot(f, abs(y_fft))
title('Espectro do áudio filtrado')
xlabel('Frequência [Hz]')
ylabel('Amplitude [V]')

audio_recuperado = audioplayer(x_n_rec', Fs_in);
play(audio_recuperado)
pause(3) % intervalo para reproducao do audio
audio_recuperado_filtrado = audioplayer(y_n', Fs_in);
play(audio_recuperado_filtrado)
pause(4) % intervalo para reproducao do audio
