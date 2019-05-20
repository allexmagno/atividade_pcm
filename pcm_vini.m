clear all;
close all;
clc;
%%
Vpp = 5; % escala de tensao
k = 12; % quantidade de bits do quantizador

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
pause(5) % intervalo para reproducao do audio

f_cut = 8e3;
filtro_PB = fir1(100, (f_cut*2)/Fs_in);

fig = 2;
passo_q = Vpp/(2^k); % passo de quantizacao

x_quant_1 = x_n + ((Vpp/2)-(passo_q/2)); % deslocamento pra cima
x_quant_2 = x_quant_1/passo_q; % joga pro nivel de tensao
x_quant_3 = round(x_quant_2); % arredonda
aux1 = x_quant_3 == 2^k; % gera vetor de "bool" pra 8
aux2 = x_quant_3 == -1; % gera vetor de "bool" pra -1
x_quant_4 = x_quant_3 - aux1; % ajusta max ate 2^k-1
x_quant = x_quant_4 + aux2; % ajusta min em 0

x_bin = de2bi(x_quant);
[x, y] = size(x_bin);
y = y+1; % 13 bits
%     x_bin = reshape(transpose(x_bin_1), 1, linha*coluna);
%%
% compressao !!
i = 1;
x_comp = zeros(x, 8);
while i<=x
    j = 1;
    while j <= 7
       if x_bin(i, j) == 1
           break;
       end
       j = j + 1;
    end
   if j == 1
       b = j+1;
%        c = [x_bin(i, 1) 1 1 1 x_bin(i, b:b+3)];
       c = [0 1 1 1 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 2
       b = j+1;
%        c = [x_bin(i, 1) 1 1 0 x_bin(i, b:b+3)];
       c = [0 1 1 0 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 3
       b = j+1;
%        c = [x_bin(i, 1) 1 0 1 x_bin(i, b:b+3)];
       c = [0 1 0 1 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 4
       b = j+1;
%        c = [x_bin(i, 1) 1 0 0 x_bin(i, b:b+3)];
       c = [0 1 0 0 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 5
       b = j+1;
%        c = [x_bin(i, 1) 0 1 1 x_bin(i, b:b+3)];
       c = [0 0 1 1 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 6
       b = j+1;
%        c = [x_bin(i, 1) 0 1 0 x_bin(i, b:b+3)];
       c = [0 0 1 0 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   elseif j == 7
       b = j+1;
%        c = [x_bin(i, 1) 0 0 1 x_bin(i, b:b+3)];
       c = [0 0 0 1 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   else
       b = j+1;
%        c = [x_bin(i, 1) 0 0 0 x_bin(i, b:b+3)];
       c = [0 0 0 0 x_bin(i, b:b+3)];
       x_comp(i,:) = c;
   end
   i = i + 1;
end


% end comp ...
[linha_comp, coluna_comp] = size(x_comp);
x_send = reshape(transpose(x_comp), 1, linha_comp*coluna_comp);

%%
% descompressao
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

%%
% voltando ...
%     x_dec_1 = reshape(x_bin, coluna, linha);
%     x_dec_2 = transpose(x_dec_1);
% x_dec_1 = bi2de(x_expand(:,2:13), 'left-msb');
% i = 1;
% while i<=linha_comp
%     if x_expand(i, 1) == 1
%         x_dec_1(i) = -x_dec_1(i);
%     end
%     i = i + 1;
% end
x_dec_1 = bi2de(x_expand, 'left-msb');
x_dec = transpose(x_dec_1);
x_n_rec_1 = x_dec * passo_q;
x_n_rec = x_n_rec_1 - ((Vpp/2)-(passo_q/2));

x_n_rec_fft = fftshift(fft(x_n_rec)/length(x_n_rec));
y_n = filter(filtro_PB, 1, x_n_rec);
y_fft = fftshift(fft(y_n)/length(y_n));


figure(fig)
suptitle(['Quantização de ',num2str(k),' bits'])
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
pause(5) % intervalo para reproducao do audio
audio_recuperado_filtrado = audioplayer(y_n', Fs_in);
% play(audio_recuperado_filtrado)
% pause(5) % intervalo para reproducao do audio
