function fnl = extractionWithAttack(wav, wtk, counterFig, extFile, lenFramNum, timer_embed_finished, embed_wtk)

[main, f] = audioread (wav);
% main = resample(main, 11000, f);% for resample attack
% main = upsample(main, 2);% for upsample attack
% main = downsample(main, 2);% for downsample attack
[host, f] = audioread (wtk);
wtmHost = (host);

% -------------- Subsitiution Attack ----------- %

% host(21001: 21500) = main(21001: 21500);
% % wtmHost(21001: 24000) = main(21001: 24000);
% 
% figure;
% subplot(1, 1, 1)
% plot(host);
% % title('Encrypted Watermarked Signal');
% axis([0 (length(host)+10^4) -1.2 1.2 ])

% ----------------- End Attack ----------------- %

% --------------- Mute Attack ------------------ %

% host(28001: 28500) = 0;
% wtmHost(28001: 28500) = 0;
% 
% figure;
% subplot(1, 1, 1)
% plot(host);
% % title('Encrypted Watermarked Signal');
% axis([0 (length(host)+10^4) -1.2 1.2 ])
% ----------------- End Attack ----------------- %

% ------------- Insertion Attack --------------- %

% a = 0.1;
% b = -0.1;
% r = a + (b-a).*rand(50001, 1);
% 
% host(2000: 52000) = r;
% wtmHost(2000: 52000) = r;
% 
% figure;
% subplot(1, 1, 1)
% plot(host);
% % title('Encrypted Watermarked Signal');
% axis([0 (length(host)+10^4) -1.2 1.2 ])

%--------------End Attack------------------------%
%--------------Mixed Attack----------------------

% host(22657: 23056) = 0;
% wtmHost(22657: 23056) = 0;
% % Mixed By
% a = 0.5;
% b = -0.5;
% r = a + (b-a).*rand(600, 1);
% 
% host(36625: 37224) = r;
% wtmHost(36625: 37224) = r;
% 
% h = figure;
% subplot(1, 1, 1)
% plot(host);
% % title('Encrypted Watermarked Signal');
% axis([0 (length(host)+10^4) -1.2 1.2 ])

%--------------End Attack------------------------%

host = floor((host + 1) * (32768));
c = host';
c = de2bi(c, 'left-msb')';
audioBit = reshape(c, 1, size(c, 1) * size(c, 2));
[~, bitX] = lajesticMap(length(audioBit));
encriptS = xor(bitX, audioBit);
for i = 1 : length(encriptS) / 16
    S_decript(i)  = bin2dec(num2str(encriptS(16 * (i - 1) + 1 : 16 * i))); % watermarked host
end

wtk_decrypt = S_decript/32768 - 1;
audiowrite('wtk_dec.wav', wtk_decrypt, f);
[wtk_decrypt, ~] = audioread( 'wtk_dec.wav');


%--------------------MP3 Attack-------------------

S_decript = S_decript / 32768 - 1;
de_wtk = 'wav-file/mp3/de_wtk1.wav';
audiowrite(de_wtk, S_decript, f);
[d,sr] = audioread(de_wtk);
mp3write(d,sr,'wav-file/mp3/1.mp3');
% Read it back again
[S_decript,f] = mp3read('wav-file/mp3/1.mp3');
S_decript = floor((S_decript + 1) * (32768));

%---------------------End Attack------------------
%----------------Gaussian Noise Attack 1------------

% S_decript = S_decript / 32768 - 1;
% de_wtk = 'wav-file/noisy/de_wtk1.wav';
% audiowrite(de_wtk, S_decript, f);
% [d,~] = audioread(de_wtk);
% noise = 0.01*randn(1,length(d))';
% S_decript = d + noise;
% S_decript = floor((S_decript + 1) * (32768));
%----------------Gaussian Noise Attack 2------------
% S_decript = S_decript / 32768 - 1;
% host = awgn(S_decript, 25);
% S_decript = floor((host + 1) * (32768));
%---------------------End Attack------------------
%----------------Resample Attack---------------- %
% S_decript = S_decript / 32768 - 1;
% host = resample(S_decript, 11000, f);
% S_decript = floor((host + 1) * (32768));
% ----------------- End Attack ----------------- %
%----------------Upsampling Attack---------------- %
% S_decript = S_decript / 32768 - 1;
% host = upsample(S_decript, 2);
% S_decript = floor((host + 1) * (32768));
% ----------------- End Attack ----------------- %
%----------------Downsampling Attack------------ %
% S_decript = S_decript / 32768 - 1;
% host = downsample(S_decript, 2);
% S_decript = floor((host + 1) * (32768));
% ----------------- End Attack ----------------- %
%------------Low-Pass Filtering Attack---------- %
% S_decript = S_decript / 32768 - 1;
% host = lowpass(S_decript, 4000, f);
% S_decript = floor((host + 1) * (32768));
% ----------------- End Attack ----------------- %
%------------High-Pass Filtering Attack 1---------- %
% S_decript = S_decript / 32768 - 1;
% host = highpass(S_decript, 75, f);
% S_decript = floor((host + 1) * (32768));
% ----------------- End Attack ----------------- %
%------------High-Pass Filtering Attack 2--------- %
% S_decript = S_decript / 32768 - 1;
% d = fdesign.highpass('Fst,Fp,Ast,Ap',0.15,0.25,60,1);
% designmethods(d);
% Hd = design(d,'equiripple');
% fvtool(Hd);
% host = filter(Hd,S_decript);
% plot([S_decript  host]);
% S_decript = floor((host + 1) * (32768));
% ----------------- End Attack ----------------- %
Fn = 512;
L = floor(length(host)/ Fn);
extra = S_decript(L * Fn + 1 : end)';
alpha = 5;

for i = 1 : L
    F(i, :) = S_decript((i - 1) * Fn + 1 : i * Fn);
end

[ext, BER, ext_wtk] = extraction(F, L, alpha, counterFig);
% ext = [ext, extra'];%%Normal situation
ext = [ext, extra];%% Just for MP3 & noise 
ext = ext / 32768 - 1;
audiowrite(extFile, ext, f);
[ext, f] = audioread(extFile);

% ----------------------- Analyz ----------------------- %
% SNR1=Cal_SNR(wtk_decrypt, ext);
SNR=Cal_SNR(main, ext);
SegSNR = Cal_SegSNR(main, ext, L, Fn);
SSIM = ssim(main, ext);
BER = (sum(abs(embed_wtk - ext_wtk)) / length(embed_wtk));
NC = Cal_NC(main, ext);
PayLoad = (lenFramNum*2 + Fn/4) * L/timer_embed_finished;
fnl = [SNR, SSIM, BER, SegSNR, NC, PayLoad, timer_embed_finished];
disp('   SNR    SSIM    BER   SegSNR   NCC   PayLoad   timer_embed');
disp(fnl);

%% 

subplot 312
plot(wtmHost, 'g');
title('Encrypted and Watermaked signal');
axis([0 (length(wtmHost)+10^4) -1.2 1.2 ])
subplot 313
plot(ext, 'k');
title('Extracted signal');
axis([0 (length(ext)+10^4) -1.2 1.2 ]);

% nameFIG = ['figure/InsertionAttack/org_', num2str(counterFig)...
%     , '_enct_', num2str(counterFig),'.jpg'];
% saveas(h, nameFIG);
end

%% Extraction and Authentication Phase
function [ext, BER, ext_wtk] = extraction(f, l , alpha, counterFig)

authFlag = 0;
ext_wtk = [];
load lenFN;
for i = 1 : l
    Fi = f(i, :);
    A = zeros(1, floor(size(f, 2) / 4));
    [ca, cd] = lwt(Fi, 'lazy');
    [ca1, cd1] = lwt(ca, 'lazy');
    [ext_w1, extFrameNum, compersionFlag(i), authFlag(i), BER(i)] = ...
        extractionFarmeNum(cd1, i, lenFramNum);
    
        miu = mean(ca1);
        gama = miu + alpha;
        ca1_1 = find(ca1 > gama);
        A(ca1_1) = 1;
        [~, bitX] = lajesticMap(length(A));
        encriptA = xor(A, bitX);
        extCa = extractionCaDetail(cd, length(encriptA));
        ext_w2 = extCa;
        ext_wtk = [ext_wtk, ext_w1, ext_w2];
        diff = find(encriptA ~= extCa);
       if authFlag(i) ~= 1
        if length(diff) == 0
            checkCaDetailAuth(i) = 0;
           
        else
            checkCaDetailAuth(i) = 1;
            BER(i) = BER(i) + length(diff);
        end
        
    else
        checkCaDetailAuth(i) = 2;
        
    end
     ca = ilwt(ca1, cd1, 'lazy');
     fiFinal(i, 1 : size(f, 2)) = ilwt(ca, cd, 'lazy');
end
ext = reshape(fiFinal', 1, size(fiFinal, 1) * size(fiFinal, 2));
h = figure;
subplot 311
plot(authFlag, 'r');
axis([0 (length(authFlag)+100) -0.2 1.2 ]);
xlabel('Frame Number');
ylabel('Amplitude');
% nameFIG = ['figure/noAttack/Faild_Frame_Number_', num2str(counterFig), '.jpg'];
% saveas(h, nameFIG);
end

%% extraction frame number from cd1
function [ext_w1, extFrameNum, compersionFlag, authFlag, BER] = extractionFarmeNum(c, mainFramNum, lenFramNum)
counter = 1;
j = 1;
BER = 0;
for i = 1 : lenFramNum * 2
    frame_num(j, counter) = bitget(c(i), 1);
    counter = counter + 1;
    if counter == lenFramNum + 1
        counter = 1;
        j = 2; 
    end
end
if sum(frame_num(1, :) - frame_num(2, :)) == 0
    compersionFlag = 0;
    extFrameNum = bi2de(frame_num(1, :), 'left-msb');
    BER = sum(abs(frame_num(1, :) - frame_num(2, :)));
else
    extFrameNum = 0;
    compersionFlag = 1; 
end
if extFrameNum == mainFramNum
    authFlag = 0;
else
    authFlag = 1;
end
counter = 1;
for i=1 : 2
    for j=1 : lenFramNum
        ext_w1(1, counter) = frame_num(i,j);
        counter = counter +1;
    end
end
end

%% insertion CA Detail in CD
function [extCa] = extractionCaDetail(c, n)
counter = 1;
for i = 1 : n
    extCa(i) = bitget(c(i), 1);
    counter = counter + 1;
end
end

%% Chotic Map
function [indX, bitX] = lajesticMap(n)
r = 3.99;
alpha = 0.5;
x(1) = 0.67;
bitX(1) = 1;
for i = 2 : 3 * n
    x(i) = r * x(i - 1) * (1 - x(i -1));
    if x(i) > 0.5
        bitX(i) = 1;
    else
        bitX(i) = 0;
    end
end
[vX, indX] = sort(x(2 * n  + 1 : end));
bitX = bitX(2 * n  + 1 : end);
end
%% Calculating NC
function [NC] = Cal_NC(x, y)
sum1 = 0
sum2 = 0;
sum3 = 0;
main = x;
ext_1 = y;
mc = mean(mean(main));
ms = mean(mean(ext_1));
for i=1 : length(main)
       sum1 = sum1 + abs((main(i,1))-mc) * abs((ext_1(i,1))-ms);
       sum2 = sum2 + (main(i,1)-mc)^2;
       sum3 = sum3 + (ext_1(i,1)-ms)^2;
end

NC = sum1 / ((sqrt(sum2) * sqrt(sum3))); 
end
%% Calculating SegSNR
function [SegSNR] = Cal_SegSNR(x, y, I, J)

L1 = floor(length(x)/J);
extra1 = x(L1*J + 1 : end);
for i=1 : L1
  F1(i, :) = x((i-1)*J + 1 : i*J);
end

L2 = floor(length(y)/J);
extra2 = y(L1*J + 1 : end);
for i=1 : L2
  F2(i, :) = y((i-1)*J + 1 : i*J);
end

sum1=0;
sum2=0;
sum3=0;
for i=1 : I
    for j=1 : J
       sum1 = sum1 + F1(i,j)^2;
       sum2 = sum2 + (F1(i,j) - F2(i,j))^2;
    end
    if (sum2 ~= 0) && (sum1 ~= 0)
       sum3 = sum3 + 10 * log10(sum1/sum2);
    end
    sum1 = 0;
    sum2 = 0;
end
SegSNR = sum3 / I;
% for i=1 : I
%     for j=1 : J
%         if F1(i,j) ~= F2(i,j)
%             sum1 = sum1 + F1(i,j)^2/(F1(i,j) - F2(i,j))^2;
%         end
%     end
%     if sum1 ~= 0
%         sum2 = sum2 + log10(sum1);
%     end
%     sum1 = 0;
% end
% SegSNR = (10/I) * sum2;
end

%% Calculating SNR
function [SNR] = Cal_SNR(x, y)

L1 = length(x);
L2 = length(y);

sum1=0;
sum2=0;

for i=1 : L1
   
       sum1 = sum1 + x(i,1)^2;
       sum2 = sum2 + (x(i,1) - y(i,1))^2;
    
end
SNR = 10 * log10(sum1/sum2);

end
