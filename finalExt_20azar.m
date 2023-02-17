function fnl = finalExt_20azar(wav, wtk, counterFig, extFile, lenFramNum, timer_embed_finished, embed_wtk)

[main, ~] = audioread (wav);
[host, f] = audioread (wtk);%wtk is encrypted-watermarked signal file
wtmHost = (host);%host is encrypted-watermarked signal audio
% y = chiSquare(host);
host = floor((host + 1) * (32768));

c = host';
c = de2bi(c, 'left-msb')';
audioBit = reshape(c, 1, size(c, 1) * size(c, 2));
[~, bitX] = lajesticMap(length(audioBit));
decrypted_bit = xor(bitX, audioBit);
% result = reshape();
for i = 1 : length(decrypted_bit) / 16
    decrypted_S(i)  = bin2dec(num2str(decrypted_bit(16 * (i - 1) + 1 : 16 * i)));% watermarked decrypted host
    encrypted_S(i) = bin2dec(num2str(audioBit(16 * (i - 1) + 1 : 16 * i)));% watermarked encrypted host
end

%-----------------------------------------
Fn = 512;
L = floor(length(host)/ Fn);
extra1 = decrypted_S(L * Fn + 1 : end)';
extra2 = encrypted_S(L * Fn + 1 : end)';

alpha = 5;
for i = 1 : L
    F_d(i, :) = decrypted_S((i - 1) * Fn + 1 : i * Fn);% 2D array of decrypted signal
    F_e(i, :) = encrypted_S((i - 1) * Fn + 1 : i * Fn);% 2D array of encrypted signal
end

[ext, BER, ext_wtk] = extraction(F_d, L, alpha, counterFig);
[ext_e, BER_e, ext_wtk_e] = extraction(F_e, L, alpha, counterFig);

ext = [ext, extra1'];
ext = ext / 32768 - 1;
% y = chiSquare(ext);
audiowrite(extFile, ext, f);

%% Difference_Signal between Original and decrypted watermarked Signals
% diff_S = abs(main - ext');
% plot(main);
% plot(ext');
% plot(diff_S);
% title('Difference signal');
%axis([0 (length(main)+10^4) -1.2 1.2 ]);
%% Analyze (between Main and Encrypted watermarked signals
% SNR=Cal_SNR(main, wtmHost);
% SegSNR = Cal_SegSNR(main, wtmHost, L, Fn);
% SSIM = ssim(main, wtmHost);
% berSize = floor(Fn / 4) * L + 22 * L;
% BER_1 = (sum(abs(audioBit - (decrypted_bit))) / length(audioBit));
% BER_2 = (sum(abs(main - (wtmHost))) / length(main));
% BER = (sum(abs(embed_wtk - ext_wtk_e)) / length(embed_wtk));
% NC = Cal_NC(main, wtmHost);
% PayLoad = (lenFramNum*2 + Fn/4) * L/timer_embed_finished;
% fnl = [SNR, SSIM, BER, SegSNR, NC, PayLoad, timer_embed_finished, length(embed_wtk)];
% disp('   SNR        SSIM        BER        SegSNR        NCC          PayLoad       timer_embed        len_embed_wtk');
% disp(fnl);
%% Analyze (between Main and Decrypted watermarked signals
[main, ~] = audioread (wav);
SNR=Cal_SNR(main, ext');
SegSNR = Cal_SegSNR(main, ext', L, Fn);
SSIM = ssim(main, ext');
berSize = floor(Fn / 4) * L + 22 * L;
BER_1 = (sum(abs(audioBit - (decrypted_bit))) / length(audioBit));
BER_2 = (sum(abs(main - (ext'))) / length(main));
BER = (sum(abs(embed_wtk - ext_wtk)) / length(embed_wtk));
NC = Cal_NC(main, ext');
PayLoad = (lenFramNum*2 + Fn/4) * L/timer_embed_finished;
fnl = [SNR, SSIM, BER, SegSNR, NC, PayLoad, timer_embed_finished];
disp('   SNR    SSIM    BER   SegSNR   NC   PayLoad   timer_embed');
disp(fnl);
%-----------------------------------------------------------------
h = figure;
subplot(3, 1, 1);
plot(main);
title('Orginal signal');
axis([0 (length(main)+10^4) -1.2 1.2 ])
subplot(3, 1, 2);
plot(wtmHost, 'g');
title('Encrypted and watermaked signal');
axis([0 (length(wtmHost)+10^4) -1.2 1.2 ])
subplot(3, 1, 3);
plot(ext, 'k');
title('Extracted signal');
axis([0 (length(ext)+10^4) -1.2 1.2 ])

nameFIG = ['figure/noAttack/org_', num2str(counterFig)...
    , '_enct_', num2str(counterFig),'.jpg'];
saveas(h, nameFIG);

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
plot(authFlag, 'r');
title('Faild Frame Number')
axis([0 (length(authFlag)+100) -0.2 1.2 ])
nameFIG = ['figure/noAttack/Faild_Frame_Number_', num2str(counterFig), '.jpg'];
saveas(h, nameFIG);
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

%% extraction CA Detail from CD
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
%%
% function [indX, bitX] = lajesticMap(n)
% a = 0.35;
% alpha = 0.5;
% x(1) = 0.25;
% bitX(1) = 1;
% for i = 2 : 3 * n
%     if (x(i-1)>0 && x(i-1)<a)
%       x(i) = x(i-1)/a;
%     end
%     if (x(i-1)>=a && x(i-1)<0.5)
%       x(i) = (x(i-1)-a)/(0.5-a);
%     end
%     if (x(i-1)>=0.5 && x(i-1)<1)
%       x(i) = 1-x(i-1);
%     end
%     if x(i) > alpha
%         bitX(i) = 1;
%     else
%         bitX(i) = 0;
%     end
% end
% [vX, indX] = sort(x(2 * n  + 1 : end));
% bitX = bitX(2 * n  + 1 : end);
% end
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
%% Chi-square test
% function [y]=chiSquare(x)
% m = mean(x);
% l = length(x);
% sum = 0;
% for i = 1 : l
%     sum = sum + (x(i)-m)^2/m 
% end
% y = sum;
% end
