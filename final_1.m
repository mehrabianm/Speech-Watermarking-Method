function final_1()

clear, close all, clc; 

R = {'wav-file/original/2.wav'};

W = {'wav-file/watermark/w2.wav'};

E = {'wav-file/extract/attack/w2.wav'};


%res = {'index' 'SNR' 'SSIM' 'BER' 'NCC'};
for i = 1 : 1
    fnl(i, :) = resultSave(cell2mat(R(i)), cell2mat(W(i)),...
        cell2mat(E(i)), i);
   % res(i+1,:) = [i, num2cell(fnl(i,:))];
end
%xlswrite('xls/noAttackResult', res)
end


%% Result
function fnl = resultSave(wav, wtk, extFile, counterFig)
% Henon_map(1.4,0.3);
[host, f] = audioread (wav);
original = host;
% z = var(host);
x = host';
y = chiSquare(host);
plot(host);
h = figure;
subplot(1, 1, 1);
plot(host);
title('Original signal');
h = figure;
subplot(1, 1, 1);
% [h,p] = chi2gof(host,'Alpha',0.01)
% [m,v]=chi2stat(host);
hist(host, 1000);
title('Histogram of the original signal');
host = (host + 1) * (32768);

Fn = 512;

L = floor(length(host)/ Fn);
extra = host(L * Fn + 1 : end, 1)';
scrbl = host(:, 1);         %scramble audio bits
% B = abs(sum(sum(host - scrbl)));

alpha = 5;
for i = 1 : L
    F(i, :) = scrbl((i - 1) * Fn + 1 : i * Fn);
end

%% Insertion
tic
timer_embed_watermark_gen_started = tic

[S, lenFramNum, embed_wtk] = insertion(F, L, alpha);

toc
timer_embed_watermark_gen_finished = toc

tic
timer_encrypt_started = tic
S = [S, extra];
c = S';
c = de2bi(c, 'left-msb')';
audioBit = reshape(c, 1, size(c, 1) * size(c, 2));
[~, bitX] = lajesticMap(length(audioBit));
encriptS = xor(bitX, audioBit);
% result = reshape();
for i = 1 : length(encriptS) / 16
    S_new(i)  = bin2dec(num2str(encriptS(16 * (i - 1) + 1 : 16 * i)));% watermarked host
end
S = S_new / 32768 - 1;

toc
timer_encrypt_finished = toc

h = figure;
subplot(1, 1, 1);
plot(S, 'g');
title('Encrypted watermarked signal');
axis([0 8*10^4 -1.2 1.2])
h = figure;
subplot(1, 1, 1);
hist(S, 1000, 'g');
title('Histogram of Encrypted watermarked signal');
axis([-1.2 1.2 0 1600])

nameFIG = ['figure/noAttack/host_', num2str(counterFig)...
    , '_watermark_', num2str(counterFig),'.jpg'];
saveas(h, nameFIG);

audiowrite(wtk, S, f);
tic
timer_decrypt_authenticate_started = tic

fnl = finalExt_20azar(wav, wtk, counterFig, extFile, lenFramNum, timer_embed_watermark_gen_finished, embed_wtk);
% fnl = extractionWithAttack(wav, wtk, counterFig, extFile, lenFramNum, timer_embed_watermark_gen_finished, embed_wtk);
toc
timer_decrypt_authenticate_finished = toc

time_alarm=['timer_embed_watermark_gen_finished' 'timer_encrypt_finished' 'timer_decrypt_authenticate_finished'];

time_s=[timer_embed_watermark_gen_finished timer_encrypt_finished timer_decrypt_authenticate_finished];
total_time = timer_embed_watermark_gen_finished + timer_encrypt_finished + timer_decrypt_authenticate_finished;
disp('total_time');
disp(total_time);
% disp(fnl(6));
end
%% Insertion Function
function [s, lenFramNum, wtk] = insertion(f, l, alpha)

wtk = [];

lenFramNum = length(de2bi(l));
save lenFN lenFramNum;
for i = 1 : l
    Fi = f(i, :);
    A = zeros(1, floor(size(f, 2) / 4));
    [ca, cd] = lwt(Fi, 'lazy');
    [ca1, cd1] = lwt(ca, 'lazy');
    fram_num = de2bi(i, lenFramNum, 'left-msb');
    cd1 = insertionFarmeNum(cd1, fram_num);
    miu = mean(ca1);
    gama = miu + alpha;
    ca1_1 = find(ca1 > gama);
    A(ca1_1) = 1;
    [~, bitX] = lajesticMap(length(A));
    encriptA = xor(A, bitX);
    cd = insertionCaDetail(cd, encriptA);
    
    ca = ilwt(ca1, cd1, 'lazy');
    fiFinal(i, :) = ilwt(ca, cd, 'lazy');
    
    wtk = [wtk, fram_num, fram_num, encriptA];
    
end
s = reshape(fiFinal', 1, size(fiFinal, 1) * size(fiFinal, 2));


end

%% insertion frame number in cd1
function [cd] = insertionFarmeNum(c, f)
counter = 1;
mess = f;
cd = c;
for i = 1 : length(mess) * 2
    cd(i) = bitset(c(i), 1, mess(counter));
    counter = counter + 1;
    if counter == length(mess) + 1
        counter = 1;
    end
end

end

%% insertion CA Detail in CD
function [cd] = insertionCaDetail(c, mess)
counter = 1;
cd = c;
if length(c) >= length(mess)
    for i = 1 : length(mess)
        cd(i) = bitset(c(i), 1, mess(counter));
        counter = counter + 1;
    end
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
    if x(i) > alpha
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
%% Chi-square test
function [y]=chiSquare(x)
m = mean(x);
l = length(x);
sum = 0;
for i = 1 : l
    sum = sum + (x(i)-m)^2/m 
end
y = sum;
end
%%
% function Henon_map(a,b)
% %This function takes in the alpha and beta values for the Henon map and
% %iterates (0.1,0) 6000 times.  It disregards the first 50 iterates and
% %graphs the rest in the Cartesian plane.
% N=6000;
% x=zeros(1,N);
% y=zeros(1,N);
% x(1)=0.1;
% y(1)=0;
% for i=1:N
%     x(i+1)=1+y(i)-a*(x(i))^2;
%     y(i+1)=b*x(i);
% end
% axis([-1,2,-1,1])
% plot(x(50:N),y(50:N),'.','MarkerSize',1);
% fsize=15;
% set (gca,'xtick',[-1:1:1],'FontSize',fsize)
% set (gca,'ytick',[-1:1:2],'FontSize',fsize)
% xlabel('\itx','FontSize',fsize)
% ylabel('\ity','FontSize',fsize)
% end


