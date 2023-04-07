clear all, close all, clc
%ucitavanje signala:
[x, fs] = audioread('./recenica 10.wav');

%Obrada u vremenskom domenu (normmalizacija signala i uklanjanje
%jednosmerne komponente)

x = x./max(abs(x));
DC = mean(x);
x = x - DC;

%prikaz u vremenskom domenu
t = (0:length(x)-1)*(1/fs);
T=1/fs;
op= T:T:length(x)*T;

%parametri
tStep = 0.03;
count = round(tStep*fs);
p0 = 2*10^-5;
eThreshold = 55;
zcrThreshold = 70;

%nizovi
ZCR=[];
RMS = [];
f0 = [];
L=[];

%kreiranje filtra za izvlacenje zvucnih delova zbog nalazenja ZCR-a
wn1=100/(fs/2);
wn2=650/(fs/2);
M=50;
N=2*M;
h=fir1(N, [wn1 wn2], rectwin(N+1));
tI = 10*10^-3;

for i = 1:(count/2):length(x)-count
    y = x(i:i+count-1);
    RMS(end+1) = 20*log(rms(y)/p0);
    ZCR(end+1) = zcr(y);
    if RMS(end) >= eThreshold && ZCR(end) <= zcrThreshold
        L(end+1) = 1; %beleze se zvucni signali
        y1 = filter(h, 1, y);
        yCorrelated = xcorr(y1);
        [p, locs] = findpeaks(yCorrelated, 'MinPeakDistance', 200);
        pmax = find(p == max(p));
        f0(end+1) = 1/((locs(pmax)-locs(pmax-1))*(tStep/count)); %nalazenje osnovne frekvencije
    else
        L(end+1) = 0;
        f0(end+1) = 0;

    end

    
end

%Plotovanje RMS, ZCR i osnovne frekvencije
tNormal = tStep/2;
tOsa = (0:length(RMS)-tNormal)*tNormal;
maxRMS = max(abs(RMS));
RMSNovo = RMS/maxRMS;
maxF0 = max(abs(f0));
f0Novo = f0/maxF0;
figure, plot(tOsa, L, 'r', t, x, 'b', tOsa, RMSNovo, 'g', tOsa, f0Novo, 'black'), ylim([-1.1, 1.1]);

r = 0:tStep/2:(length(RMS)-1)*tStep/2;
%iscrtavanje
figure, plot(r, RMS)
title("RMS")
figure, plot(r, ZCR)
title("ZCR")
figure, plot(op, x, r, L)
ylim([-1.1 1.1])
title("Zvucnost")
figure,plot(r, f0)
title("f0")

%ZCR Funkcija
function zcr = zcr(x)
    zcr = 0;
    for br = 1:length(x)-1
        zcr = zcr + abs(sign(x(br)) - sign(x(br+1)))/2;
    end
end