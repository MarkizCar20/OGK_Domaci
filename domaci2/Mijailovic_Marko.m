clear all, close all, clc
%ucitavanje signala:
[x, fs] = audioread('./domaci2/recenica 10.wav');

%Obrada u vremenskom domenu (normmalizacija signala i uklanjanje
%jednosmerne komponente)
x = x./max(abs(x));
DC = mean(x);
x = x - DC;

%prikaz u vremenskom domenu
t = (0:length(x)-1)*(1/fs);
%figure, plot(t, x, 'b'), ylim([-1.1 1.1])

%parametri
tStep = 0.03;
count = round(tStep*fs);
p0 = 2*10^-5;
eThreshold = 64;
zcrThreshold = 70;
%nizovi
RMS = [];
ZCR=[];
L=[];
f0=[];

%kreiranje filtra za izvlacenje zvucnih delova zbog nalazenja ZCR-a
wn1=80/(fs/2);
wn2=600/(fs/2);
M=50;
N=2*M;
h=fir1(N, [wn1 wn2], rectwin(N+1));

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

%ZCR Funkcija
function zcr = zcr(x)
    zcr = 0;
    for br = 1:length(x)-1
        zcr = zcr + abs(sign(x(br)) - sign(x(br+1)))/2;
    end
end