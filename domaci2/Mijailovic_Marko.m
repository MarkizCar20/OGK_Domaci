clear all, close all, clc
%ucitavanje signala:
[x, fs] = audioread('./recenica 10.wav');

%Obrada u vremenskom domenu (normmalizacija signala i uklanjanje
%jednosmerne komponente)
x = x./max(abs(x));
DC = mean(x);
x = x - DC;

%prikaz u vremenskom domenu i spektrogram
t = (0:length(x)-1)*(1/fs);

%parametri
tStep = 0.03;
count = round(tStep*fs);
p0 = 2*10^-5;
eThreshold = 78;
zcrThreshold = 120;

%nizovi
ZCR=[];
RMS = [];
f0 = [];
L=[];

%kreiranje filtra za izvlacenje zvucnih delova zbog nalazenja ZCR-a
wn1=100/(fs/2);
wn2=650/(fs/2);
M=47;
N=2*M;
h=fir1(N, [wn1 wn2], rectwin(N+1));

for br = 1:(count/2):length(x)-count
    y = x(br:br+count-1);
    RMS(end+1) = 20*log10(rms(y)/p0);
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
maxZCR = max(abs(ZCR));
ZCRNovo = ZCR/maxZCR;
figure, plot(tOsa, L, 'r', t, x, 'b', tOsa, RMSNovo, 'g', tOsa, f0Novo, 'black'), ylim([-1.1, 1.1]);
figure, plot(tOsa, L, 'r', 'LineWidth', 1.3),title('Promena nivoa signala'), ylim([-0.5, 1.1]);
figure, plot(tOsa, RMSNovo, 'r', 'LineWidth', 1.3),title('RMS'), ylim([0.4, 1.1]);
figure, plot(tOsa, f0Novo, 'r', 'LineWidth', 1.3),title('Promena osnovne ucestanosti'), ylim([0, 1.1]);
figure, plot(tOsa, ZCRNovo, 'r', 'LineWidth', 1.3),title('Promena ZCR'), ylim([0, 1.1]);    

%ZCR Funkcija
function zcr = zcr(x)
    zcr = 0;
    for br = 1:length(x)-1
        zcr = zcr + 1/2*abs(sign(x(br))-sign(x(br+1)));
    end
end