clear all; close all; clc;

%ucitavanje signala i f0 iz fajla:
[x, fs] = audioread('recenica 10.wav');
x1 = load('f0_recenica 10.mat');
f0 = x1.f0;
tOsa = 0:length(x)-1;

%Parametri zadatka:
tStep = 0.03;
count = round(tStep*fs);
R = 0.5*count; %prozor za R koji je polovina originalnog prozora

%Generisanje pobudnog signala
e = [];
location = 0;

for i = 1:length(f0)
    pobuda=zeros(1,R);
  if location>R %provera za mesto sledeceg odbirka
      location=location-R;
  else %dobro mesto dalja obrada
    if isnan(f0(i))
      if i>1
        if isnan(f0(i-1))
          pobuda=0.01*randn(1,R);
        else
          pobuda(location)=1; 
          pobuda((location+1):R)=0.01*randn(1,R-location);
        end
      else %popunjavanje prvog elementa
        pobuda=0.01*randn(1,R);
      end
    else
        if i>1
            if isnan(f0(i-1))
                pobuda(1:fix(fs/f0(i)):R)=1;
                location= fix(fs/f0(i) - R + find(pobuda, 1, 'last'));
            else 
                pobuda(location:fix(fs/f0(i)):R)=1;
                location= fix(fs/f0(i) - R + find(pobuda, 1, 'last'));
            end
        else %popunjavanje prvog elementa
            pobuda(location:fix(fs/f0(i)):R)=1;
            location= fix(fs/f0(i) - R + find(pobuda, 1, 'last'));
        end
    end
  end
    e=[e pobuda];
end
eOsa = 0:length(e)-1;
figure, plot(eOsa,e); %plotovanje pobudnog signala

%Rekonstrukcija signala
signal = [];
for i = 1:(count/2):length(x)-count
    y = x(i:i+count-1); %Prozorovanje
    yW = y.*hamming(count); %Haming mnozenje
    
    p = 35;
    [A, G] = autolpc(yW, p); %Racunanje koeficijenata LPC
    
    Gain = G/(sqrt(sum(e(i:i+R-1).^2))+0.01);
    s(i:i+R-1) = filter(Gain, A, e(i:i+R-1));
    signal = [signal s(i:i+R-1)];
end

sOsa = 0:length(signal) - 1;
figure, plot(sOsa, signal, tOsa, x); %plotovanje signala origala i rekonstruisanog
sound(signal, fs); %pustanje signala




