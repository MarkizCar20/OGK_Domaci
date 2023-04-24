clear all; close all; clc;

%ucitavanje signala i f0 iz fajla:
[x, fs] = audioread('recenica 10.wav');
x1 = load('f0_recenica 10.mat');
f0 = x1.f0;

%Parametri zadatka:
tStep = 0.03;
count = round(tStep*fs);
R = 0.5*count; %prozor za R koji je polovina originalnog prozora

%Generisanje pobudnog signala
e = [];
location = 0;

for i = 1:length(f0)
    pobuda=zeros(1,R);
  if location>R
      location=location-R;
  else
    if isnan(f0(i))
      if i>1
        if isnan(f0(i-1))
          pobuda=0.01*randn(1,R);
        else
          pobuda(location)=1;
          pobuda((location+1):R)=0.01*randn(1,R-location);
        end
      else
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
        else
            pobuda(location:fix(fs/f0(i)):R)=1;
            location= fix(fs/f0(i) - R + find(pobuda, 1, 'last'));
        end
    end
  end
    e=[e pobuda];
end

e_osa = 0:length(e)-1;
figure, plot(e_osa,e);




