%Cetvrti Domaci iz Osnova Govornih komunikacija: Marko Mijailovic
clear all; close all; clc;

%Ucitavanje fajlova:
C = 6;
frameCount = 20;
mfccMatrix = [];

%%Proba sa malom bazom:
a = 5;
b = 205;

for cnt1 = 1:a
    for cnt2 = 1:b
        fileName = sprintf('Baza_Snimaka/broj_%d_%d.wav', cnt1, cnt2 );
        [y, fs] = audioread(fileName);
        Tw = (length(y)*1000)/(frameCount*fs);
        mfccVector = testMFCC(Tw, C, y, fs);
        mfccMatrix = [mfccMatrix; mfccVector];
    end
end

mfccCompare = [];
Threshold = [];

for i = 1:a*b
    for j = 1:a*b
        d = norm(mfccMatrix(i,:) - mfccMatrix(j, :));
        mfccCompare(i, j) = d;
        if d < 205;
            Threshold(i, j) = 1;
        else
            Threshold(i, j) = 0;
        end
    end
end

%Uspeh u odnosu na prag :)
successRate = zeros(1, round(max(mfccCompare, [], 'all')));

for k = 1:round(max(mfccCompare, [], 'all'));
    for i = 1:a*b
        for j = 1:a*b
            d = norm(mfccMatrix(i, :) - mfccMatrix(j, :));
            mfccCompare(i, j) = d;
            if d < k
                successRate(k) = successRate(k) + 1;
            end
        end
    end
end

successRate = successRate/(a*b);
figure(1), plot(1 - successRate);
figure(2), surf(mfccCompare, 'EdgeColor','none'), colorbar, view(2), colormap jet;
figure(3), surf(Threshold, 'EdgeColor','none'), colorbar, view(2), colormap jet;

%MFCC Funkcija
function mfccVector = testMFCC(Tw, C, speech, fs)
    
    Tw; %duzina prozora u ms
    Ts = Tw/2; %preklapanje (u nasem slucaju 50%)
    alpha = 0.97; %preemphasis koef.
    R = [300 3700]; % frekv. opseg
    M = 30; %broj filtara u banci
    C = 12; %broj kepstralnih koef.
    L = 22; % cepstral sine lifter param.

    %Hamming prozor
    hamming = @(N)(0.54 - 0.46*cos(2*pi*[0:N-1].'/(N-1)));
    
    %Izracunavanje MFCC koef.
    [MFCCs, FBEs, frames] = mfcc(speech, fs, Tw, Ts, alpha, hamming, R, M, C, L);
    [m, n] = size(MFCCs);
    if n == 38
        MFCCs(:, n+1) = [mean(MFCCs.').'];
    end

    %Preklapanje MFCC matrice
    mfccVector = reshape(MFCCs, 1, size(MFCCs,1)*size(MFCCs,2));

end

