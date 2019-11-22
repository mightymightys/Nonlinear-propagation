function [spcgrm,spcgrmtime] = myspectrogram(data,t,spcgrmwinwidth,winchoice)

nPoints = length(data);
winwidth = round(spcgrmwinwidth /(t(2)-t(1)));

if winchoice == 1
    window = gausswin(winwidth);
elseif winchoice ==2
    window = hanning(winwidth);
else
    window = hanning(winwidth);
end

steps=(0:ceil(winwidth/8):nPoints-winwidth);

spcgrm = zeros(length(steps),nPoints);
spcgrmtime = zeros(length(steps),1);

n=1;
for i=steps
    filter=[zeros(i,1); window; zeros(nPoints-winwidth-i,1)];
    dmp = ift(data.*filter');
    spcgrm(n,:) = (abs(dmp)).^2;
    spcgrmtime(n) = t(i+round(winwidth/2));
    n=n+1;
end



%figure; imagesc(spectrogramtime,nu,spectrogram'.*conj(spectrogram'))
%ylim([0,1])
