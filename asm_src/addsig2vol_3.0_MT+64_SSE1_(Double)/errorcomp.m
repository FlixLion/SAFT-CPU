clear all()
close all
% Threads
threads = 4;
% Bildgröße x-Achse (Immer durch 4 teilbar. Vollständiges Bild kubisch)
x=40;
% Anzahl der AScans
blocksize= 10;
AscanLength = 3000;
% Anzahl der Messungen pro Konfiguration
samples = 300;
sArray = 1:1:samples;
% Zeitmessung
timesBlocked = zeros(samples,1);
timesBlockedOrig = zeros(samples,1);
errorArray = zeros(samples,1);

% Fixe Anzahl an threads
addsig2vol_3_mex(threads);
addsig2vol_3_orig(threads);

% Rechne blocked-AScans
count=blocksize; senderPos = 0.01.*rand(3,count); receiverPos = 0.01.*rand(3,count); IMAGE_STARTPOINT = [0,0,0]; IMAGE_RESOLUTION= 0.001; TimeInterval=1e-7;
DataLength=AscanLength;
Data=zeros(AscanLength,count);
Speed=1500+rand(1,count);
%Data(floor(DataLength.*rand(count-1,1))+1,1:count)=1;
Data = rand(AscanLength,count);
imagesum=zeros([x,x,x]);

[solution, buf]= addsig2vol_3_orig(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),imagesum);

for s=1:1:samples
    while 1
        tic;
        [n1, n2]= addsig2vol_3_mex(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),imagesum);
        timesBlocked(s)=toc;,  if timesBlocked(s)<10^8 break; end %%workaround for ugly times
    end
    errorArray(s) = abs(sum(sum(sum(n1-solution))));
end

for s=1:1:samples  
    while 1
        tic;
        [o1, o2]= addsig2vol_3_orig(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),imagesum);
        timesBlockedOrig(s)=toc;,  if timesBlockedOrig(s)<10^8 break; end %%workaround for ugly times
    end
end

%%%SAFT von ascan unblocked, averaging time measurements
figure; 
subplot(2,1,1);
hold on;
plot(1:1:samples,(blocksize*x.^3)./(timesBlocked), 'b'); 
plot(1:1:samples,(blocksize*x.^3)./(timesBlockedOrig), 'r');
legend("experimental", "original");
hold off;
subplot(2,1,2);
plot(1:1:samples,errorArray);
legend("Fehlerabweichung");