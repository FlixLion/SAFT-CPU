clear all()
close all
% Bildgröße x-Achse (Immer durch 4 teilbar. Vollständiges Bild kubisch)
voxelgenerator=unique(round(2.^(2:0.1:7)))
voxel=voxelgenerator.+(4.-mod(voxelgenerator,4));
voxel=unique(voxel)
voxel=160;
% Anzahl der AScans
blocksize=1;
% Anzahl der Messungen pro Konfiguration
samples = 80;
AscanLength = 3000;
% timetable
timesBlocked = zeros(length(blocksize),length(voxel),1:1:samples);
timesBlockedOrig = zeros(length(blocksize),length(voxel),1:1:samples);

% Fixe Anzahl an threads
addsig2vol_3_mex(4);

% Rechne blocked-AScans
for i=1:length(blocksize)
    for j=1:length(voxel)
        count=blocksize(i); senderPos = 0.01.*rand(3,count); receiverPos = 0.01.*rand(3,count); IMAGE_STARTPOINT = [0,0,0]; IMAGE_RESOLUTION= 0.001; TimeInterval=1e-7;
        DataLength=AscanLength;
        Data=zeros(AscanLength,count);
        Speed=1500+rand(1,count);
        %Data(floor(DataLength.*rand(count-1,1))+1,1:count)=1;
        x=voxel(j);
        Data = rand(AscanLength,count);
        imagesum=zeros([x,x,x]);

        for s=1:1:samples
            while 1
                tic;
                [bild, buffer]= addsig2vol_3_mex(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),imagesum);
                timesBlocked(s)=toc;,  if timesBlocked(s)<10^8 break; end %%workaround for ugly times
            end
        end
    end
end

addsig2vol_3_orig(4);
% Rechne blocked-AScans (original Algorithmus)
for i=1:length(blocksize)
    for j=1:length(voxel)
        count=blocksize(i); senderPos = 0.01.*rand(3,count); receiverPos = 0.01.*rand(3,count); IMAGE_STARTPOINT = [0,0,0]; IMAGE_RESOLUTION= 0.001; TimeInterval=1e-7;
        DataLength=AscanLength;
        Data=zeros(AscanLength,count);
        Speed=1500+rand(1,count);
        %Data(floor(DataLength.*rand(count-1,1))+1,1:count)=1;
        x=voxel(j);
        Data = rand(AscanLength,count);
        imagesum=zeros([x,x,x]);

        for s=1:1:samples
            while 1
                tic;
                [bild, buffer]= addsig2vol_3_orig(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),imagesum);
                timesBlockedOrig(s)=toc;,  if timesBlockedOrig(s)<10^8 break; end %%workaround for ugly times
            end
        end
    end
end

%%%SAFT von ascan unblocked, averaging time measurements
figure; 
hold on;
plot(1:1:samples,(blocksize*voxel.^3)./(timesBlocked));
plot(1:1:samples,(blocksize*voxel.^3)./(timesBlockedOrig));
legend("experimental", "original");
