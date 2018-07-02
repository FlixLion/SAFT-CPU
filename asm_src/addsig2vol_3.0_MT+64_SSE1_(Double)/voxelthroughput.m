clear all()
close all
% Bildgröße x-Achse (Immer durch 4 teilbar. Vollständiges Bild kubisch)
voxelgenerator=unique(round(2.^(2:0.08:7)))
voxel=voxelgenerator.+(4.-mod(voxelgenerator,4));
voxel=unique(voxel)
% Anzahl der AScans
blocksize=[1 30];
% Anzahl der Messungen pro Konfiguration
samples = 200;
AscanLength = 3000;
% timetable
timesBlocked = zeros(length(blocksize),length(voxel),1:1:samples);
timesBlockedOrig = zeros(length(blocksize),length(voxel),1:1:samples);
timesBlockedInterlaced = zeros(length(blocksize),length(voxel),1:1:samples);

% Fixe Anzahl an threads
addsig2vol_3_interlaced(4);
addsig2vol_3_orig(4);
addsig2vol_3_sticky(4);


% Rechne blocked-AScans (interlaced)
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
                [bild, buffer]= addsig2vol_3_interlaced(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),imagesum);
                timesBlockedInterlaced(i,j,s)=toc,  if timesBlockedInterlaced(i,j,s)<10^8 break; end %%workaround for ugly times
            end
        end
    end
end


% Rechne blocked-AScans (orig)
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
                timesBlockedOrig(i,j,s)=toc,  if timesBlockedOrig(i,j,s)<10^8 break; end %%workaround for ugly times
            end
        end
    end
end

% Rechne blocked-AScans (sticky)
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
                [bild, buffer]= addsig2vol_3_sticky(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),imagesum);
                timesBlocked(i,j,s)=toc,  if timesBlocked(i,j,s)<10^8 break; end %%workaround for ugly times
            end
        end
    end
end




stickTimes = timesBlocked;
orig = timesBlockedOrig;
interlaced = timesBlockedInterlaced;

%%%SAFT von ascan unblocked, averaging time measurements
%figure; imagesc(voxel,blocksize,(repmat(blocksize',[1 length(voxel)]).*voxel.^3)./median(timesBlocked,3)); colorbar;title('blocked')
%figure; plot(voxel,(repmat(blocksize',[1 length(voxel)]).*voxel.^3)./median(timesBlocked,3));title('blocked')
%figure; plot(voxel,(repmat(blocksize',[1 length(voxel)]).*voxel.^3)./median(timesBlockedOrig,3));title('blocked')


%B = median(timesBlocked,3);
%figure; plot(voxel,B');title('blocked')

%%%SAFT von ascan unblocked, averaging time measurements
%figure; imagesc(voxel,1:length(blocksize),voxel.^3./mean((times(1:length(blocksize),:,2))./repmat((1:length(blocksize))',[1 length(voxel)]))  ); colorbar, title('unblocked')
%figure; plot(voxel,(voxel.^3)./mean(timesUnblocked,2)');title('unblocked')

% Abgeschnittener Plot
%a= 30;
%figure; plot(voxel(1:a),(repmat(blocksize',[1 length(voxel(1:a))]).*voxel(1:a).^3)./mean(timesBlocked,3)(:,1:a));title('blocked')
%figure; imagesc(voxel(1:a),blocksize,(repmat(blocksize',[1 length(voxel(1:a))]).*voxel(1:a).^3)./mean(timesBlocked,3)(:,1:a)); colorbar;title('blocked')

%B = mean(timesBlocked,3)(:,1:a);
%figure; plot(voxel(1:a),B');title('blocked')