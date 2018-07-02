clear all()
close all
% Bildgröße x-Achse (Immer durch 4 teilbar. Vollständiges Bild kubisch)
%voxelgenerator=unique(round(2.^(2:0.1:7)))
%voxel=voxelgenerator.+(4.-mod(voxelgenerator,4));
% Anzahl der AScans
blocksize=1;
% Anzahl der Messungen pro Konfiguration
samples = 5;
AscanLength = 3000;
% timetable

% Fixe Anzahl an threads
addsig2vol_3_mex(3);

%Insgesamt: 10 MB
%Eine Scheibe: 10kB
x=10;
y=125;
z=3000;
vol = x*y*z;

% segmentation 
segraw = 1:10:100;
segment=segraw.*x.*y; %Z Schicht

timesBlocked = zeros(length(segment),1:1:samples);
 
% Vorlauf

% Rechne blocked-AScans

for l=1:length(segment)
  count=blocksize; senderPos = 0.01.*rand(3,count); receiverPos = 0.01.*rand(3,count); IMAGE_STARTPOINT = [0,0,0]; IMAGE_RESOLUTION= 0.001; TimeInterval=1e-7;
  DataLength=AscanLength;
  Data=zeros(AscanLength,count);
  Speed = 1500;
  %Data(floor(DataLength.*rand(count-1,1))+1,1:count)=1;
  Data = rand(AscanLength,count);
  imagesum=zeros([x,y,z]);
  % set segmentation
  addsig2vol_3_mex(1, segment(l));
      
  for s=1:1:samples
      while 1
          tic;
          [bild, buffer]= addsig2vol_3_mex(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),imagesum);
          timesBlocked(l,s)=toc;,  if timesBlocked(l,s)<10^8 break; end %%workaround for ugly times
      end
  end
end


figure; hold on;
plot(segraw,blocksize*vol./mean(timesBlocked,2), 'k');
plot(segraw,blocksize*vol./median(timesBlocked,2), 'b');
plot(segraw,blocksize*vol./max(timesBlocked,[],2), 'r');
plot(segraw,blocksize*vol./min(timesBlocked,[],2), 'g');
title('Segmentierung Z-Achse');
legend('mean', 'median', 'worst', 'best')
