clear all()
close all

count=5;
senderPos = 0.01.*rand(3,count);
receiverPos = 0.01.*rand(3,count);
IMAGE_STARTPOINT = [0,0,0];
IMAGE_RESOLUTION= 0.001;
Speed=1500;
TimeInterval=1e-7;
DataLength=1000;
Data=floor(rand(DataLength,count).*2);
x=120;


 
addsig2vol_3_mex(3);
tic;
[n1,n2,n3,n4,n5,n6]=addsig2vol_3_mex(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),zeros([x,x,x]));
t1 = toc
%n1: Rekonstruiertes Bild
%n2: Scanbuffer
%n3: Zeitstart für Tasks
%n4: Zeitende für Tasks
%n5: Tasknummern
%n6: Threadnummern

% Plotte Reihenfolge der tasks 
figure;
plot(n3(1:count,:), n5(1:count,:), 'x'); hold on; %starts
%plot(n4(1:count,:), n5(1:count,:), 'x'); hold on; %ends
xlabel('time[milisec]');
ylabel('tasknumber[nacheinanderliegend im Speicher]');

% Trennlinie pro fertigem AScan zeichnen
maxTask = max(n5(1,:));
for i=1:1:count
  %plot([n4(count,maxTask),n4(count,maxTask)], [0,maxTask]);
end
%finishedTask = max(n4(count,:));
% Linie für Programmstart
plot([0,0], [0,maxTask], 'r');
% Linie für Threads-Ende
%plot([finishedTask,finishedTask], [0,maxTask]);
% Linie für Programmende
%plot([t1*1000,t1*1000], [0,maxTask], 'r');


% thread nummer labels
for j=1:1:maxTask+1
  for i=1:1:count
    b = num2str(n6(i,j));
    c = cellstr(b);
    dx = 0.01; dy = 0.1; % displacement so the text does not overlay the data points
    text(n3(i,j)+dx, n5(i,j)+dy, c);
  end
end

hold off;