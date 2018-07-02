
clear();

% parameters
x = 57;
thread =  4;
multipleAscans = 1000;


count=1;
senderPos = 0.01.*rand(3,count);
receiverPos = 0.01.*rand(3,count);
IMAGE_STARTPOINT = [0,0,0];
IMAGE_RESOLUTION= 0.001;
Speed=1500;
TimeInterval=1e-7;
DataLength=1000;
Data=floor(rand(DataLength,count).*2);

addsig2vol_3_mex(thread);
addsig2vol_3_orig(thread);


tic;
[n1,n2]=addsig2vol_3_mex(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),zeros([x,x,x]));
t1 = toc
tic;
[o1,o2]=addsig2vol_3_orig(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),zeros([x,x,x]));
t2 = toc

a=n1-o1;
disp("1 Scan--------------")
diff_1 = (sum(sum(sum((n1-o1)))))/(x^3)
disp("--------------------")



count=multipleAscans;
senderPos = 0.01.*rand(3,count);
receiverPos = 0.01.*rand(3,count);
Data=floor(rand(DataLength,count).*2);

tic;
[m1,m2]=addsig2vol_3_mex(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),zeros([x,x,x]));
t1 = toc
tic;
[p1,p2]=addsig2vol_3_orig(Data,single(IMAGE_STARTPOINT),single(receiverPos),single(senderPos),single(Speed),single(IMAGE_RESOLUTION),single(TimeInterval),uint32([x,x,x]),zeros([x,x,x]));
t2 = toc

b=m1-p1;
disp("Mehrere Scans--------------")
diff_n = (sum(sum(sum((m1-p1)))))/(x^3)
disp("---------------------------")