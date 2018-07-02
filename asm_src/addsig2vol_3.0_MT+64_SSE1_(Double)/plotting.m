figure; 
hold on;


plot(voxel,(repmat(blocksize',[1 length(voxel)]).*voxel.^3)./mean(orig,3), '.-.');
plot(voxel,(repmat(blocksize',[1 length(voxel)]).*voxel.^3)./mean(stickTimes,3),'-');
plot(voxel,(repmat(blocksize',[1 length(voxel)]).*voxel.^3)./mean(interlaced,3),'-');

title('blocked');
legend('orginal', 'experimental', 'interlaced');



%a=2;
%plot(voxel,(40*a.*voxel.^3)./(mean(orig,3))(a,:), '.-.r');
%plot(voxel,(40*a.*voxel.^3)./(mean(stickTimes,3))(a,:), '-r');

%a=3;
%plot(voxel,(a.*voxel.^3)./(mean(orig,3))(a,:), '.-.g');
%plot(voxel,(a.*voxel.^3)./(mean(stickTimes,3))(a,:), '-g');
%
%a=4;
%plot(voxel,(a.*voxel.^3)./(mean(orig,3))(a,:), '.-.m');
%plot(voxel,(a.*voxel.^3)./(mean(stickTimes,3))(a,:), '-m');
%
%a=5;
%plot(voxel,(a.*voxel.^3)./(mean(orig,3))(a,:), '.-.k');
%plot(voxel,(a.*voxel.^3)./(mean(stickTimes,3))(a,:), '-k');
%
