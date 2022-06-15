
filename = 'C:\Users\pyeung\Downloads\HighwayLong3.csv_10.10.10.143_17038.csv';
T = readtable(filename, 'PreserveVariableNames', true);

x = T.("x (m)")
y = T.("y (m)")
z = T.("z (m)")
intens = T.intensity
xyz = [x, y, z];

writematrix(xyz,'coordinates.csv')
pc = pointCloud(xyz, 'Intensity', intens);

lower = min([pc.XLimits pc.YLimits]);
upper = max([pc.XLimits pc.YLimits]);
  
xlimits = [lower upper];
ylimits = [lower upper];
zlimits = pc.ZLimits;

% player = pcplayer(xlimits,ylimits,zlimits);
% 
% xlabel(player.Axes,'X (m)');
% ylabel(player.Axes,'Y (m)');
% zlabel(player.Axes,'Z (m)');
% 
% for i = 1:360      
%     view(player,pc);     
% end

normals = pcnormals(pc);
x1 = pc.Location(1:100:end,1);
y1 = pc.Location(1:100:end,2);
z1 = pc.Location(1:100:end,3);
u = normals(1:100:end,1);
v = normals(1:100:end,2);
w = normals(1:100:end,3);

writematrix(normals,'normals.csv');

figure
pcshow(pc)
title('Estimated Normals of Point Cloud')
hold on
quiver3(x1,y1,z1,u,v,w);
hold off

%% Divide point cloud into cells

index = 1;
incx = 50;
incy = 50;
incz = 5;

r = zeros(size(x, 1),1);

for i = min(pc.XLimits)-1:incx:max(pc.XLimits)
    for j = min(pc.YLimits)-1:incy:max(pc.YLimits)
        for k = min(pc.ZLimits)-1:incz:max(pc.ZLimits)
            indices = x>i & x<=i+incx & y>j & y<=j+incy & z>k & z<=k+incz;
            %size(indices(indices==1))
            r(indices) = index;
            index = index + 1;
        end
    end
end

writematrix(r,'indexes.csv');
cmap = jet(index);

%% Create new point cloud
newPoints = zeros(index,3);

for i = 1:index
    if size(x(r==i),1) > 0
        newPoints(i,1) = mean(x(r==i));
        newPoints(i,2) = mean(y(r==i));
        newPoints(i,3) = mean(z(r==i));
    end
end
figure(), scatter3(x,y,z,10,cmap(r,:),'filled')

newPc = pointCloud(newPoints);

lower = min([newPc.XLimits newPc.YLimits]);
upper = max([newPc.XLimits newPc.YLimits]);
  
xlimits = [lower upper];
ylimits = [lower upper];
zlimits = newPc.ZLimits;

player = pcplayer(xlimits,ylimits,zlimits);

xlabel(player.Axes,'X (m)');
ylabel(player.Axes,'Y (m)');
zlabel(player.Axes,'Z (m)');

for i = 1:360      
    view(player,newPc);     
end

newNormals = pcnormals(newPc);
x1 = newPc.Location(1:100:end,1);
y1 = newPc.Location(1:100:end,2);
z1 = newPc.Location(1:100:end,3);
u = newNormals(1:100:end,1);
v = newNormals(1:100:end,2);
w = newNormals(1:100:end,3);

figure
pcshow(newPc)
title('Estimated Normals of New Point Cloud')
hold on
quiver3(x1,y1,z1,u,v,w);
hold off

%dividednormals = horzcat(normals,r);

%% New normals
newNormals = zeros(index, 3);

for i = 1:index
    if size(x(r==i),1) > 0
        newNormals(i,1) = 100*mean(normals(r==i,1));
        newNormals(i,2) = 100*mean(normals(r==i,2));
        newNormals(i,3) = 100*mean(normals(r==i,3));
    end
end

x1 = newPc.Location(1:1:end,1);
y1 = newPc.Location(1:1:end,2);
z1 = newPc.Location(1:1:end,3);
u = newNormals(1:1:end,1);
v = newNormals(1:1:end,2);
w = newNormals(1:1:end,3);

writematrix(newNormals, 'newNormals.csv')

figure
pcshow(newPc)
title('Estimated Normals of New Point Cloud (2)')
hold on
quiver3(x1,y1,z1,u,v,w);
hold off
