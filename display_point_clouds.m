
filename = 'C:\Users\pyeung\Downloads\HighwayLong3.csv_10.10.10.143_17038.csv';
T = readtable(filename, 'PreserveVariableNames', true);

x = T.("x (m)")
y = T.("y (m)")
z = T.("z (m)")
intens = T.intensity
xyz = [x, y, z];
pc = pointCloud(xyz, 'Intensity', intens);

lower = min([pc.XLimits pc.YLimits]);
upper = max([pc.XLimits pc.YLimits]);
  
xlimits = [lower upper];
ylimits = [lower upper];
zlimits = pc.ZLimits;

player = pcplayer(xlimits,ylimits,zlimits);

xlabel(player.Axes,'X (m)');
ylabel(player.Axes,'Y (m)');
zlabel(player.Axes,'Z (m)');

for i = 1:360      
    view(player,pc);     
end

normals = pcnormals(pc);
x = pc.Location(1:10:end,1);
y = pc.Location(1:10:end,2);
z = pc.Location(1:10:end,3);
u = normals(1:10:end,1);
v = normals(1:10:end,2);
w = normals(1:10:end,3);

figure
pcshow(pc)
title('Estimated Normals of Point Cloud')
hold on
quiver3(x,y,z,u,v,w);
hold off

index = 2;
incx = 50;
incy = 50;
incz = 10;

%r = zeros(350/inc,175/inc,20/inc);

r = ones(size(normals, 1));

for i = 0:350:incx
    for j = -100:175:incy
        for k = -5:15:incz
            r(x>i & x<=i+incx & y>j & y<=y+incy & z>k & z<=k+incz) = index;
            index = index + 1;
        end
    end
end

writematrix(r,'indexes.csv');

cmap = jet(20000);
figure(), scatter3(x,y,z,10,cmap(r,:),'filled')

%dividednormals = horzcat(normals)