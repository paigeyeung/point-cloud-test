
filename = '~/Documents/MATLAB/highway_data/HighwayLong3.csv_10.10.10.143_17038.csv';
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
%     pc = pctransform(pc,tform);     
    view(player,pc);     
end