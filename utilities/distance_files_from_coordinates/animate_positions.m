# Copyright 2018 Oticon A/S
# SPDX-License-Identifier: Apache-2.0

function [] = animate_positions(device, outputfileprefix, SpeedRatio)

%generate distances and print to output files
N = length(device);

TimeResolution = 0.1; %50 ms time resolution (max error of 7cm at 10km/h)

supersettimes = device{1}.time;
for i = 2:N,
  supersettimes = union(supersettimes, device{i}.time);
end
supersettimes = union(supersettimes, min(supersettimes):TimeResolution:max(supersettimes));

x = zeros(N,numel(supersettimes));
y = zeros(N,numel(supersettimes));

for i = 1:N,
    x(i,:) = interpolate(device{i}.time, device{i}.x, supersettimes);
    y(i,:) = interpolate(device{i}.time, device{i}.y, supersettimes);
end


clf;
h = plot(x(:,1), y(:,1),'o');
axis([min(min(x)),max(max(x)),min(min(y)),max(max(y))])
axis equal
%ttt=tic; 

t1=tic;
for k = 2:length(supersettimes)
    set(h,'XData',x(:,k))
    set(h,'YData',y(:,k))
    title(sprintf('%s: %03.3f seconds',outputfileprefix, supersettimes(k)));
    CurrentExpectedTime = (supersettimes(k)-supersettimes(1))/SpeedRatio;
    ActualWait = toc(t1);
    
    pause(max( CurrentExpectedTime - ActualWait , 0));
    drawnow
end

%toc(ttt)
%(supersettimes(end)-supersettimes(1))/SpeedRatio
