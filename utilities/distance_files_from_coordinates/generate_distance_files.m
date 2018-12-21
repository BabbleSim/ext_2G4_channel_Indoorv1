# Copyright 2018 Oticon A/S
# SPDX-License-Identifier: Apache-2.0

function [] = generate_distance_files(device, outputfileprefix, OutputDir, RelativePath)
%From a set of device times and coordinates, generate the distances files used by
%the 2G4Indoorv1 channel.
%Note that the device{}.times are in seconds and device{}.x and .y
%coordinates are in meters

%generate distances and print to output files
N = length(device);

MatrixF = fopen([OutputDir outputfileprefix '.matrix'],'w');
fprintf(MatrixF,'#<Txnbr> <Rxnbr> : {<distance>|"<distance_file_name>"}\n');

Acceptable_error = 0.05; %5 cm (we dont care about errors smaller than 5 cm)
% previous script had an acceptable error of 0.001 == 1mm

TimeResolution = 0.05; %50 ms time resolution (max error of 7cm at 10km/h)

plotdistance = 0; %debug

for tx = 1:N,
  for rx = 1:tx-1,
    
    assert(issorted(device{tx}.time) && issorted(device{rx}.time),'Times must be in ascending order');
    assert( all( device{tx}.time == unique(device{tx}.time))...
           && all( device{rx}.time == unique(device{rx}.time)),...
           'Times cannot repeat (a device cannot be in 2 places at the same time)');
         
    supersettimes = union(device{tx}.time, device{rx}.time);
    supersettimes = union(supersettimes, min(supersettimes):TimeResolution:max(supersettimes));
    
    c_tx.x = interpolate(device{tx}.time, device{tx}.x, supersettimes);
    c_tx.y = interpolate(device{tx}.time, device{tx}.y, supersettimes);
    
    c_rx.x = interpolate(device{rx}.time, device{rx}.x, supersettimes);
    c_rx.y = interpolate(device{rx}.time, device{rx}.y, supersettimes);
    
    distance = sqrt((c_tx.x - c_rx.x).^2 + (c_tx.y - c_rx.y).^2);
    
    %remove consecutive distances which are linear interpolations of each surrounding ones
    keep = true(size(distance));
    last_kept=1;
    
    for i = 2:numel(distance)-1,
      Ok_to_remove = 1;
      j = (last_kept+1):i; %for all points which would get interpolated if we remove i
      interp =   ( distance(i+1) - distance(last_kept) ) / ( supersettimes(i+1) - supersettimes(last_kept) ) ...
               * ( supersettimes(j) - supersettimes(last_kept) ) ...
               + distance(last_kept);
       
      Error = abs(distance(j) - interp);
      if ( max(Error) >= Acceptable_error ),
        Ok_to_remove = 0;
      end 
      
      if Ok_to_remove,
        %if (tx == 3) && (rx == 2),
        %  fprintf('Removed point %i, time %f, distance %f~=%f (error = %i)\n',i, supersettimes(i), distance(i), interp,Error);
        %end
        keep(i) = false;
      else
        last_kept=i;
      end
    end %for i
    
    %remove the first distance if it is equal to the next
    if ( length(distance) > 1 ) && (keep(2) == true) && ( distance(1) == distance(2) ),
      keep(1) = false;
    end
    %remove the last distance if it is equal to the previous
    if ( length(distance) > 1 ) && (keep(end-1) == true) && ( distance(end) == distance(end-1) ),
      keep(end) = false;
    end

    if plotdistance,
      figure(); clf;
      plot(supersettimes, distance);
      hold on;
    end
    distance = distance(keep);
    supersettimes = supersettimes(keep);
    if plotdistance,
      plot(supersettimes, distance,'or');
      title([num2str(tx) ' to ' num2str(rx) ]);
      pause;
    end
    
    if numel(distance)==1,
      fprintf(MatrixF,'%i %i : %i\n',tx-1, rx-1, distance);
    else
      filenameopen = [OutputDir outputfileprefix '.' num2str(tx-1) '_' num2str(rx-1) '.dist'];
      filenameprint = [RelativePath outputfileprefix '.' num2str(tx-1) '_' num2str(rx-1) '.dist'];
      fprintf(MatrixF,'%i %i : "%s"\n',tx-1, rx-1, filenameprint);
      tmpfile = fopen(filenameopen,'w');
      fprintf(tmpfile,'#<time> <distance>\n');
      fprintf(tmpfile,'%i %e\n',[round(supersettimes'*1e6) distance']');
      fclose(tmpfile);
    end
  end
end

fclose(MatrixF);
