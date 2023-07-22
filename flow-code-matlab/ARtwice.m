close all
clear all

downhalf_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\DhalfR\';
% upsampled_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwice\';
% upsampled_BI_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwiceBI\';
% upsampled_BC_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwiceBC\';

upsampled_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwice2\';
upsampled_BI_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwiceBI2\';
upsampled_BC_path = 'C:\Users\23839\Desktop\RAFT-master\RAFT-master\uptwiceBC2\';


k = 50;
for i = 1:(k-1) 
     flo_file_name = strcat(downhalf_path, 'flow_', sprintf('%05d', i), '_', sprintf('%05d', i+1), '.flo');
     disp(['Processing ' flo_file_name]);
     flow_data = readFlowFile(flo_file_name);

    %% Duplication
    % Upsample the flow data
      upsampled_flow_data = upsample_twice(flow_data);
%  
%      % Save the upsampled_twice_duplication flow data 
      %upsampled_2 = strcat(upsampled_path, 'flow_', sprintf('%05d', i), '_', sprintf('%05d', i+1), '_upsampled.flo');
      upsampled_2 = strcat(upsampled_path, 'flow_', sprintf('%05d', i), '_', sprintf('%05d', i+1), '.flo');
      writeFlowFile(upsampled_flow_data, upsampled_2);

    %% Bilinear Interpolate
%        upsampled_BI_flow_data = bilineartwice(flow_data);
%    
%        % Save the upsampled_twice_duplication flow data 
%        %upsampled_BC = strcat(upsampled_BC_path, 'flow_', sprintf('%05d', i), '_', sprintf('%05d', i+1), '_BIupsampled.flo');
%        upsampled_BI = strcat(upsampled_BI_path, 'flow_', sprintf('%05d', i), '_', sprintf('%05d', i+1), '.flo');
%        writeFlowFile(upsampled_BI_flow_data, upsampled_BI);

    %% Bicubic Interplate
%        upsampled_BC_flow_data = Bcubicupsample(flow_data);
%    
%        % Save the upsampled_twice_duplication flow data 
%        %upsampled_BC = strcat(upsampled_BC_path, 'flow_', sprintf('%05d', i), '_', sprintf('%05d', i+1), '_BCupsampled.flo');
%        upsampled_BC = strcat(upsampled_BC_path, 'flow_', sprintf('%05d', i), '_', sprintf('%05d', i+1), '.flo');
%        writeFlowFile(upsampled_BC_flow_data, upsampled_BC);


end

