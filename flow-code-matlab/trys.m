%% try some examples



%% numl ---Calculate the absolute error (AE) between the forward and ground truth flows
%      ae_forward = abs(truth - festflo_forward);
%     % Calculate the mean absolute error (MAE) over the flow image
%      mae_forward = sum(ae_forward(:)) / numel(ae_forward);
%      mae_values(i) = mae_forward;
%figure(10)
% yyaxis left
% plot(frame_numbers, mae_values, '-o')
% ylabel('MAE value')
% yyaxis right
% plot(frame_numbers, mse_values, '-x')
% xlabel('Frame number')
% ylabel('Motion compensation MSE')
%% 
% est_flo = readFlowFile('flow_00001_00002.flo');
% est_flo1 = readFlowFile('flow_00015_00016.flo');
% est_flo2 = readFlowFile('flow_00014_00015.flo');
% truth14 = readFlowFile('frame_0014.flo');
% test14 = abs(est_flo -  truth14);
% mae_14 = mean(test14(:));
% fprintf('Mean Absolute Error: %f\n', mae_14);
% 
% test141 = abs(est_flo1 -  truth14);
% mae_141 = mean(test141(:));
% fprintf('Mean Absolute Error: %f\n', mae_141);

% est_flo = readFlowFile('flow_00001_00002.flo');
% truth15 = readFlowFile('frame_0014.flo'); %groundtruth of frame15
% truth16 = readFlowFile('frame_0015.flo'); %groundtruth of frame16
% truth17 = readFlowFile('frame_0016.flo'); %groundtruth of frame17
% 
% % Upsample estimated flow to the same size of ground truth flow
% est_flow_upsample = imresize(est_flo, [size(truth17, 1), size(truth17, 2)]);
% 
%
%% MAE
% test15 = abs(est_flo -  truth15);
% mae_15 = mean(test15(:));
% fprintf('Mean Absolute Error: %f\n', mae_15);
% test16 = abs(est_flo -  truth16);
% mae_16 = mean(test16(:));
% fprintf('Mean Absolute Error: %f\n', mae_16);
% test17 = abs(est_flo -  truth17);
% mae_17 = mean(test17(:));
% fprintf('Mean Absolute Error: %f\n', mae_17);


%     % Calculate the magnitude of the ground truth flow
%     truth_magnitude = sqrt(truth(:,:,1).^2 + truth(:,:,2).^2);
%     truth_magnitudes(i) = mean(truth_magnitude(:));
%plot(frame_numbers, truth_magnitudes, 'LineWidth', 2);
% xlabel('Frame Number');
% ylabel('Ground Truth Magnitude');
% title('Ground Truth Magnitude vs. Frame Number');
% grid on;





%% EPE 
% Compute endpoint error (EPE) between ground truth and estimated flow
% epe15 = sqrt(sum((truth15 - est_flow_upsample).^2, 3));
% mean_epe15 = mean(epe15(:));
% fprintf('Mean Endpoint Error: %f\n', mean_epe15);
% 
% epe16 = sqrt(sum((truth16 - est_flow_upsample).^2, 3));
% mean_epe16 = mean(epe16(:));
% fprintf('Mean Endpoint Error: %f\n', mean_epe16);
% 
% epe17 = sqrt(sum((truth17 - est_flow_upsample).^2, 3));
% mean_epe17 = mean(epe17(:));
% fprintf('Mean Endpoint Error: %f\n', mean_epe17);
% 
% test15(1,1,:);
% %test 15 ans(:,:,1) = 0.0610, ans(:,:,2) = 0.0152
% test16(1,1,:);
% %test 16 ans(:,:,1) = 0.0498, ans(:,:,2) = 0.0138
% test17(1,1,:);
% %test 17 ans(:,:,1) = 0.0383, ans(:,:,2) = 0.0124
% 
% out17 = flowToColor(test17);
% figure(5)
% imshow(out17);
% out16 = flowToColor(test16);
% figure(6)
% imshow(out16);

%% (no downsample) Visualize flow fields as arrows overlaid on the image
% figure(9); imshow(im); hold on;
% quiver(truth17(:,:,1), truth17(:,:,2), 'y');
% title('Ground Truth Flow');
% figure(10); imshow(im); hold on;
% quiver(est_flow_upsample(:,:,1), est_flow_upsample(:,:,2), 'r');
% title('Estimated Flow');



% figure(5)
% imshow(flowToColor(truth16))
% figure(6)
% imshow(flowToColor(truth17))
% figure(7)
% imshow(flowToColor(truth15))
% figure(8)
% imshow(flowToColor(est_flo))