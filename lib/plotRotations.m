function plotRotations(R_gt_list, R_est_list)

origin = zeros(3,1);

%% plot
figure
for i=1:size(R_est_list,3)
   my_ref_frame(R_est_list(:,:,i),origin,['r','g','b'],0.9,0.01,0.8) 
end

%% plot 
for i=1:size(R_gt_list,3)
   my_ref_frame(R_gt_list(:,:,i),origin,['r','g','b'],0.75,0.05,0.5)  
end
set(gca,'fontsize', 18)