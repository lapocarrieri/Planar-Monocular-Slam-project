
source "./ba_utils.m"


l=loadworld("dataset_for_one_person/world.dat");

K = camera.camera_matrix


%load('variables.mat',x_world)
xi=0.1;
num_iterations = 3;
damping = 0.1;
num_landmarks = 999;
num_poses = 200;


%for id = 1: 1000
%	if(x_world(id,3)<0 || x_world(id,3)>5)
%		x_world(id,:)= 0;
%	endif
%endfor

# landmarks in a matrix, one per column


# check if all the landmarks are visible is already done in functions.m where the landmark with less than 
# 2 measurements give as result   0 0 0; 

# compute associations on those landmarks

t=loadtrajectory("dataset_for_one_person/trajectoy.dat");
load("variables.mat","x_world");


XR_guess=zeros(4,4,num_poses);
for i= 1:1000
    XL_real(1:3,i)=[l(i).x_pose;l(i).y_pose;l(i).z_pose];
    
endfor
XL_guess=transpose(x_world);
pose_num = 1;

for (pose_num=1:num_poses)
	
    Xr=v2T(t(pose_num).odometryPose);
    XR_guess(:,:,pose_num)=Xr;
endfor;

#   Z:  the measurements (3xnum_measurements)
#   associations: 2xnum_measurements. 
#                 associations(:,k)=[p_idx,l_idx]' means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx

 

kernel_threshold = 10.0;

pose_dim = 6;
landmark_dim = 3;
  
chi_stats=zeros(1,num_iterations);
num_inliers=zeros(1,num_iterations);
    # size of the linear system
    Z = zeros(2,0);
system_size=pose_dim*num_poses+landmark_dim*num_landmarks;
for (iteration=1:num_iterations)
  printf("iteration number %d\n",iteration)
  H=zeros(system_size, system_size);
  b=zeros(system_size,1);
  chi_stats(iteration)=0;
      for pose = 0:199
        k=num2str(pose,'%05.f');
        pose=pose+1;
        T1 = inv(camera.cam_transform);
        T2 = inv([cos(t(pose).odometryPose(3)) -sin(t(pose).odometryPose(3)) 0 t(pose).odometryPose(1); sin(t(pose).odometryPose(3)) cos(t(pose).odometryPose(3))  0 t(pose).odometryPose(2); 0 0 1 0; 0 0 0 1]);
        T=T1*T2;%trasformation from world to camera
        filename = strcat("dataset_for_one_person/meas-", k,  '.dat');
        observations_t=loadmeas(filename);
        for j = 1: length(observations_t)

           z =  [observations_t(j).observation.x_pose ; observations_t(j).observation.y_pose];
           Xr = T;
           if(observations_t(j).observation.id == 0)
              continue;
           endif
           Xl = [l(observations_t(j).observation.id).x_pose;l(observations_t(j).observation.id).y_pose;l(observations_t(j).observation.id).z_pose];
          [is_behind, e,Jr,Jl] = errorAndJacobian(Xr, Xl, z, K);

          if (is_behind)
            continue;
          endif
          chi=e'*e;
          # apply the robust kernel as a weight to the linear system
          w = 1;
          if (chi>kernel_threshold)
            w = 1/kernel_threshold;
            chi=kernel_threshold;
          else
            num_inliers(iteration)++;
          endif;
          chi_stats(iteration)+=chi;
          omega = w * eye(2);
          Hrr = Jr' * omega * Jr;
          Hrl = Jr' * omega * Jl;
          Hll = Jl' * omega * Jl;
          br = Jr' * omega * e;
          bl = Jl' * omega * e;
    
          pose_matrix_index=poseMatrixIndex(pose, num_poses, num_landmarks);
          landmark_matrix_index=landmarkMatrixIndex(observations_t(j).observation.id , num_poses, num_landmarks);
    
          H(pose_matrix_index:pose_matrix_index+pose_dim-1,
      pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrr;
    
          H(pose_matrix_index:pose_matrix_index+pose_dim-1,
      landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hrl;
    
          H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
      landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=Hll;
    
          H(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,
      pose_matrix_index:pose_matrix_index+pose_dim-1)+=Hrl';
    
          b(pose_matrix_index:pose_matrix_index+pose_dim-1)+=br;
          b(landmark_matrix_index:landmark_matrix_index+landmark_dim-1)+=bl;
        endfor
      endfor
      size_H=size(H,1);
      rank_H=rank(H);
      H+=eye(system_size)*damping;
      dx=zeros(system_size,1);
      
      % we solve the linear system, blocking the first pose
      % this corresponds to "remove" from H and b the locks
      % of the 1st pose, while solving the system
  
      dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
      [XR_guess, XL_guess]=boxPlus(XR_guess,XL_guess,num_poses, num_landmarks, dx)

endfor

sum_sq_errors1 = 0;
sum_sq_errors2 = 0;
  for pose = 0:198
    k=num2str(pose,'%05.f');
    pose=pose+2;
    rel_T = XR_guess(:,:,pose)*XR_guess(:,:,pose-1);
    rel_GT = inv(v2T(t(pose).GroundTruthPose))* v2T(t(pose-1).GroundTruthPose);
%    - compute the relative motion rel_T = inv(T_0)*T_1
%- compute the relative gt motion rel_GT = inv(GT_0)*GT_1
    error_T = inv(rel_T)*rel_GT;
%- compute the SE(2) error error_T = inv(rel_T)*rel_GT
%- rotation part:
    Angle(pose) = atan2(error_T(2, 1), error_T(1, 1));
%- translation part
%	- compute the RMSE on translations error_T(1:3, 4) 
    translations_error= error_T(1:3, 4);
    
    sum_sq_errors1 = sum_sq_errors1 + norm(translations_error,'fro')^2;
    sum_sq_errors2 = sum_sq_errors2 + norm(error_T, 'fro')^2;
  endfor
rmseT = sqrt(sum_sq_errors2 / (4 * pose));
rmse_translations = sqrt(sum_sq_errors1 / (pose*size(translations_error,1)));


% Display the result
fprintf('The RMSE  for the homogeneus matrixs is %f\n', rmseT);


% Display the result

fprintf('The RMSE for the translational vector is %f\n', rmse_translations);

figure;
x1=transpose(XL_guess(1,:));
y1=transpose(XL_guess(2,:));
z1=transpose(XL_guess(3,:));
x2=transpose(XL_real(1,:));
y2 = transpose(XL_real(2,:));
z2=transpose(XL_real(3,:));
scatter3(x1,y1 ,z1, 'b', 'filled');
scatter3(x2,y2 ,z2, 'r', 'filled');

xlabel('X-axis');
ylabel('Y-axis');
zlabel('Z-axis');

for i=1:length(x1)
  text(x1(i), y1(i), z1(i), num2str(i), 'FontSize', 1, 'FontWeight', 'bold');
  text(x2(i), y2(i), z2(i), num2str(i), 'FontSize', 1, 'FontWeight', 'bold');
end
axis equal;
title('Landmarks');

saveas(gcf, 'my_figure.png');
