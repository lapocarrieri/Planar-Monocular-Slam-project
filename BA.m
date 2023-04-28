
source "./ba_utils.m"
# This is the Bundle Adjustment file, it loads x_world that contains the initialization
# for the landmarks done in Triangulation.m

l=loadworld("dataset_for_one_person/world.dat");
# check if all the landmarks are visible is already done in functions.m where the landmark with less than 
# 2 measurements give as result   0 0 0; 
# compute associations on those landmarks


t=loadtrajectory("dataset_for_one_person/trajectoy.dat");
load("variables2.mat","x_world");
x_world1=x_world;
load("variables4.mat","x_world");
x_world2=x_world;
load("variables5.mat","x_world");
x_world3=x_world;


x_world=(x_world1+x_world2+x_world3)/3;


# preprocessing over x_world when the datas are too strange
for i=1:length(x_world(:,1))
  if abs(x_world(i,1))>15 || abs(x_world(i,2))>15 || abs(x_world(i,3))>3
    x_world(i,:)=[0 0 0 ];
    a=3
  endif
endfor

K = camera.camera_matrix;

#Hyperparameters that can be changed
xi=0.1;
num_iterations = 2
damping = 0.1;
kernel_threshold = 10.0;
#parameters and initializations
num_landmarks = 999;
num_poses = 200;
sum=0;
pose_dim = 6;
landmark_dim = 3;
chi_stats=zeros(1,num_iterations);
num_inliers=zeros(1,num_iterations);
Z = zeros(2,0);
system_size=pose_dim*num_poses+landmark_dim*num_landmarks;
e1=0;
XR_guess=zeros(4,4,num_poses);
XL_guess=transpose(x_world);


size(XL_guess);
pose_num = 1;


#In XR_guess there are all the GroundTruthPose that are necessary for the final testing 
for (pose_num=1:num_poses)
    Xr_OdometryPose=v2T(t(pose_num).odometryPose);
    XR_guess(:,:,pose_num)=Xr_OdometryPose;
endfor

# in l are stored all the real landmarks always to be tested at the end
for i= 1:999
    XL_real(1:3,i)=[l(i+1).x_pose;l(i+1).y_pose;l(i+1).z_pose];
endfor



fid=1;

#Here start the Bundle Adjustment, the BA is done num:iterations times and at every iteration the Xr_guess and the XL_guess are updated

for (iteration=1:num_iterations)
  printf("iteration number %d\n",iteration)
  H=zeros(system_size, system_size);
  b=zeros(system_size,1);
  
  chi_stats(iteration)=0;
      for pose = 1:200
        fclose('all');
        Kpose=pose-1;
        k=num2str(Kpose,'%05.f');
        %pose=pose+1;
        T1 = inv(camera.cam_transform); %transform from robot to camera frame

        % position of the robot in respect to the world that is inverted in order to transform from world frame to robot frame
        %T2 = inv([cos(t(pose).GroundTruthPose(3)) -sin(t(pose).GroundTruthPose(3)) 0 t(pose).GroundTruthPose(1); sin(t(pose).GroundTruthPose(3)) cos(t(pose).GroundTruthPose(3))  0 t(pose).GroundTruthPose(2); 0 0 1 0; 0 0 0 1]);
        T2 = inv([cos(t(pose).odometryPose(3)) -sin(t(pose).odometryPose(3)) 0 t(pose).odometryPose(1); sin(t(pose).odometryPose(3)) cos(t(pose).odometryPose(3))  0 t(pose).odometryPose(2); 0 0 1 0; 0 0 0 1]);
        
        T=T1*T2;%trasformation from world to camera

        filename = strcat("dataset_for_one_person/meas-", k,  '.dat');
        
        
        [observations_t,fid]=loadmeas(filename);
        
        
        for j = 1: length(observations_t)
          
           z =  [observations_t(j).observation.x_pose ; observations_t(j).observation.y_pose];
           Xr = T;
           
           if(observations_t(j).observation.id == 0);
            continue;
          endif  
          %Xl = [l(observations_t(j).observation.id+1).x_pose,l(observations_t(j).observation.id+1).y_pose,l(observations_t(j).observation.id+1).z_pose];
          Xl = [XL_guess(1,observations_t(j).observation.id),XL_guess(2,observations_t(j).observation.id),XL_guess(3,observations_t(j).observation.id)];
          Xl=transpose(Xl);
          
          %fclose(fid)
       
          
          
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
          
          landmark_matrix_index=landmarkMatrixIndex(observations_t(j).observation.id, num_poses, num_landmarks);
          
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
      [XR_guess, XL_guess]=boxPlus(XR_guess,XL_guess,num_poses, num_landmarks, dx);
endfor


sum_sq_errors1 = 0;
sum_sq_errors2 = 0;


printf("Now Evaluation")

  for pose = 2:num_poses
    
    
    rel_T = inv(XR_guess(:,:,pose))*XR_guess(:,:,pose-1);
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
    
    sum_sq_errors1 = sum_sq_errors1 + norm(translations_error,'fro');
    sum_sq_errors2 = sum_sq_errors2 + norm(error_T, 'fro');
  endfor
rmseT = sqrt(sum_sq_errors2 / (4 * pose));
rmse_translations = sqrt(sum_sq_errors1 / (pose*size(translations_error,1)));

sum_sq_errors1=0;
sum_sq_errors2=0;
sum_sq_errors3=0;

#RMSE of the landmarks(Map)

for i=1:998
    ErrorLandmarks(:,i)=XL_guess(:,i)-XL_real(:,i);
endfor
for id = 1:num_landmarks-1
    sum_sq_errors1 = sum_sq_errors1 + ErrorLandmarks(1,id)^2;
    sum_sq_errors2 = sum_sq_errors2 + ErrorLandmarks(2,id)^2;
    sum_sq_errors3 = sum_sq_errors3 + ErrorLandmarks(3,id)^2;
endfor
rmseLandmarks1 = sqrt(sum_sq_errors1/ num_landmarks);
rmseLandmarks2 = sqrt(sum_sq_errors2/ num_landmarks);
rmseLandmarks3 = sqrt(sum_sq_errors3/ num_landmarks);
rmseLandmarks = [rmseLandmarks1, rmseLandmarks2, rmseLandmarks3];
%rmse=rms(ErrorLandmarks)

% Display the result
fprintf('The RMSE  for the homogeneus matrixs is %f\n', rmseT);


% Display the result

fprintf('The RMSE for the translational vector is %f\n', rmse_translations);
fprintf('The RMSE  for the map so considering every landmarks in the system is %f\n', rmseLandmarks);


save('Results10.mat','XL_guess','XL_real','XR_guess','t')




