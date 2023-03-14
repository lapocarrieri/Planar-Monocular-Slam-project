
source "./ba_utils.m"

t=loadtrajectory("/home/lapo/ProbRobLapo/dataset_for_one_person/trajectoy.dat");


%SLAM COMINCIA QUI
for c = 2:200
    transitions(c-1).v=t(c).odometryPose-t(c-1).odometryPose;
    transitions(c-1).id_from = t(c-1).id;
    transitions(c-1).id_to = t(c).id;
    
end
#set initial pose at the origin - we don't know the map and neither our location
mu = [0;  #x coordinate
      0;  #y coordinate
      0;
      0]; #orientation theta (yaw angle)
printf("initial pose: [%f, %f, %f, %f]\n", mu(1), mu(2), mu(3),mu(4));
h=0;
#initialize covariance: high value means high uncertainty
sigma = eye(3);

#bookkeeping: to and from mapping between robot pose (x,y, theta) and landmark indices (i)
#all mappings are initialized with invalid value -1 (meaning that the index is not mapped)
#since we do not know how many landmarks we will observe, we allocate a large enough buffer
id_to_state_map = ones(10000, 1)*-1;
state_to_id_map = ones(10000, 1)*-1;

#------------------------------------------ VISUALIZATION ONLY ------------------------------------------
#initialize GUI with initial situation
figure("name", "ekf_slam",    #figure title
       "numbertitle", "off"); #remove figure number
trajectory = [mu(1), mu(2),mu(3)];

measurements = eye();


x=1;
K = camera.camera_matrix;
for id = 1 : 1000%length(l)
  x=0;
  A = [];
  b = [];
    for pose = 0:199
    k=num2str(pose,'%05.f');
    pose=pose+1;
    T1 = inv(camera.cam_transform);
    T2 = inv([cos(t(pose).odometryPose(3)) -sin(t(pose).odometryPose(3)) 0 t(pose).odometryPose(1); sin(t(pose).odometryPose(3)) cos(t(pose).odometryPose(3))  0 t(pose).odometryPose(2); 0 0 1 0; 0 0 0 1]);
    T=T1*T2;%trasformation from world to camera
    filename = strcat("/home/lapo/ProbRobLapo/dataset_for_one_person/meas-", k,  '.dat');
    observations_t=loadmeas(filename);
      for j = 1: length(observations_t)
        x=observations_t(j).observation.x_pose;
        y=observations_t(j).observation.y_pose;
        if observations_t(j).observation.id == id   
          AA= eye(4);
          AA(1:3,1:3) = K*T(1:3,1:3);
          AA(1:3,4) = K*T(1:3,4);
          %A(x:x+2,:) = camera.camera_matrix*WorldToCameraRotation;
          %b(x:x+2,:) = T1(1:3,1:3)* [observations_t(j).observation.x_pose ; observations_t(j).observation.y_pose; 1]+T1(1:3,4)-camera.camera_matrix*WorldToCameraTranslation ;
          %A_i = [x*K(1,1)+K(1,3), x*K(2,1)+K(2,3), R(1,:)+x*R(3,:); y*K(1,2)+K(2,3), y*K(2,2)+K(2,3), R(2,:)+y*R(3,:)];
        
           % Append Ai and bi to A and b
          A_i = [AA(1,1)-x*AA(3,1), AA(1,2)-x*AA(3,2), AA(1,3)-x*AA(3,3);
               AA(2,1)-y*AA(3,1), AA(2,2)-y*AA(3,2), AA(2,3)-y*AA(3,3)];
          A = [A; A_i];
          b_i = [x* AA(3,4)-AA(1,4);y * AA(3,4)-AA(2,4)];
          %b_i = -t(3)*R(3,:)';
          b = [b; b_i];
          x=x+1;
        endif

      endfor
  
    endfor
  size(pinv(A)*b)

  if x>1 && size(pinv(A)*b,1)==3 && size(pinv(A)*b,2)==1
        x_world(id,1:3) = pinv(A)*b;
        h=h+1
  endif %x_world3D(id,:) = x_world;
endfor

x_world
save('variables.mat','x_world')





