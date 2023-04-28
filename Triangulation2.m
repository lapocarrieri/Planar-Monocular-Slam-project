
source "./ba_utils.m"
t=loadtrajectory("dataset_for_one_person/trajectoy.dat");


  
  % triangulates a point, passing through two lines
  % one passing through the origin, and having
  % direction vector d1
  % one passing through a point p2, and having
  % direction d2
  function [success, p, e]=triangulatePoint(p2, d1, d2)
    p=zeros(3,1);
    success=false;
                   
    D=[-d1, d2];         #assemble system matrix to find ascissa 
    s=-(D'*D)\(D'*p2);          #s: ascissa of closest point on p1 and p2
    if (s(1)<0 || s(2)<0)
      return;
    endif;
    success=true;
    p1_triangulated=d1*s(1);   # point on 1st line
    p2_triangulated=d2*s(2)+p2; # point on 2nd line
    
    p=0.5*(p1_triangulated+p2_triangulated);          #midpoint
  endfunction;
% se

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
Index=zeros(199,1000);
#------------------------------------------ VISUALIZATION ONLY ------------------------------------------
#initialize GUI with initial situation
figure("name", "ekf_slam",    #figure title
       "numbertitle", "off"); #remove figure number
trajectory = [mu(1), mu(2),mu(3)];

measurements = eye();


x=1;
K = camera.camera_matrix;
% these ten lines are used to create a vector to associate the landmark to the pose so the system is much faster during the iterations
for pose = 1:199
  k=num2str(pose,'%05.f');
  filename = strcat("dataset_for_one_person/meas-", k,  '.dat');
  fclose('all');
  observations_t=loadmeas(filename);
  
  for j = 1: length(observations_t) 
    if(Index(observations_t(j).observation.id+1)==0)
      Index(observations_t(j).observation.id+1)=pose;
    endif
  endfor
endfor

printf("start the triangulation")
for id = 1 : 1000 %length(l)
  id
    firstPoint = false;
    n_success = false;
    success=false;
    if Index(id+1)==0
        continue;
    endif
    for pose = Index(id+1):3:199
      
        if (success)
            break;
        endif
    k=num2str(pose,'%05.f');
    pose=pose+1;


    filename = strcat("dataset_for_one_person/meas-", k,  '.dat');
    fclose('all');
    observations_t=loadmeas(filename);
      for j = 1: length(observations_t)
        if (success)
          break;
        endif
        

        


        if observations_t(j).observation.id == id  
          x=observations_t(j).observation.x_pose;
          y=observations_t(j).observation.y_pose;
          T1 = camera.cam_transform;

          %T2 = [cos(t(pose).GroundTruthPose(3)) -sin(t(pose).GroundTruthPose(3)) 0 t(pose).GroundTruthPose(1); sin(t(pose).GroundTruthPose(3)) cos(t(pose).GroundTruthPose(3))  0 t(pose).GroundTruthPose(2); 0 0 1 0; 0 0 0 1];
          T2 = [cos(t(pose).odometryPose(3)) -sin(t(pose).odometryPose(3)) 0 t(pose).odometryPose(1); sin(t(pose).odometryPose(3)) cos(t(pose).odometryPose(3))  0 t(pose).odometryPose(2); 0 0 1 0; 0 0 0 1];
   
               
            
            if (firstPoint)
                X21 = T1;
                X22 = T2;
                P2_img =  [x;y];

                iK=inv(K);
                #inverse rotation * inverse camera matix
                iX=(((inv(X11)*inv(X12))*X22)*X21);
                iRiK=iX(1:3,1:3)*iK;
              
                #express the points in camera coordinates
                #rotate the direction vector of P2 in other point frame
                p1_cam=iK*[P1_img;1];
                p2_cam=iRiK*[P2_img;1];
                p2=iX(1:3,4);
                n_success= 0;




                [success, p]=triangulatePoint(p2, p1_cam, p2_cam);
                % The triangulation is done using the function given byProf. Grisetti and considering as the center of the graph the camera 1
                % then the point is transformed in the world frame:


                PP=X12*(X11*[p;1])
                if abs(PP(1))>10 || abs(PP(2))>10 || abs(PP(3))>2
                  success=false;
                  break;
                else
                  
                  x_world(id,1:3)=PP(1:3);
                endif
               
                
                
               
                

                
                
                
                            
            
            endif  
            
            firstPoint=true;
            X11 = T1;
            X12 = T2;
            P1_img = [x;y];
            
             
        
        endif

      endfor
  
    endfor
endfor

x_world
save('variables5.mat','x_world')

# The results are saved 

