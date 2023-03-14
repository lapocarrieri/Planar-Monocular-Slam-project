camera.camera_matrix = [180 0 320;0 180 240; 0 0 1];
camera.inverse_camera_matrix = inv(camera.camera_matrix);
camera.cam_transform =   [0   0   1 0.2; -1   0   0   0;  0  -1   0   0;  0   0   0   1];
camera. z_near = 0;
camera.z_far = 5;
camera.width = 640;
camera.height = 480;
camera
#rotation matrix around x axis
function R=Rx(rot_x)
  c=cos(rot_x);
  s=sin(rot_x);
  R= [1  0  0;
      0  c  -s;
      0  s  c];
 endfunction
 
 #rotation matrix around y axis
 function R=Ry(rot_y)
  c=cos(rot_y);
  s=sin(rot_y);
  R= [c  0  s;
      0  1  0;
      -s  0 c];
 endfunction
 
 #rotation matrix around z axis
 function R=Rz(rot_z)
  c=cos(rot_z);
  s=sin(rot_z);
  R= [ c  -s  0;
       s  c  0;
       0  0  1];
 endfunction
 
 #derivative of rotation matrix around z
 function R=Rx_prime(rot_x)
  dc=-sin(rot_x); #derivative of cos(rot(x)
  ds=cos(rot_x);  #derivative of sin(rot(x)
  R= [0  0  0;
      0  dc  -ds;
      0  ds  dc];
 endfunction
 
 #derivative of rotation matrix around y
 function R=Ry_prime(rot_y)
  dc=-sin(rot_y); #derivative of cos(rot(x)
  ds=cos(rot_y);  #derivative of sin(rot(x)
  R= [dc  0 ds;
      0  0  0;
      -ds  0 dc];
 endfunction
 
 #derivative of rotation matrix around z
 function R=Rz_prime(rot_z)
  dc=-sin(rot_z); #derivative of cos(rot(x)
  ds=cos(rot_z);  #derivative of sin(rot(x)
  R= [ dc  -ds  0;
       ds  dc  0;
       0  0  0];
 endfunction
 
 function R=angles2R(a)
   R=Rx(a(1))*Ry(a(2))*Rz(a(3));
 endfunction;
 
 #from 6d vector to homogeneous matrix
 function T=v2t(v)
     T=eye(4);
     T(1:3,1:3)=angles2R(v(4:6));
     T(1:3,4)=v(1:3);
 endfunction;
 
 #from 7d vector to similiarity matrix Sim(3)
 function S=v2s(v)
   S=eye(4);
   S=v2t(v(1:6));
   S(4,4) = exp(v(7));
 endfunction
 
 function S=skew(v)
   S=[0,    -v(3), v(2);
      v(3),  0,    -v(1);
      -v(2), v(1), 0];
 endfunction
 
 function v=flattenIsometry(T)
 v=zeros(12,1);
 v(1:9)=reshape(T(1:3,1:3)',9,1);
 v(10:12)=T(1:3,4);
 endfunction
 
 function T=unflattenIsometry(v)
   T=eye(4);
   T(1:3,1:3)=reshape(v(1:9),3,3)';
   T(1:3,4)=v(10:12);
 endfunction
 
 function v=flattenIsometryByColumns(T)
 v=zeros(12,1);
 v(1:9)=reshape(T(1:3,1:3),9,1);
 v(10:12)=T(1:3,4);
 endfunction
 
 function T=unflattenIsometryByColumns(v)
   T=eye(4);
   T(1:3,1:3)=reshape(v(1:9),3,3);
   T(1:3,4)=v(10:12);
 endfunction
 
 function M=flatTransformationMatrix(v)
   T=unflattenIsometry(v);
   R=T(1:3,1:3);
   t=T(1:3,4);
   M=eye(12);
   M(1:3,1:3)=R';
   M(4:6,4:6)=R';
   M(7:9,7:9)=R';
   M(10,1:3)=t';
   M(11,4:6)=t';
   M(12,7:9)=t';
 endfunction;
 
 #derivative of rotation matrix w.r.t rotation around x, in 0
 global  Rx0=[0 0 0;
        0 0 -1;
        0 1 0];
 
 #derivative of rotation matrix w.r.t rotation around y, in 0
 global  Ry0=[0 0 1;
        0 0 0;
        -1 0 0];
 
 #derivative of rotation matrix w.r.t rotation around z, in 0
 global  Rz0=[0 -1 0;
        1  0 0;
        0  0 0];
 
 #homogeneous division 
 function p_img = hom(p)
   p_img=p(1:2)/p(3)
 endfunction;
 
 #jacobian of homogeneous division
 function J = J_hom(p)
   x=p(1);
   y=p(2);
   iw=1./p(3);
   iw2=iw^2;
   J = [ iw, 0,  -x*iw2;
         0,  iw, -y*iw2];
 endfunction;
 
 %derivative of icp, euler, left mult, inverse transform and manifold
 function J = J_icp(p)
   J=zeros(3,6);
   J(1:3,1:3)=eye(3);
   J(1:3,4:6)=-skew(p);
 endfunction
 
 #jacobian of the norm function
 function J = J_norm(p)
   n=norm(p);
   J = (1/n) * p';
 endfunction;
 
 #p/norm(p)
 function p2=normalize(p)
   p2=p/norm(p);
 endfunction;
 
 # transform 2 vector
 function v=r2a(R)
   v = [0;0;0];
   s = norm(diag(R));
   if(s < 1e-6)
     v(1) = atan2(R(3,2),R(3,3));
     v(2) = atan2(-R(3,1),s);
     v(3) = atan2(R(2,1),R(1,1));
   else
     v(1) = atan2(-R(2,3),R(2,2));
     v(2) = atan2(-R(3,1), s);
   endif
 endfunction
 
 # TODO check t2v
 function v=t2v(T)
   v = zeros(6,1);
   v(1:3) = T(1:3,4);
   v(4:6) = r2a(T(1:3,1:3));
 endfunction
 
 #from 3D similiarity 3D to 7d vector
 function v=s2v(S)
   v = zeros(7,1);
   v(1:6) = t2v(S); # not working as expected
   v(7) = log(S(4,4));
 endfunction
 
 # jacobian of the normalize function (it's square)
 function J = J_normalize(p)
   n=norm(p);
   in=1./n;
   in3=in^3;
   J = eye(size(p,1))*in - in3*p*p';
 endfunction;
 
 # cartesian2polar:
 # p: coordinates in cartesian space
 # p_polar: coordinates in polar space (range, azimuth, elevation)
 function p_polar = c2p(p)
   n2_xy=p(1)^2+p(2)^2;  # x^2+y^2
   n2 = n2_xy+p(3)^2;    # x^2+y^2+z^2
   n = sqrt(n2);         # norm(p)
   n_xy = sqrt(n2_xy);   # norm(p(1:2)
   p_polar =  [n;
               atan2(p(2),p(1));
               atan2(p(3), n_xy)
              ];
 endfunction;
 
 function J = J_c2p(p)
   n2_xy=p(1)^2+p(2)^2;
   n2 = n2_xy+p(3)^2;
   n = sqrt(n2);
   n_xy = sqrt(n2_xy);
   n32_xy = n_xy^3;
   
   J=zeros(3,3);
   J(1,:)=(1/n)*p';
   J(2,1:2)= 1./n2_xy * [-p(2) p(1)];
   J(3,:) = (n2_xy/n2) * [
                            -p(3)*p(1)/n32_xy;
                            -p(3)*p(2)/n32_xy;
                            1./n_xy];
 endfunction;
 
%(minimal) size of pose and landmarks
global pose_dim=6;
global landmark_dim=3;


# retrieves the index in the perturbation vector, that corresponds to
# a certain pose
# input:
#   pose_index:     the index of the pose for which we want to compute the
#                   index
#   num_poses:      number of pose variables in the state
#   num_landmarks:  number of pose variables in the state
# output:
#   v_idx: the index of the sub-vector corrsponding to 
#          pose_index, in the array of perturbations  (-1 if error)
function v_idx=poseMatrixIndex(pose_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;

  if (pose_index>num_poses)
    v_idx=-1;
    return;
  endif;
  v_idx=1+(pose_index-1)*pose_dim;
endfunction;


# retrieves the index in the perturbation vector, that corresponds to
# a certain landmark
# input:
#   landmark_index:     the index of the landmark for which we want to compute the
#                   index
#   num_poses:      number of pose variables in the state
#   num_landmarks:  number of pose variables in the state
# output:
#   v_idx: the index of the perturnation corrsponding to the
#           landmark_index, in the array of perturbations
function v_idx=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks)
  global pose_dim;
  global landmark_dim;
  if (landmark_index>num_landmarks)
    v_idx=-1;
    return;
  endif;
  v_idx=1 + (num_poses)*pose_dim + (landmark_index-1) * landmark_dim;
endfunction;

function [is_behind, z_hat] = project(Xr, Xl, K)
  p_camera_frame = Xr(1:3,1:3) * Xl + Xr(1:3,4);
  p_projected = K * p_camera_frame;
  is_behind = true;
  z_hat = [-1;-1];
  if(p_camera_frame(3) > 0)
    z_hat = p_projected(1:2)/p_projected(3);
    is_behind = false;
  endif
endfunction




# error and jacobian of a measured landmark
# input:
#   Xr: the robot pose (4x4 homogeneous matrix)
#   Xl: the landmark pose (3x1 vector, 3d pose in world frame)
#   z:  measured position of landmark
# output:
#   e: 2x1 is the difference between prediction and measurement
#   Jr: 2x6 derivative w.r.t a the error and a perturbation on the
#       pose
#   Jl: 2x3 derivative w.r.t a the error and a perturbation on the
#       landmark
function [is_behind,e,Jr,Jl]=errorAndJacobian(Xr, Xl, z, K)
   R=Xr(1:3,1:3);
   t=Xr(1:3,4);
   [is_behind,z_hat]=project(Xr,Xl,K); #prediction
   e = zeros(2,1);
   Jl = zeros(2,3);
   Jr = zeros(2,6);
   if (is_behind)
    return
   endif
   p_camera_frame = R * Xl + t;
   p_projected = K * p_camera_frame;
   e=z_hat-z;
   Jr=zeros(2,6);
   Jl=zeros(2,3);
   # TODO: fill the jacobians
   inverse_z = 1/p_projected(3);
   inverse_square_z = inverse_z * inverse_z;
   J_proj = [inverse_z, 0, -p_projected(1) * inverse_square_z; 
             0, inverse_z, -p_projected(2) * inverse_square_z];

   Jicp = J_icp(p_camera_frame);

   Jr = J_proj * K * Jicp;
   Jl = J_proj * K * R;
endfunction;


# implementation of the boxplus
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the robot poses (4x4xnum_poses: array of homogeneous matrices)
#   XL: the landmark pose (3xnum_landmarks matrix of landmarks)
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   dx: the perturbation vector of appropriate dimensions
#       the poses come first, then the landmarks
# output:
#   XR: the robot poses obtained by applying the perturbation
#   XL: the landmarks obtained by applying the perturbation
function [XR, XL]=boxPlus(XR, XL, num_poses, num_landmarks, dx)
  global pose_dim;
  global landmark_dim;
  for(pose_index=1:num_poses)
    pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
    dxr=dx(pose_matrix_index:pose_matrix_index+pose_dim-1);
    XR(:,:,pose_index)=v2t(dxr)*XR(:,:,pose_index);
  endfor;
  for(landmark_index=1:num_landmarks)
    landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);
    dxl=dx(landmark_matrix_index:landmark_matrix_index+landmark_dim-1,:);
    XL(:,landmark_index)+=dxl;
  endfor;
endfunction;

# implementation of the optimization loop with robust kernel
# applies a perturbation to a set of landmarks and robot poses
# input:
#   XR: the initial robot poses (4x4xnum_poses: array of homogeneous matrices)
#   XL: the initial landmark estimates (3xnum_landmarks matrix of landmarks)
#   Z:  the measurements (3xnum_measurements)
#   associations: 2xnum_measurements. 
#                 associations(:,k)=[p_idx,l_idx]' means the kth measurement
#                 refers to an observation made from pose p_idx, that
#                 observed landmark l_idx
#   num_poses: number of poses in XR (added for consistency)
#   num_landmarks: number of landmarks in XL (added for consistency)
#   num_iterations: the number of iterations of least squares
#   damping:      damping factor (in case system not spd)
#   kernel_threshod: robust kernel threshold

# output:
#   XR: the robot poses after optimization
#   XL: the landmarks after optimization
#   chi_stats: array 1:num_iterations, containing evolution of chi2
#   num_inliers: array 1:num_iterations, containing evolution of inliers
function [XR, XL, chi_stats, num_inliers]=doBundleAdjustment(XR, XL, Z, K,
							associations, 
							num_poses, 
							num_landmarks, 
							num_iterations, 
							damping, 
							kernel_threshold)
  global pose_dim;
  global landmark_dim;

  chi_stats=zeros(1,num_iterations);
  num_inliers=zeros(1,num_iterations);
  # size of the linear system
  system_size=pose_dim*num_poses+landmark_dim*num_landmarks; 
  for (iteration=1:num_iterations)
    H=zeros(system_size, system_size);
    b=zeros(system_size,1);
    chi_stats(iteration)=0;
    for (measurement_num=1:size(Z,2))
      pose_index=associations(1,measurement_num);
      landmark_index=associations(2,measurement_num);
      z=Z(:,measurement_num);
      Xr=XR(:,:,pose_index);
      Xl=XL(:,landmark_index);
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

      pose_matrix_index=poseMatrixIndex(pose_index, num_poses, num_landmarks);
      landmark_matrix_index=landmarkMatrixIndex(landmark_index, num_poses, num_landmarks);

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
    size_H=size(H,1)
    rank_H=rank(H)
    H+=eye(system_size)*damping;
    dx=zeros(system_size,1);
    
    % we solve the linear system, blocking the first pose
    % this corresponds to "remove" from H and b the locks
    % of the 1st pose, while solving the system

    dx(pose_dim+1:end)=-(H(pose_dim+1:end,pose_dim+1:end)\b(pose_dim+1:end,1));
    [XR, XL]=boxPlus(XR,XL,num_poses, num_landmarks, dx);
  endfor

  % visualize the sparsity pattern of the approximate hessian
  visual_H = ones(system_size,system_size);
  for r=1:system_size
    for c=1:system_size
      if(abs(H(r,c)) > 0)
        visual_H(r,c) = 0;
      endif
    endfor
  endfor
  imshow(visual_H)
endfunction

function i = plotState(XL, XL_guess, XL_gt)
#plot landmarks
hold on;
plot3(XL(1,:),XL(2,:),XL(3,:),'b*',"linewidth",2);
hold on;
plot3(XL_guess(1,:),XL_guess(2,:),XL_guess(3,:),'ro',"linewidth",2);
hold on;
plot3(XL_gt(1,:),XL_gt(2,:),XL_gt(3,:),'g*',"linewidth",2);
hold on;
legend("estimate","initial guess","ground truth")
i = 1;
endfunction
function out = landmarkss(varargin)
  %There is always an id
  out.id = varargin{1}; %id
  Nvarargin = length(varargin);
  %if we call the fx with 3 argument, we are populating
  %the landmark absolute position structures.
  if(Nvarargin == 2)
    var2 = varargin{2};
    if(length(var2) == 1)
      out.bearing = var2; %bearing only
    end
    if(length(var2) == 3)
    out.x_pose = var2(1); %x-pose
    out.y_pose = var2(2);
    out.z_pose = var2(3); %y-pose
    end
  end  

end
function out = extractLandmark(elements)
  id = str2double(elements{1});
  x_pose = str2double(elements{2});
  y_pose = str2double(elements{3});
  z_pose = str2double(elements{4});
  out = landmarkss(id,[x_pose,y_pose,z_pose]);
end



function [landmarks] = loadworld(filepath)
    fid = fopen(filepath, 'r');
    i_vert_xy = 0;
    while true
		%get current line
		c_line = fgetl(fid);

		%stop if EOF
		if c_line == -1
			break;
		endif

		%Split the line using space as separator
		elements = strsplit(c_line,' ');
        if exist('landmarks')
        				landmarks(end+1) = extractLandmark(elements);
		else
						landmarks(1) = extractLandmark(elements);
		endif
	    i_vert_xy = i_vert_xy + 1;

    end

end
function out = landmark(varargin)
  %There is always an id
  out.id = varargin{1}; %id
  Nvarargin = length(varargin);
  %if we call the fx with 3 argument, we are populating
  %the landmark absolute position structures.
  if(Nvarargin == 2)
    var2 = varargin{2};
    if(length(var2) == 1)
      out.bearing = var2; %bearing only
    end
    if(length(var2) == 2)
    out.x_pose = var2(1); %x-pose
    out.y_pose = var2(2); %y-pose
    end
  end  

end
function out = observation(from_id, land_id, obs)
  out.pose_id = from_id;
  out.observation = landmark(land_id, obs);
end
function out = extractPoint(elements)
  id = str2double(elements{2});
  land_id = str2double(elements{3});
  x_p = str2double(elements{4});
  y_p = str2double(elements{5});
  out = observation(id,land_id, [x_p; y_p]);
end
function [observations] = loadmeas(filepath)
  fid = fopen(filepath, 'r');
  i_vert_xy = 0;
  curr_id = -1;
  point = 'point';
  observations = [];
  while true
  %get current line
  c_line = fgetl(fid);

  %stop if EOF
  if c_line == -1
    break;
  endif
  
  %Split the line using space as separator
  elements = strsplit(c_line,' ');
  switch(elements{1})
    case point
      current_obs = extractPoint(elements);
      if exist('observations')			
        if current_obs.pose_id == curr_id
          observations(end).observation = [observations(end).observation; current_obs.observation]; 
        else
          observations = [observations; current_obs];
          curr_id = observations(end).pose_id;
        end
      else
        if current_obs.pose_id == curr_id
          observations(1).observation(2) = current_obs.observation;
        else
            observations(2) = current_obs;
          curr_id = observations(1).pose_id;
        end
      endif
    end
  end
  
end
function [is_behind, z_hat] = project(Xr, Xl, K)
  p_camera_frame = Xr(1:3,1:3) * Xl + Xr(1:3,4);
  p_projected = K * p_camera_frame;
  is_behind = true;
  z_hat = [-1;-1];
  if(p_camera_frame(3) > 0)
    z_hat = p_projected(1:2)/p_projected(3);
    is_behind = false;
  endif
endfunction
function out = v2T(v)
	out = [cos(v(3)) -sin(v(3)) 0 v(1); sin(v(3)) cos(v(3))  0 v(2); 0 0 0 0; 0 0 0 1];
	  
endfunction

function out = landmark(varargin)
  %There is always an id
  out.id = varargin{1}; %id
  Nvarargin = length(varargin);
  %if we call the fx with 3 argument, we are populating
  %the landmark absolute positionf structures.
  if(Nvarargin == 2)
    var2 = varargin{2};
    if(length(var2) == 1)
      out.bearing = var2; %bearing only
    end
    if(length(var2) == 2)
    out.x_pose = var2(1); %x-pose
    out.y_pose = var2(2); %y-pose
    end
  end  

end
%qua comincia lo slam leggendo la trajectory
function out = extractTransition(elements)
  id = str2double(elements{2});
  x = str2double(elements{3});
  y = str2double(elements{4});
  z = str2double(elements{5});
  xG = str2double(elements{6});
  yG = str2double(elements{7});
  zG = str2double(elements{8});
  out = transition(id,[x,y,z],[xG,yG,zG]);
end
function out = transition(id, o, g)
  out.id = id;
  out.odometryPose = o;
  out.GroundTruthPose = g;
end
function [transitions] = loadtrajectory(filepath)


	%open the file
	fid = fopen(filepath, 'r');
   
	curr_id = -1;

	

	while true
		%get current line
		c_line = fgetl(fid);

		%stop if EOF
		if c_line == -1
			break;
		end

		%Split the line using space as separator
		elements = strsplit(c_line,' ');

		if exist('transitions')
			        	transitions(end+1) = extractTransition(elements);
		else
						transitions(1) =  extractTransition(elements);
		endif
    end
		
end



function out = observation(from_id, land_id, obs)
  out.pose_id = from_id;
  out.observation = landmark(land_id, obs);
end
function out = extractPoint(elements)
  id = str2double(elements{2});
  land_id = str2double(elements{3});
  x_p = str2double(elements{4});
  y_p = str2double(elements{5});
  out = observation(id,land_id, [x_p; y_p]);
end
curr_id = -1;
function [observations] = loadmeas(filepath)
    fid = fopen(filepath, 'r');
    i_vert_xy = 0;
    curr_id = -1;
    point = 'point';
    observations = [];
    while true
		%get current line
		c_line = fgetl(fid);

		%stop if EOF
		if c_line == -1
			break;
		endif
    
		%Split the line using space as separator
		elements = strsplit(c_line,' ');
    switch(elements{1})
			case point
				current_obs = extractPoint(elements);
				if exist('observations')			
					if current_obs.pose_id == curr_id
						observations(end).observation = [observations(end).observation; current_obs.observation]; 
					else
						observations = [observations; current_obs];
						curr_id = observations(end).pose_id;
					end
				else
					if current_obs.pose_id == curr_id
						observations(1).observation(2) = current_obs.observation;
					else
							observations(2) = current_obs;
						curr_id = observations(1).pose_id;
					end
				endif
      end
    end
    
end

function out = v2T(v)
  out = [cos(v(3)) -sin(v(3)) 0 v(1); sin(v(3)) cos(v(3))  0 v(2); 0 0 1 0; 0 0 0 1];
    
endfunction