# Planar monocular SLAM project for Sapienza
Project for the Probabilistic project course at Sapienza


## Description

Given a trajectory with 200 poses and 1000 landmarks to be observed the task is to triangulate the position of the landmarks and make Bundle adjustment using the truiangulate positions and the positions of the robot.


## Getting Started

### Dependencies
On PC terminal: 
```
git clone https://github.com/lapocarrieri/Probabilistic-Robotics-Sapienza
cd Probabilistic-Robotics-Sapienza
unzip dataset_for_one_person.zip   
```

### Installing

Install the dataset unzipping 09-PlanarMonocularSlam_one_person/dataset_for_one_person inside the actual folder

### Executing program

Firstly run Triangulation2 in order to triangulate the positions  using the triangulation between two points using the function present in Tools. in order to use multiplu points in Triangulation.m  the Pseudoinverse is used ( but it takes hours and it is not precise).
If you want to see the Bundle Adjustment it is possible to skip this part since the landmarks are already stored in variables2.mat in the variable x_world.
 Uncommenting the line 105 it is possible to see that the triangulation is perfectly matched if the ground truthposition is considered.
```
octave Triangulation2.m
```

After the previous command is finished it is possible to run Bundle Adjustment ( it is possible to change num_iterations in order to increase the performances of the system)
```
octave BA.m
```

Then to plot the real landmarks against the blue landmarks use
```
octave plot_figures.m
```
![Landmarks](https://user-images.githubusercontent.com/56505429/235121793-fb2d52df-384f-435b-94b9-37afcb73fbb9.png)


Finally to plot the real landmarks against the blue landmarks use
```
octave PlotTrajectory.m
```
![Trajectory](https://user-images.githubusercontent.com/56505429/235121749-c354a0cf-f4bf-447a-8f70-62ad11f7100a.png)


### Results

The results are not particularly interesting since the error is still evident after some iterations. Anyway the problem is the triangulated points that are not very precise due to the error in the odometry. Making the BA using real landamrks gives error 0 so improving the initialization it will be possible to improve the system.
For the Robot position the results are good since the RMSE is small while for the landmarks it can be improved

![image](https://user-images.githubusercontent.com/56505429/235121047-fa3f8e6c-34a7-4753-8359-b106d2709a39.png)


### Future improvements

A way to simply improve the system is triangulate the points using more than two points, this require much time computing but, since in the project the triangulation is done between two subsequent point is not so precise in respect to 2 far positions.




## Authors

Project made by Lapo Carrieri
