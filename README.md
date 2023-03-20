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

Firstly run TriangulationLapo in order to triangulate the position of the landmarks using the Pseudoinverse
```
octave TriangulationLapo.m
```

After the previous command is finished the in variables.mat the landmarks are saved in x_world, now it is possible to run Bundle Adjustment
```
octave BA.m
```
Then to plot the real landmarks against the blue landmarks use
```
octave Plot_figures.m
```



## Authors

Project made by Lapo Carrieri
