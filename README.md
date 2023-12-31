# Matlab_Sim_robot
This repository is created to calculate data for hexapod robot, like weight, torques, inverse kinematics and forward kinematics of each leg but can be used to simulate work of robotics arm as well.

## Files:

### Kinematyka robota FUNKCJ 
![git_1](https://github.com/Vendicus/Matlab_Sim_robot/assets/119676540/d9ae6a1e-7f25-4a2b-8178-43aca92c24bc)
<p align = "center">
example of inverse kinematics point map generated by alghoritm.
</p>

This is inverse kinematics program used to calculate your robot arm / pedibulator move.
- Program is created by theory of Jacobian matrixes and its usage with inverse kinematics.
- It has weights which can be quiet helpfull if you would like to move more in some directions and less in another, and it used the least squares method to generates next points with more predictibable way but it shouldn't be used when your effector ending points will be beyond workspace of your robotic arm. In that situation your point map could be generated randomly because algorithm will look for possibilites to get to the points but it can't do it in real world.
- you should define required end position and starting angles of your robot arm before start
- you can change the lengths of all three joints

### kinematyka hexa
![git2](https://github.com/Vendicus/Matlab_Sim_robot/assets/119676540/93cac8ac-1137-4d7b-8798-9d0b8f004711)
<p align = "center">
example of forward kinematics solution to arm.
</p>

This is forward kinematics algorithm used to see where in space is your robot arm, you can use it with your pedibulator with 3 joints too.
- to use it, just specify angles and lengths of all three parts.
- algorithm has some testing process of generating jacobian, you can easily see here, how Jacobian is created.

### krzywa Beziera
![git3](https://github.com/Vendicus/Matlab_Sim_robot/assets/119676540/96c58a2c-81e4-4fe1-8940-ed8a97eb23f1)
<p align = "center">
example of Bezier Curve.
</p>

Algorithm created to calculate Bezier Curves.
- you can specify the floating points by changing n but you need to redifine positions of those points and del/create new ones.
- you can change resolution of curve by specifying "s" in first for loop

### robot_fizyka
![git4](https://github.com/Vendicus/Matlab_Sim_robot/assets/119676540/fe75763f-569b-4210-b4f5-7963efff76e2)
<p align = "center">
The front view of hexapod in simulation.
</p>

![up view](https://github.com/Vendicus/Matlab_Sim_robot/assets/119676540/4a9426e1-0863-491f-867c-2e60dd7623e0)
<p align = "center">
The up view of hexapod in simulation. You can see legs which touch the ground (-line), and not touching ground (.-line). The support polygon is draw to see stablisation of robot during move.
</p>

This is advanced program created to calculate masses, weight, reactions forces, joint points, torque in each joint (where servo are located), possible speed of servos and even its accelerations.
- The program has two modes, one with value "symulacja" set to 1 will generate forces, front view and up view of robot.
- The value 0 will generates figures of torques depending on angles, speed and accelerations for each join point.
- you can specify masses of electronics, links, main hull and lengths, and angles too.

![Torques](https://github.com/Vendicus/Matlab_Sim_robot/assets/119676540/9b52f675-f338-4031-bbe8-1b7ff906c0bb)
<p align = "center">
Torques figures for all joints.
</p>

![speed](https://github.com/Vendicus/Matlab_Sim_robot/assets/119676540/2d5cddd4-d595-4fa1-aee1-df6294971dee)
<p align = "center">
Speed figures for B joint with acceleration.
</p>

