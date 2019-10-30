# Probabilistic_Planners_for_n_DOF_arm
PRM

**Compilation Instructions**

>>mex planner.cpp

**Run Instructions**

>>startQ = [pi/2 pi/4 pi/2 pi/4 pi/2];   
>>goalQ = [pi/8 3*pi/4 pi 0.9*pi 1.5*pi]; 
>>planner_id = 0;     ( 0:RRT 1:RRT-Connect 2:RRT* 3:PRM)
>>runtest('map1.txt',startQ, goalQ, planner_id);
