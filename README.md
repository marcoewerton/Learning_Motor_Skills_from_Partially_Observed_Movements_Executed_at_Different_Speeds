# Learning_Motor_Skills_from_Partially_Observed_Movements_Executed_at_Different_Speeds
Here is a code that learns to predict trajectories after training with partially observed movements executed at different speeds. This version uses one single phase parameter.

Please, run the script "training.m" first. When 100 EM iterations are over, you can see the approximations that the EM algorithm makes for each of the demonstrations with missing data points. You can press a key on the keyboard to show the next demonstration depicted with the blue asterisks and the approximation depicted in red.
Afterwards, please run the script "test_with_partial_observations.m". Again, by pressing a key on the keyboard, you can see the next test case. In those tests, the algorithm is predicting the rest of the trajectory given the observation of the first 50% of the trajectory.

"recordTrajectory.m" allows the user to draw trajectories with the mouse. You can use this script to generate your own training and test trajectories.

For more information about this work, please refer to http://www.ausy.tu-darmstadt.de/uploads/Team/MarcoEwerton/ewerton_iros_2015_proceedings.pdf
