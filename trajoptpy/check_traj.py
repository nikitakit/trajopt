import trajoptpy.math_utils as mu
import openravepy as rave
import numpy as np

def traj_collisions(traj, robot, n=None):
    """
    Returns the set of collisions. 
    manip = Manipulator or list of indices
    """
    if n is None:
        n = len(traj) * 5
    traj_up = mu.interp2d(np.linspace(0,1,n), np.linspace(0,1,len(traj)), traj)
    _ss = rave.RobotStateSaver(robot)
    
    env = robot.GetEnv()
    col_times = []
    for (i,row) in enumerate(traj_up):
        robot.SetActiveDOFValues(row)
        col_now = env.CheckCollision(robot)
        if col_now: 
            col_times.append(i)
    return col_times
    

def traj_is_safe(traj, robot, n=None):
    return traj_collisions(traj, robot, n) == []
