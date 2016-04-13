#!/usr/bin/env python
import roslib
roslib.load_manifest('rover_traj_follower')

''' member variables of the rrt class '''
import rospy
import time
import sys
import math
import random as rd
from random import *
from numpy import *
import matplotlib.collections as plt_col

import matplotlib.pyplot as plt
from matplotlib import collections, transforms
from docutils.parsers.rst.directives import path
import time
from geometry_msgs.msg import TwistStamped
from geometry_msgs.msg import PoseStamped

class rand_wpt_generator:
    def __init__(self):
        rd.seed(1)
        self.cmdPub = rospy.Publisher('/RR01s/cmd_dest', PoseStamped, latch=True, queue_size=1)
        self.cmdSub = rospy.Subscriber('/RR01s/pose', PoseStamped, self.poseCall)
        self.vehSub = rospy.Subscriber('/RR01s/vel', TwistStamped, self.updateVel)
        self.pose = PoseStamped()
        self.goal = PoseStamped()
        self.vel = TwistStamped()
        self.bnd = array([-7.5, 1.1, -6.9, 1.0])
        self.total_wpt = 10
        self.tol = 0.5
        self.timeLimit = 20.0
        self.start_time = time.time() - 20
        
         # main loop
        r = rospy.Rate(10)
        
        self.count = 0
        while not rospy.is_shutdown() and self.count<self.total_wpt:
            r.sleep()
            
    def generateRndWpt(self):
        valid = False
        while valid == False:
            valid = True
            x = rd.uniform(self.bnd[0]+0.5, self.bnd[1]-0.5)
            y = rd.uniform(self.bnd[2]+0.5, self.bnd[3]-0.5)
            if x < -2.0 and y >-1.8:
                valid = False
            if x < -2.6 and x > -3.9 and y < -3.7 and y > -5.0:
                valid = False
            
        self.goal.pose.position.x = x
        self.goal.pose.position.y = y
        self.goal.pose.position.z = 1.0
        
    def updateVel(self, msg):
        self.vel = msg
        
    def poseCall(self, msg):
        self.pose = msg
        tmp = time.time()
        x_diff = msg.pose.position.x - self.goal.pose.position.x
        y_diff = msg.pose.position.y - self.goal.pose.position.y
        dist = math.sqrt(x_diff **2 + y_diff **2)
        speed = math.sqrt(self.vel.twist.linear.x **2.0 + self.vel.twist.linear.x ** 2.0)
        if (tmp - self.start_time > self.timeLimit and speed < 0.2) or dist < self.tol:
            self.start_time = tmp
            self.publishWpt()
    
    def publishWpt(self):
        self.generateRndWpt()
        self.cmdPub.publish(self.goal)
        print 'published a new goal location'
        print self.goal.pose.position
        print str(self.count) + 'th number of goal location'
        self.count = self.count + 1
        
        
''' for debugging purpose'''
if __name__ == '__main__':
    S = "************ hello world by rnd wpt generator *******************"
    print(S)
    rospy.init_node('rnd_wpt_generator')
    rospy.loginfo('running rand_wpt_generator node')
    rnd_wpt_generator = rand_wpt_generator()
    #rrt_ros_wrapper.rrt_solver.findPath(array([0, 0]), array([-7, -4]))
    #rrt_ros_wrapper.rrt_solver.plot_env()
    #rrt_ros_wrapper.postProcessPath()