#!/usr/bin/env python
import roslib
from reportlab.lib.units import toLength
roslib.load_manifest('rover_traj_follower')
from std_msgs.msg import Float64MultiArray
from std_msgs.msg import MultiArrayDimension

''' member variables of the rrt class '''
import rospy
import time
import sys

import random as rd
from random import *
from numpy import *
from rrt_var import *
from rrt import *
import matplotlib.collections as plt_col

import matplotlib.pyplot as plt
from matplotlib import collections, transforms
from docutils.parsers.rst.directives import path
import time
from geometry_msgs.msg import PoseStamped
import threading
import tf
import math

class rrt_ros:
    def __init__(self, turning_radius, tol, veh_radius):
        self.rrt_solver = rrt(turning_radius, tol, veh_radius)
        
        
        self.pubPath = rospy.Publisher('path', Float64MultiArray)
        self.path_ros = Float64MultiArray()    
        
        self.init_raven()       
        self.rand_sph_obs(20)
        self.cmdSub = rospy.Subscriber('cmd_dest', PoseStamped, self.cmdCall)
        self.poseSub = rospy.Subscriber('pose', PoseStamped, self.poseCall)
        self.pose = PoseStamped()
        self.ifSolved = True
        self.start_time = time.time()
        self.bestCost = None
        self.timeLimit = 20
        self.tol = 0.001
        self.lock = threading.Lock()
        
         # main loop
        self.last_update_time = time.time()
        r = rospy.Rate(1)
        while not rospy.is_shutdown():
            if self.lock.acquire(False) == True:
                self.run()
                self.lock.release()
            
    def run(self):
        #print 'in run function'
        #print 'ifSolved', self.ifSolved
        if self.ifSolved == False:
            self.rrt_solver.updatePath(50)
            tmp = time.time()
            if tmp - self.last_update_time > 1.5:
                self.last_update_time = tmp
                self.bestCostUpdate()
                self.postProcessPath()
                self.rrt_solver.plot_env()

    
    def bestCostUpdate(self):
        if self.rrt_solver.bestNode !=None:
            self.bestCost = self.rrt_solver.bestNode.cost
            if self.bestCost != None:
                cost_reduced = self.bestCost - self.rrt_solver.bestNode.cost
                #if cost_reduced < self.tol:
                #    self.ifSolved = True
                #    print 'problem solved and rrt stopped running because cost change below tol'
                
        timeElasped = time.time() - self.start_time   
        print 'has been solving this problem for ' + str(timeElasped) + 's' 
        if timeElasped > self.timeLimit:
            self.ifSolved = True
            print 'rrt stopped running because run time exceeded ' + str(self.timeLimit) + 's'
                
            
            
    def poseCall(self, msg):
        self.pose = msg  
        
    def cmdCall(self, msg):
        if msg.pose.position.z == 1: # replan
            print '---- got new pose cmd message -----'
            self.lock.acquire()
            x = self.pose.pose.position.x
            y = self.pose.pose.position.y
            x_des = msg.pose.position.x
            y_des = msg.pose.position.y
            
            quat = (self.pose.pose.orientation.x, self.pose.pose.orientation.y, self.pose.pose.orientation.z, self.pose.pose.orientation.w)
            euler = tf.transformations.euler_from_quaternion(quat)
            psi = euler[2]
            
            self.ifSolved = False
            self.bestCost = None
            self.start_time = time.time()
            self.last_update_time = time.time()
            self.rrt_solver.findPath(array([x,y]), array([x_des, y_des]), psi)
            self.lock.release()
            

            
        
        
    
    ''' MIT Raven environment'''
    def init_raven(self):
        self.rrt_solver.bnd = array([-7.5, 1.1, -6.9, 1.0])
        self.rrt_solver.ini_pos = array([0, 0])
        self.rrt_solver.fin_pos = array([-4, -6])
        obs_poly_corner = obs_poly(4, array([[-7.5, 1.0], [-2.0, 1.0], [-2.0, -1.8], [-7.5, -1.8]]))
        obs_poly_pole = obs_poly(4, array([[-2.6, -3.7], [-2.6, -5.0], [-3.9, -5.0], [-3.9, -3.7]]))

        obs_poly_array = []
        obs_poly_array.append(obs_poly_pole)
        obs_poly_array.append(obs_poly_corner)
        obs_sph_array = []
        #obs_sph_array.append(obs_sph1)
        #bs_sph_array.append(obs_sph2)
        #obs_sph1 = obs_sphere(array([0,-3.0]), 0.5)
        #obs_sph2 = obs_sphere(array([-7,-3]), 0.4)
        self.rrt_solver.updateObs(obs_poly_array, obs_sph_array)
        
        # for visualization debugging
        '''
        self.vertices.append(vertex(0,array([0,0])))
        self.vertices.append(vertex(1,array([0,-1])))
        self.vertices[1].parent = self.vertices[0]
        self.vertices.append(vertex(2.15,array([-1,-2])))
        self.vertices[2].parent = self.vertices[1]
        self.vertices.append(vertex(3,array([0.5,-3])))
        self.vertices[3].parent = self.vertices[0]
        self.vertices.append(vertex(4,array([0,-4])))
        self.vertices[4].parent = self.vertices[3]
        self.ini_pos = array([0, 0])
        self.fin_pos = array([-1, -2])
        self.obs_sph_array = []
        self.bestNode = self.vertices[2]
        self.solnFound = True
        '''
        
    
    def rand_sph_obs(self, num_obs):
        rd.seed(1)
        for i in range(num_obs):
            x = rd.uniform(self.rrt_solver.bnd[0]+0.3, self.rrt_solver.bnd[1]-0.3)
            y = rd.uniform(self.rrt_solver.bnd[2]+0.3, self.rrt_solver.bnd[3]-0.3)
            radius = rd.uniform(0.2, 0.5)
            if ((linalg.norm(self.rrt_solver.ini_pos-array([x,y])) < (radius + self.rrt_solver.veh_radius)) \
                or (linalg.norm(self.rrt_solver.fin_pos-array([x,y])) < (radius + self.rrt_solver.veh_radius))):
                i = i-1
                continue              
            
            new_obs_sph = obs_sphere(array([x,y]),radius)
            #print new_obs_sph.pos
            #print new_obs_sph.radius
            self.rrt_solver.obs_sph_array.append(new_obs_sph)
        
        pass
    
    def postProcessPath(self):
        # create a path for the pure-pursuit controller to follow
        # calculate speed at each point in the trajectory
        segment_length = 0.05
        ifFoundPath, best_path = self.rrt_solver.curBestPath() # reversed order
        if ifFoundPath == False:
            return
        nominal_speed = 0.6
        turning_speed = 0.5 * nominal_speed
        accel = 1
        
        speed = []
        for i in range(best_path.shape[0]):
            if i == 0 or i == len(best_path) - 1:
                speed.append(nominal_speed)
            else:
                theta, R = self.findAngle(best_path[i-1,:], best_path[i,:], best_path[i+1,:])
                #print pi-theta
                
                speed.append(turning_speed + (nominal_speed - turning_speed) * saturate((1 - (pi - theta)/(1.57/2.0)), 0.0, 1.0) )
                
        # pad the best path such that inter-point distance is no more than 0.1 m apart
        int_speed = []
        int_speed.append(nominal_speed)
        int_path = best_path[0,:]
        last_point = best_path[0,:]
        last_speed = nominal_speed
        i = 1
        while i<len(speed):
            
            next_point, next_speed = self.nextPoint(last_point, best_path[i,:], \
                                                    last_speed, speed[i], nominal_speed, accel, segment_length)
            if next_point == None:
                last_speed = speed[i]
                last_point = best_path[i,:]
                i = i + 1            
            else:
                last_speed = next_speed
                last_point = next_point
            
            int_speed.append(last_speed)
            int_path = vstack((int_path, last_point))   
        self.convertToRosPath(int_path, int_speed)
        self.publishPath()
        #print 'best path', best_path
        #print 'speed', speed
        #print 'ini_path', int_path
        #rint 'ini_speed', int_speed

    def convertToRosPath(self, processedPath, processedSpeed):
        numPoints = len(processedSpeed)
        # construct path to be published
        d1 = MultiArrayDimension()
        d2 = MultiArrayDimension()
        
        self.path_ros.layout.dim = [d1,d2]
        
        self.path_ros.layout.dim[0].size = numPoints
        self.path_ros.layout.dim[0].stride = numPoints * 4
        self.path_ros.layout.dim[1].size = 4
        self.path_ros.layout.dim[1].stride = 4
        
        self.path_ros.data = 4 * numPoints*[None]
        
        # reverse order
        print 'numPoints', numPoints
        for i in range(numPoints):
            self.path_ros.data[i*4] = processedPath[numPoints - 1 - i, 0]
            self.path_ros.data[i*4+1] = processedPath[numPoints - 1 - i, 1]
            self.path_ros.data[i*4+2] = 0
            self.path_ros.data[i*4+3] = processedSpeed[numPoints - 1 - i]
    
    def nextPoint(self, x1, x2, x1_speed, x2_speed, nominal_speed, accel, segment_len):
        if linalg.norm(x2-x1) < segment_len:
            return None, 0
        else:
            l = x2 - x1;
            slope = l[1] / l[0]
            dx = math.sqrt(segment_len**2.0 / ( 1.0 + slope **2.0))
            if l[0] > 0:
                x3 = array([x1[0] + dx, x1[1] + slope * dx])
            else:
                x3 = array([x1[0] - dx, x1[1] - slope * dx])
                
            return x3, min(nominal_speed, x1_speed + linalg.norm(x3-x2)*accel, x2_speed + linalg.norm(x3-x2)* accel )
    
    def findAngle(self, x1, x2, x3):
        l1 = x1 - x2
        l2 = x3 - x2
        l1_len = linalg.norm(l1)
        l2_len = linalg.norm(l2)
        cos_theta = saturate(dot(l1,l2)/(l1_len*l2_len), -0.9999 , 0.9999)
        theta = math.acos(cos_theta)
        R = 0.5 * min(l1_len,l2_len) / (math.cos(theta/2))
        return theta, R
    
    def publishPath(self):
        self.pubPath.publish(self.path_ros)
        print 'rrt just published a path'

        
        
        
''' for debugging purpose'''
if __name__ == '__main__':
    S = "************ hello world by rrt *******************"
    print(S)
    rospy.init_node('car_rrt')
    rospy.loginfo('running rrt node')
    rrt_ros_wrapper = rrt_ros(0.5, 0.1, 0.5)
    
    
    # debugging
    #rrt_ros_wrapper.rrt_solver.findPath(array([0.0, 0.0]), array([-7, -2]), 1.57)
    #for i in range(20):
    #    rrt_ros_wrapper.rrt_solver.updatePath(50)
    #    rrt_ros_wrapper.rrt_solver.plot_env()
    #rrt_ros_wrapper.postProcessPath()
    #plt.show()