#!/usr/bin/env python

''' member variables of the rrt class '''
#import rospy
import time
import sys

import random as rd
from random import *
from numpy import *
from rrt_var import *
import matplotlib.collections as plt_col

import matplotlib.pyplot as plt
from matplotlib import collections, transforms
from docutils.parsers.rst.directives import path
import time
import math  


class rrt:
    def __init__(self, turning_radius, tol, veh_radius):
        self.bnd = array([-5,5,-5,5])   # [x_min, x_max, y_min, y_max] default
        self.obs_poly_array = list()    # list of polygon obstacles 
        self.obs_sph_array = list()     # list of spherical obstacles
        self.vertices = list()
        self.solnFound = False
        self.bestNode = None
        self.ini_pos = None
        self.fin_pos = None
        self.feasibleProblem = False
        self.turning_radius = turning_radius
        self.tol = tol
        self.veh_radius = veh_radius
        self.init_raven()  
        # for plotting
        self.fig, self.ax = plt.subplots()
        self.draw = False
        self.startObs = False
        
        # for speed testing
        #self.x_random = random.random((1,10000)) * (self.bnd[1] - self.bnd[0]) + self.bnd[0]
        #self.y_random = random.random((1,10000)) * (self.bnd[3] - self.bnd[2]) + self.bnd[2]
        #self.counter = 0;
        # ros
        
    def pruneTree(self):
        if self.solnFound == False:
            return
        bestCost = self.bestNode.cost
        
        i = 0
        while (i<len(self.vertices)):
            if (distance(self.vertices[i].pos, self.ini_pos) + distance(self.vertices[i].pos,self.bestNode.pos)) > bestCost+0.01:
                self.vertices.pop(i)
            else:
                i = i + 1

                

    
    def findPath(self, ini_pos, fin_pos, ini_angle): 
        self.start_time = time.time()
        self.num_update = 0
        print ' ------------- new goal location received -------------' 
        self.solnFound = False
        self.bestNode = None
        self.vertices = []
        self.ini_pos = ini_pos
        self.fin_pos = fin_pos
        self.vertices.append(vertex(0, self.ini_pos))
        self.vertices[0].angle = ini_angle
        # add virtual obstacles confirm to starting orientation
        self.virtual_obs()
        if self.ifFeasible(self.ini_pos) == False or self.ifFeasible(self.fin_pos) == False:
            self.feasibleProblem = False
            print 'infeasible problem!!!'
            return 
        

        
        self.feasibleProblem = True
        self.updatePath(200)
    
    ''' MIT Raven environment'''
    def init_raven(self):
        self.bnd = array([-7.5, 1.1, -6.9, 1.0])
        self.ini_pos = array([0, 0])
        self.fin_pos = array([-4, -6])
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
        self.updateObs(obs_poly_array, obs_sph_array)
        
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

    def virtual_obs(self):
        ini_pos = self.ini_pos
        ini_angle = self.vertices[0].angle
        turn_rad = self.turning_radius
        if self.startObs == True:
            self.obs_poly_array.pop()
            self.obs_sph_array.pop()
            self.obs_sph_array.pop()
            self.startObs = False
        self.startObs = True
        self.obs_sph_array.append(obs_sphere((ini_pos[0]+1.1*turn_rad*math.cos(ini_angle+pi/2.0), \
                                   ini_pos[1]+1.1*turn_rad*math.sin(ini_angle+pi/2.0)), turn_rad))
        
        self.obs_sph_array.append(obs_sphere((ini_pos[0]+1.1*turn_rad*math.cos(ini_angle-pi/2.0), \
                                   ini_pos[1]+1.1*turn_rad*math.sin(ini_angle-pi/2.0)), turn_rad))         
        poly_vertices = array( [[ini_pos[0] + 0.01 * math.cos(ini_angle+pi) - 0.2 * turn_rad * math.sin(ini_angle+pi), 
                                 ini_pos[1] + 0.01 * math.sin(ini_angle+pi) + 0.2 * turn_rad * math.cos(ini_angle+pi)], 
                                [ini_pos[0] + 0.1 * math.cos(ini_angle+pi) - 0.2 * turn_rad * math.sin(ini_angle+pi), 
                                 ini_pos[1] + 0.1 * math.sin(ini_angle+pi) + 0.2 * turn_rad * math.cos(ini_angle+pi)], 
                                [ini_pos[0] + 0.1 * math.cos(ini_angle+pi) + 0.2 * turn_rad * math.sin(ini_angle+pi), 
                                 ini_pos[1] + 0.1 * math.sin(ini_angle+pi) - 0.2 * turn_rad * math.cos(ini_angle+pi)],
                                [ini_pos[0] + 0.01 * math.cos(ini_angle+pi) + 0.2 * turn_rad * math.sin(ini_angle+pi), 
                                 ini_pos[1] + 0.01 * math.sin(ini_angle+pi) - 0.2 * turn_rad * math.cos(ini_angle+pi)],])
        self.obs_poly_array.append(obs_poly(4,poly_vertices))
        
        
    def updatePath(self, max_iter):
        if self.feasibleProblem == False:
            return
        self.num_update = self.num_update + 1;
        start_time = time.time()
        print str(self.num_update) + 'th updatePath() function call'
        eta = 10.0
        gamma_rrg = 4.0 * 2.0 * (1.0 + 1.0/2.0) ** (1.0/2.0)
        # print 'gamma_rrg', gamma_rrg
        sample_goal = 50
        sample_prune = 30
        
        for kk in range(max_iter):
            if len(self.vertices) < 20:
                ifFree, v_rand = self.sampleFree_ini()
            else:
                ifFree, v_rand = self.sampleFree()
            # sample goal once in a while
            goalNode = False
            if (self.solnFound == True and kk % sample_prune == 0):
                self.pruneTree()
            
            if (self.solnFound == False and kk % sample_goal == 0):
                ifFree = True
                goalNode = True
                v_rand = vertex(9999, self.fin_pos)
            if ifFree == False:
                continue

            #v_nearest = self.nearest(v_rand)
            n = len(self.vertices)
            r = min(eta, gamma_rrg * (log(n)/float(n))**(1/2.0))
            v_near, v_nearest = self.near_and_nearest(v_rand, r)
            v_new = self.steer(v_nearest, v_rand)
            # print v_rand.pos
            
            if self.obsFree(v_nearest, v_new) == True:
                
                #print n
                #print gamma_rrg * (log(n)/float(n))**(1/2.0)
                #n = len(self.vertices)
                #r = min(eta, gamma_rrg * (log(n)/float(n))**(1/2.0))
                #v_near = self.near(v_new, r)
                
                self.vertices.append(v_new)
                if goalNode == True:
                    self.solnFound = True
                    self.bestNode = v_new
                    #print 'feasible soln found in ' + str( time.time()-self.start_time) + 's'
                
                v_new.parent = v_nearest
                v_new.cost = v_nearest.cost + self.cost(v_nearest, v_new)
                v_new.angle = math.atan2((v_nearest.pos[1] - v_new.pos[1]) , (v_nearest.pos[0] - v_new.pos[0]))
                
                #print v_new.pos
                for j in range(len(v_near)):
                    if self.obsFree(v_new, v_near[j]) == True:
                        tmp_cost = v_near[j].cost + self.cost(v_near[j],v_new)
                        if tmp_cost < v_new.cost:
                            v_new.cost = tmp_cost
                            v_new.parent = v_near[j]
                            v_new.angle = math.atan2((v_new.pos[1] - v_near[j].pos[1]) , (v_new.pos[0] - v_near[j].pos[0]))
                            
                for k in range(len(v_near)):
                    if self.obsFree(v_new, v_near[k]) == True:
                        tmp_cost = v_new.cost + self.cost(v_new, v_near[k])
                        if tmp_cost < v_near[k].cost:
                            v_near[k].cost = tmp_cost
                            v_near[k].parent = v_new
                            v_near[k].angle = math.atan2((v_near[k].pos[1] - v_new.pos[1]), (v_near[k].pos[0] - v_new.pos[0]))          
        print str(self.num_update) + 'th updatePath() function call returned in ' + str(time.time()-start_time) \
            + 's, tree now has ' + str(len(self.vertices)) + ' vertices, solnFound: ' + str(self.solnFound) #+ 'with cost: ' + str(self.bestNode.cost)
        #if self.solnFound == False:
        #    print 'did not find feasible solution'
        #else:
        #    self.pruneTree()    
        #    print 'feasible soln is found, best path has cost', self.bestNode.cost
            
        #print 'run time: ' + str( time.time()-self.start_time) + 's'
        #print str(len(self.vertices))+ ' vertices in current tree'
        
    
    def obsFree(self, x_nearest, x_new):
        for obs_poly in self.obs_poly_array:
            #print 'x_nearest', x_nearest
            #print 'x_new', x_new
            if obs_poly.ifCollide(x_nearest.pos, x_new.pos) == True:
                return False
            
        for obs_sph in self.obs_sph_array:
            if obs_sph.ifCollide(x_nearest.pos, x_new.pos) == True:
                return False
            
        return True
            
            
    
    def updateObs(self, obs_poly_array, obs_sph_array):
        self.obs_poly_array = list(obs_poly_array) # making a deep copy
        self.obs_sph_array = list(obs_sph_array)   # making a deep copy
        pass
    
    def curBestPath(self):
        if self.solnFound == False:
            return False, array([0,0])
        bestPath = self.fin_pos
        v_tmp = self.bestNode
        while v_tmp.parent != None:
            #print 'best path', bestPath
            #print 'parent pos', v_tmp.parent.pos
            v_tmp = v_tmp.parent
            bestPath = vstack((bestPath, v_tmp.pos))
            #print v_tmp.pos, v_tmp.angle
            
        return True, bestPath
    
    def curTree(self):
        tree = array([0,0,0,0])
        for i in range(1,len(self.vertices)):         
            if self.vertices[i].parent != None:
                tmp = hstack((self.vertices[i].pos, self.vertices[i].parent.pos))
                tree = vstack((tree, tmp))
                
        if tree.ndim == 2:
            tree = tree[1:len(self.vertices),:]
        return tree
            
    def plot_env(self):
        if self.draw == False:
            plt.show(block=False)
            self.draw = True
        plt.cla()
        
        
        ''' bounding box '''
        self.ax.set_xlim(self.bnd[0:2]+array([-0.5, 0.5])) 
        self.ax.set_ylim(self.bnd[2:4]+array([-0.5, 0.5]))    
        self.ax.plot([self.bnd[0], self.bnd[0], self.bnd[1], self.bnd[1], self.bnd[0]], \
                [self.bnd[2], self.bnd[3], self.bnd[3], self.bnd[2], self.bnd[2]], \
                color='k', linewidth=3.0)
        
        ''' polygon obstacles '''
        col_poly_vert = []
        for i in range(len(self.obs_poly_array)):
            col_poly_vert.append(self.obs_poly_array[i].vertices)
            
        col_poly = plt_col.PolyCollection(col_poly_vert, edgecolor='none',facecolor='k')
        self.ax.add_collection(col_poly)
        
        ''' spherical obstacles '''
        for i in range(len(self.obs_sph_array)):
            cir_tmp = plt.Circle((self.obs_sph_array[i].pos), self.obs_sph_array[i].radius, color='k')
            self.ax.add_artist(cir_tmp)

        
        ''' initial and final pos '''
        circle_init = plt.Circle((self.ini_pos), 0.2, color='r')
        circle_end = plt.Circle((self.fin_pos), 0.2, color='g')
        self.ax.add_artist(circle_init)
        self.ax.add_artist(circle_end)
        
        ''' rrt tree '''
        segments = []
        for i in range(len(self.vertices)):
            # root
            if self.vertices[i].parent == None:
                continue
            else:
                segments.append((self.vertices[i].pos, self.vertices[i].parent.pos))
        col_tree = plt_col.LineCollection(segments, linewidths=1, colors='b')
        self.ax.add_collection(col_tree)
        
        ''' best path so far '''
        isSolnFound, bestPath = self.curBestPath()
        self.curTree()
        if isSolnFound == True:
            plt.plot(bestPath[:,0], bestPath[:,1], linewidth=2, color='r', marker='o' )
        plt.draw()
        
        pass
    
    
    def sampleFree_ini(self):
        # loop until found a feasible trajectory
        max_iter = 10
        x_rand = array([0.0,0.0])
        for kk in range(max_iter):
            # sample a new point
            x_rand[0] = saturate(self.ini_pos[0] + self.turning_radius * 2.0 * rd.random() * math.cos(self.vertices[0].angle) \
                                        + self.turning_radius * 4.0 * (rd.random()-0.50) * math.sin(self.vertices[0].angle), \
                                        self.bnd[0], self.bnd[1]) 
            
            x_rand[1] = saturate(self.ini_pos[1] + self.turning_radius * 2.0 * rd.random() * math.sin(self.vertices[0].angle) \
                                        - self.turning_radius * 4.0 * (rd.random()-0.50) * math.cos(self.vertices[0].angle), \
                                        self.bnd[2], self.bnd[3]) 
                                        
                                        
            #x_rand[0] = self.x_random[0,self.counter]
            #x_rand[1] = self.y_random[0,self.counter]
            #self.counter = self.counter + 1
            #print 'x_rand', x_rand
            # test whether point is feasible
            if self.ifFeasible(x_rand) == True:
                return True, vertex(9999, x_rand)
            
            
        return False, vertex(9999, x_rand)
    
    def sampleFree(self):
        # loop until found a feasible trajectory
        max_iter = 1000
        x_rand = array([0.0,0.0])
        for kk in range(max_iter):
            # sample a new point
            x_rand[0] = rd.uniform(self.bnd[0], self.bnd[1])
            x_rand[1] = rd.uniform(self.bnd[2], self.bnd[3])
            #x_rand[0] = self.x_random[0,self.counter]
            #x_rand[1] = self.y_random[0,self.counter]
            #self.counter = self.counter + 1
            #print 'x_rand', x_rand
            # test whether point is feasible
            if self.ifFeasible(x_rand) == True:
                return True, vertex(9999, x_rand)
            
            
        return False, vertex(9999, x_rand)
    
    # whether point is in Free Space
    def ifFeasible(self, pt):
        for i in range(len(self.obs_poly_array)):
            if self.obs_poly_array[i].isInObs(pt) == True:
                return False
            
        for i in range(len(self.obs_sph_array)):
            if self.obs_sph_array[i].isInObs(pt) == True:
                return False
            
        return True

    
    def nearest(self, v_rand):
        #print 'in nearest function'
        dist = distance(self.bnd[0:4:2],self.bnd[1:4:2])
        #print 'distance ', distance
        for i in range(len(self.vertices)):
            tmp = self.cost(self.vertices[i],v_rand)
            if tmp < dist:
                dist = tmp
                v_rand.parent = self.vertices[i]
            
        return v_rand.parent
            
    
    def cost(self, v_from, v_to): 
        return distance(v_from.pos, v_to.pos)    
        dist = distance(v_from.pos, v_to.pos)
        angle = math.atan2((v_to.pos[1] - v_from.pos[1]) , (v_to.pos[0] - v_from.pos[0]))
        rel_angle = abs(angle - v_from.angle)
        if v_to.angle!=None:
            rel_angle = rel_angle + abs(angle - v_to.angle)
        return 0.1 * rel_angle**2 + dist
        #if rel_angle > pi /2.0:
        #    return rel_angle**2 + dist
        #if v_to.angle!=None:
        #    if abs(angle - v_to.angle) > pi/2.0:
        #        return rel_angle**2 + dist
            
        #return dist
        #return ((v_new.pos[0] - v_near.pos[0]) ** 2.0 + \
        #    (v_new.pos[1] - v_near.pos[1]) ** 2.0 ) ** 0.5
        '''
        return distance(v_new.pos, v_near.pos)
        l1 = v_new.pos-v_near.pos
        l1_len = seg_length(l1)
        if v_near.parent == None:
            return l1_len
        l2 =  v_near.parent.pos - v_near.pos
        l2_len = seg_length(l2)
        return l1_len
        # for turning radius
        cos_theta = saturate(dot(l1,l2)/(l1_len*l2_len), -0.9999 , 0.9999)
        theta = math.acos(cos_theta)
        R = 0.5 * min(l1_len,l2_len, self.turning_radius) / (math.cos(theta/2.0))
        if theta < pi/2.0: #or R < self.turning_radius:
            return l1_len + (pi-theta)
        return l1_len
        '''
        
        #if R > self.turning_radius:
        #    return l1_len
        #else:
        #    print theta * 180.0 / pi
        #    return l1_len + (pi-theta) / min(l1_len, l2_len, 1.0)
        #return l1_len
        #return l1_len + turning_penalty * self.turning_radius * (1/R)
        #return l1_len + self.turning_radius * (1.5*(pi-theta))**2
        

        
    
    # make sure path is feasible
    def steer(self, v_nearest, v_rand):
        return v_rand
        
    def near(self, v_new, r):
        #print 'v_new_pos', v_new.pos
        #print 'v_new_cost', v_new.cost
        #print 'r', r
        v_near = []
        for i in range(len(self.vertices)):
            if self.cost(self.vertices[i],v_new) < r and self.vertices[i] != v_new:
                v_near.append(self.vertices[i])
        return v_near
    
    def near_and_nearest(self, v_rand, r):
        v_near = []
        dist = 99999
        for i in range(len(self.vertices)):
            tmp = self.cost(self.vertices[i],v_rand)
            if tmp < r and self.vertices[i] != v_rand:
                v_near.append(self.vertices[i])
            if tmp < dist:
                dist = tmp
                v_rand.parent = self.vertices[i]
        return v_near, v_rand.parent
        
        
''' for debugging purpose'''
if __name__ == '__main__':
    S = "************ hello world by rrt *******************"
    print(S)
    #rospy.init_node('car_rrt')
    #rospy.loginfo('running rrt node')
    rrt_solver = rrt(0.5, 0.1, 0.5)
    #rrt_solver.findPath(array([0, 0]), array([1, 1]), 0)
    #rrt_solver.plot_env()

    # debugging
    rrt_solver.findPath(array([0.0, 0.0]), array([-7, -2]), 1.57)
    for i in range(20):
        rrt_solver.updatePath(50)
        rrt_solver.plot_env()
    plt.show()
    