#!/usr/bin/env python

''' member variables of the rrt class '''
#import rospy
import time
import random
import sys
from math import *
from random import *
from numpy import *


def saturate (value,min,max):
    if value < min:
        return min
    if value > max:
        return max
    return value

def distance(pt1, pt2):
    return ((pt1[0]-pt2[0]) ** 2.0 + (pt1[1]-pt2[1]) ** 2.0) **0.5

def seg_length(l):
    return (l[0] ** 2.0 + l[1] ** 2.0) **0.5

''' calculate distance between point p3 and 
   line segment p1->p2'''
def distPointToSegment(p1, p2, p3):
    #print p1
    #print p2
    #print p3
    d = p2 - p1
    #print 'd', d
    #print '(p3-p1)', (p3-p1)
    #print 'linalg.norm(d) ** 2', linalg.norm(d) ** 2.0
    if seg_length(d) == 0:
        u = 0.0
    else:
        u = dot(d, (p3-p1)) / (seg_length(d) ** 2.0)
    u = saturate(u,0.0,1.0)
    inter = p1 + u * d
    dist = distance(p3,inter)
    return dist
    


class obs_sphere:
    def __init__ (self, pos, radius):
        # pos is a 1 x 2 numpy array
        self.pos = pos
        self.radius = radius
        
    def ifCollide(self, point1, point2):
        #print self.pos
        r = distPointToSegment(point1, point2, self.pos)
        if r > self.radius:
            return False
        else:
            return True
        
    def isInObs(self, pt1):
        return distance(pt1, self.pos) < self.radius
        

''' whether A->B->C is in counter clockwise pattern '''
def if_ccw(A,B,C):
    # slope AB = u/v ; slope BC = w/x
    # tmp1 = u*x; tmp2 = w*v;
    tmp1 = (C[1] - A[1]) * (B[0] - A[0])
    tmp2 = (B[1] - A[1]) * (C[0] - A[0])
    return tmp1 > tmp2

''' whether the two line segments intersect '''
def if_intersect(A,B,C,D):
    # two line segments intersect iff the infinite extension of any line segment 
    # intersects the other line segment
    if if_ccw(A,C,D) != if_ccw(B,C,D) and \
            if_ccw(A,B,C) != if_ccw(A,B,D):
        return True
    else:
        return False


class obs_poly:
    def __init__(self, numEdges, vertices):
        self.numEdges = numEdges
        # vertices is a n x 2 numpy array
        self.vertices = vertices
    
    
    
    def ifCollide(self, point1, point2):
        # line formed by point1 and point2 cannot intersect any line segment 
        # forming the polygon
        if if_intersect(point1, point2, self.vertices[0,:], self.vertices[self.numEdges-1,:]):
            return True
        for i in range(1,self.numEdges):
            if if_intersect(point1, point2, self.vertices[i-1,:], self.vertices[i]):
                return True
            
        return False
    
    def isInObs(self, pt):
        # assume the polygon is convex, then all edges must be in the same orientation
        # wrt pt
        orientation = if_ccw(self.vertices[self.numEdges-1,:], self.vertices[0,:],pt)
        for i in range(1,self.numEdges):
            if orientation != if_ccw(self.vertices[i-1,:], self.vertices[i],pt):
                return False
        return True
            
    
        
class vertex:
    def __init__(self, cost, pos):
        self.cost = cost
        self.pos = pos
        self.parent = None
        self.angle = None
        
        
        
''' for debugging purpose'''
if __name__ == '__main__':
    S = "************ hello world by rrt_var *******************"
    print(S)
    
    v1 = vertex(0, array([0,0]))
    print 'v1', v1.pos, v1.cost
    
    print "--- spherical obstacle ---"  
    obs_s = obs_sphere(array([-0.3,-0.1]), 0.25)
    print "--- polygon obstacle ---"
    obs_poly = obs_poly(4, array([[2,2],[2,4],[4,4],[4,2]]))
    
    print "--- distance to sphere --- "
    p1 = array([-0.5,-0.25])
    p2 = array([0.2,0.2])
    print distPointToSegment(p1, p2, obs_s.pos)
    
    print "--- whether collide with sphere --- "
    print obs_s.ifCollide(p1,p2)
    
    print "--- whether counter clockwise pattern --- "
    p3 = array([0,0])
    p4 = array([2,0])
    p5 = array([1,1])
    p6 = array([1,-0.5])
    print if_ccw(p3,p4,p5)
    print if_ccw(p3,p4,p6)
    
    print "--- whether two line segments intersect --- "
    print if_intersect(p3,p4,p5,p6)
    
    print "--- whether collide with polygon --- "
    p7 = array([0,3])
    p8 = array([3,3])
    print obs_poly.ifCollide(p7,p8)
    
    print "--- whether point inside circle ---"
    p9 = array([2,3])
    print obs_s.isInObs(p9)
    
    print "--- whehter point inside polygon ---"
    p10 = array([3,4.1])
    print obs_poly.isInObs(p10)