#script that defines the objects type to be tiled
#for now, includes circles and convex polygons
import numpy as np
import copy


#helper functions

''' return distance between 2 points '''
def distance(pt1, pt2):
    return ((pt1[0]-pt2[0]) ** 2.0 + (pt1[1]-pt2[1]) ** 2.0) **0.5

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

def seg_length(l):
    return (l[0] ** 2.0 + l[1] ** 2.0) **0.5

def saturate (value,min,max):
    if value < min:
        return min
    if value > max:
        return max
    return value

def determinant(p1,p2):
    return p1[0]*p2[1] - p1[1] *p2[0]


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
        u = np.dot(d, (p3-p1)) / (seg_length(d) ** 2.0)
    u = saturate(u,0.0,1.0)
    inter = p1 + u * d
    dist = distance(p3,inter)
    return dist

#circle objects
class obj_circle:
    ''' initialize position and radius of the object'''
    def __init__ (self, pos, radius):
        # pos is a 1 x 2 numpy array
        self.pos = pos
        self.radius = radius
        self.area = np.pi * radius **2

    ''' compute distance to a circle'''
    def dist_to_circle(self,circle):
        assert(isinstance(circle,obj_circle))
        return distance(self.pos,circle.pos) - (self.radius + circle.radius) 

    '''return distance to either a circle or a polygon'''
    def calc_dist(self,item):
        if isinstance(item,obj_circle):
            return self.dist_to_circle(item)
        return self.dist_to_polygon(item)

    ''' check if collide with a circle'''
    def ifCollide_circle(self,circle):
        assert(isinstance(circle,obj_circle))
        return self.dist_to_circle(circle) < 0

    ''' compute distance to a polygon '''
    def dist_to_polygon(self,polygon):
        assert(isinstance(polygon,obj_polygon))

        min_dist = 1e10+0.0 #initialize a big number
        #print 'verts are', polygon.verts
        for i in range(len(polygon.verts)-1):
            p1 = np.array(polygon.verts[i]) + polygon.offset
            p2 = np.array(polygon.verts[i+1]) + polygon.offset

            #print self.pos, 'p1 =', p1, 'p2=', p2
            r = distPointToSegment(p1,p2,self.pos) 
            #print 'points = ', self.pos, p1,p2
            #print 'dist =', r
            if r < min_dist:
                min_dist = r

        return min_dist - self.radius

    ''' compute if collide with a polygon'''
    def ifCollide_polygon(self,polygon):
        assert(isinstance(polygon,obj_polygon))
        #print 'checking if collide polygon, distance =', self.dist_to_polygon(polygon)
        intersect_edge = self.dist_to_polygon(polygon)< 0 
        center_in_poly = polygon.isIn_poly(self.pos)
        #print 'checking if center of circle is in polygon,', center_in_poly
        #print self.pos
        #print polygon.verts, polygon.offset
        return (intersect_edge or center_in_poly)

    ''' return if the object collides with another object'''
    def ifCollide(self,item):
        if isinstance(item,obj_circle):
            return self.ifCollide_circle(item)
        return self.ifCollide_polygon(item)

#polygon objects
class obj_polygon:
    ''' initialize offset and vertices of the object'''
    def __init__ (self,  verts,offset):
        # pos is a 1 x 2 numpy array
        self.verts = verts
        self.offset = offset
        self.area = self.polygon_area(self.verts)

    ''' find the polygon area based on ccw input vertices'''
    def polygon_area(self,verts):
        area =0.0
        for i in range(len(verts)-1):
            area += determinant(np.array(verts[i]),np.array(verts[i+1]))
        return area /2.0

    ''' compute distance to a circle'''
    def dist_to_circle(self,circle):
        assert(isinstance(circle,obj_circle))
        return circle.dist_to_polygon(self)

    ''' check if collide with a circle'''
    def ifCollide_circle(self,circle):
        assert(isinstance(circle,obj_circle))
        return self.dist_to_circle(circle) < 0

    '''return distance to either a circle or a polygon'''
    def calc_dist(self,item):
        if isinstance(item,obj_circle):
            return self.dist_to_circle(item)
        return self.dist_to_polygon(item)

    ''' return if the object collides with another object'''
    def ifCollide(self,item):
        if isinstance(item,obj_circle):
            return self.ifCollide_circle(item)
        return self.ifCollide_polygon(item)

    ''' compute distance to a polygon '''
    def isIn_poly(self, pt):
        # assume the polygon is convex, then all edges must be in the same orientation
        # wrt pt
        #shift pt by offset of the polygon
        pt = pt + self.offset

        #print 'transformed pt is' ,pt
        numEdges = len(self.verts)-1
        #print 'last 2 polygon points are', np.array(self.verts[numEdges-1]),np.array(self.verts[0])
        
        orientation = if_ccw(np.array(self.verts[numEdges-1]),np.array(self.verts[0]),pt)

        #print 'orientation is,', orientation

        for i in range(1,numEdges):
            if orientation != if_ccw(np.array(self.verts[i-1]), np.array(self.verts[i]),pt):
                return False
        return True

    ''' compute if collide with a polygon'''
    def ifCollide_polygon(self,polygon):
        assert(isinstance(polygon,obj_polygon))

        #check if vertices of polygon is in the current object
        '''  
        for i in range(len(polygon.verts)):
            if self.isIn_poly(np.array(polygon.verts[i])):
                return True

        #check if object vertices are contained in the polygon
        for i in range(len(self.verts)):
            if polygon.isIn_poly(np.array(self.verts[i])):
                return True
        '''
        for i in range(len(polygon.verts)-1):
            for j in range(len(self.verts)-1):
                A1 =np.array(polygon.verts[i])+polygon.offset
                A2 =np.array(polygon.verts[i+1])+polygon.offset
                B1 =np.array(self.verts[j])+self.offset
                B2 =np.array(self.verts[j+1])+self.offset
                if if_intersect(A1,A2,B1,B2):
                    return True

        return False


    def dist_to_polygon(self,polygon):
        assert(isinstance(polygon,obj_polygon))

        #check if is collided
        if self.ifCollide_polygon(polygon):
            return 0.0

        min_dist= 1e10+0.0 

        #distance of self verts to polygon edges
        for i in range(len(self.verts)):
            for j in range(len(polygon.verts)-1):
                p1 = np.array(polygon.verts[j]) + polygon.offset
                p2 = np.array(polygon.verts[j+1]) + polygon.offset
                r = distPointToSegment(p1,p2,np.array(self.verts[i])) 
                if r < min_dist:
                    min_dist = r

        #distance of polygon verts to self edges
        for i in range(len(polygon.verts)):
            for j in range(len(self.verts)-1):
                p1 = np.array(self.verts[j]) + self.offset
                p2 = np.array(self.verts[j+1]) + self.offset
                r = distPointToSegment(p1,p2,np.array(polygon.verts[i])) 
                if r < min_dist:
                    min_dist = r

        return min_dist

#list of objects with collision relationship and distances
class object_lists:
    def __init__ (self,circle_list,polygon_list):
        # save circles and polygons
        self.circles = circle_list
        self.polygons = polygon_list
        self.num_circles = len(circle_list)
        self.num_polygons = len(polygon_list)

        #compute distance matrix of every object to other objects
        self.dist_mat = self.generate_distance_matrix()

        #compute collision matrix (0 is not collide, 1 is collide with all other objects)
        #compute collision_free vector (True or False) , whether each object is free
        self.collision_mat,self.collision_vec= self.generate_collision_matrix() 

        #for plotting circles
        self.circ_diameter, self.circ_position, self.circ_collision = self.generate_circ_info()

        #for plotting polygons
        self.poly_verts, self.poly_collision = self.generate_poly_info()

        print self.collision_mat



    def generate_distance_matrix(self):
        all_objects = self.circles + self.polygons
        num_obj= self.num_circles + self.num_polygons
        dist_mat = np.zeros((num_obj, num_obj),dtype = 'float')
        for i in range(num_obj):
            for j in range(num_obj):
                dist_mat[i,j] = all_objects[i].calc_dist(all_objects[j])

        return dist_mat

    def generate_collision_matrix(self):
        all_objects = self.circles + self.polygons
        num_obj= self.num_circles + self.num_polygons
        collision_mat = np.zeros((num_obj, num_obj),dtype = 'float')
        collision_vec = [True] * num_obj
        for i in range(num_obj):
            for j in range(num_obj):
                if all_objects[i].ifCollide(all_objects[j]):
                    collision_mat[i,j] = 1

        temp = np.sum(collision_mat,axis=0)
        for i in range(num_obj):
            if temp[i] == 1:
                collision_vec[i] = False
              
        return collision_mat,collision_vec
        
    def generate_circ_info(self):
        circ_diameter = np.zeros(len(self.circles),dtype = 'float')
        circ_position = np.zeros((len(self.circles),2),dtype = 'float')
        circ_collision = [None] * len(self.circles)

        for i in range(self.num_circles):
            circ_diameter[i] = self.circles[i].radius * 2
            circ_position[i,:] = np.array(self.circles[i].pos)
            circ_collision[i] = self.collision_vec[i]

        return circ_diameter, circ_position, circ_collision

    def generate_poly_info(self):
        poly_verts = [None] * self.num_polygons
        poly_collision = [None] * len(self.polygons)

        for i in range(self.num_polygons):
            temp_verts= copy.deepcopy(self.polygons[i].verts)
            #print 'temp_verts =', temp_verts
            poly_verts[i] = temp_verts
            for j in range(len(self.polygons[i].verts)):
                poly_verts[i][j] = tuple(np.array(temp_verts[j]) + self.polygons[i].offset)
            poly_collision[i] = self.collision_vec[i+self.num_circles]

        return poly_verts, poly_collision

    def print_items_info(self):
        print 'collision matrix =', self.collision_mat
        print 'circles collision =' ,self.circ_collision
        print 'poly collision =' ,self.poly_collision
        all_objects = self.circles + self.polygons

        for i in range(self.num_polygons+ self.num_circles):
            if isinstance(all_objects[i],obj_circle):
                print 'circle object ', i, ',area =', all_objects[i].area, 'pos = ',all_objects[i].pos
            else:
                print 'polygon object ', i, ',area =', all_objects[i].area, 'pos = ', all_objects[i].offset

    '''update delta position'''
    #dx is a numpy array of size num_objects x 2
    def update_delta_pos(self,dx):
        #update circle offsets and polygon offsets
        all_objects = self.circles + self.polygons
        for i in range(self.num_polygons + self.num_circles):
            if isinstance(all_objects[i],obj_circle):
                all_objects[i].pos += dx[i]
            else:
                all_objects[i].offset += dx[i]

        #re-compute all information due to change of positions....................
        #compute distance matrix of every object to other objects
        self.dist_mat = self.generate_distance_matrix()

        #compute collision matrix (0 is not collide, 1 is collide with all other objects)
        #compute collision_free vector (True or False) , whether each object is free
        self.collision_mat,self.collision_vec= self.generate_collision_matrix() 

        #for plotting circles
        self.circ_diameter, self.circ_position, self.circ_collision = self.generate_circ_info()

        #for plotting polygons
        self.poly_verts, self.poly_collision = self.generate_poly_info()




if __name__ == '__main__':
    print 'ckecking objects'

    #circles
    circ1 = obj_circle(np.array([0,0]),1.0)
    circ2 = obj_circle(np.array([1.5,0]),1.0)

    print circ1.ifCollide_circle(circ2)
    print circ1.dist_to_circle(circ2)

    #polygons
    poly1 = obj_polygon([(0,0),(3,0),(3,2),(0,2),(0,0)],np.array([2,0]))
    poly2 = obj_polygon([(0,0),(3,6),(0,4),(0,0)],np.array([3,0]))

    print circ1.ifCollide_polygon(poly1) 
    print circ1.dist_to_polygon(poly1)


    print poly1.ifCollide_polygon(poly2)

    print 'testing collision checking...'

    #initialize a few circles and polygon objects
    circ1 = obj_circle(np.array([0.0,0.0]),1.0)
    circ2 = obj_circle(np.array([1.5,0]),1.0)
    circ3 = obj_circle(np.array([8.0,8.0]),2.0)

    #important: polygon must be defined in ccw order!!!
    poly1 = obj_polygon([(0,0),(3,0),(3,2),(0,2),(0,0)],np.array([2.0,0.0]))
    poly2 = obj_polygon([(0,0),(3,6),(0,4),(0,0)],np.array([3.0,0.0]))
    poly3 = obj_polygon([(0,0),(2,0),(0,1),(0,0)],np.array([-6.0,0.0]))

    print 'should be False, ', circ1.ifCollide(poly3)
    #assert(0)
   
    #form arrays of objects
    circ_list= [circ1,circ2,circ3]
    poly_list= [poly1,poly2,poly3]

    #formulate a list and do all pre-processing (finding distance + collision, etc)
    item_lists = object_lists(circ_list,poly_list)
    item_lists.print_items_info()
    #print item_lists.poly_verts

    delta_x = np.array([[1,0],[1,0],[1,0],[1,0],[1,0],[1,0]],dtype='float')

    item_lists.update_delta_pos(delta_x)
    item_lists.print_items_info()

    #print item_lists.poly_verts

    #debug function
    #print if_intersect([0,0],[2,0],[1,1],[1,-1])
    #print if_intersect([0,0],[2,0],[3,1],[3,-1])
    print 'debugging function...'
    
    
    


