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

''' find intersection of two line segemets given by 4 points'''
def find_intersect(A,B,C,D):
    den = ((A[0]-B[0])*(C[1]-D[1])-(A[1]-B[1])*(C[0]-D[0]))

    px = (A[0]*B[1] - B[0]*A[1]) * (C[0]-D[0]) - (A[0]-B[0]) * (C[0]*D[1]-D[0]*C[1])
    px = px/ den

    py = (A[0]*B[1] - B[0]*A[1]) * (C[1]-D[1]) - (A[1]-B[1]) * (C[0]*D[1]-D[0]*C[1])
    py = py/ den

    return np.array([px,py])
    

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

    ''' set circle centre '''
    def set_pos(self, new_pos):
        self.pos = new_pos

    '''increment the position of circle position'''
    def increment_pos(self,disp):
        self.pos += disp

    '''print information of the current object'''
    def print_info(self):
        print 'circle: ,', 'center = ', self.pos, ',radius =', self.radius

    ''' return the height of the bounding box'''
    def get_height(self):
        return 2*self.radius

    '''return the width of the bounding box'''
    def get_width(self):
        return 2*self.radius

    ''' return the position of the circle radius'''
    def get_position(self):
        return self.pos

    '''untangle collision with an object '''
    def untangle_collision(self,obj,opt='self'):
        if isinstance(obj,obj_circle):
            untangle_collision_circle(self,obj,opt)
        elif isinstance(ob,obj_polygon):
            untangle_collision_polygon(self,obj,opt)
        else:
            assert(0)

    '''to resolve collision with a circle '''
    def untangle_collision_circle(self,circle,opt = 'self'):
        #option parameter = self, other, or equal (which to shift)
        assert(isinstance(circle,obj_circle))
        assert(self.ifCollide_circle(circle))
        
        #resolve collision in direction of centroid
        dir = (circle.pos - self.pos)/np.linalg.norm(circle.pos - self.pos)
        dist = self.dist_to_circle(circle)

        if opt == 'self':
            self.pos = self.pos + np.abs(dist) *(-dir)
        elif opt == 'other':
            circle.pos = circle.pos + np.abs(dist) * (dir)
        else:
            self.pos = self.pos + np.abs(dist) *(-dir)/2
            circle.pos = circle.pos + np.abs(dist) * (dir)/2

    ''' to resolve collision with a polygon'''
    def untangle_collision_polygon(self,polygon,opt = 'self'):
        #option parameter = self, other, or equal (which to shift)
        assert(isinstance(polygon,obj_polygon))
        assert(self.ifCollide_polygon(polygon))
        
        #resolve collision in direction of centroid
        dir = (polygon.centroid + polygon.offset - self.pos) / \
              np.linalg.norm(polygon.centroid + polygon.offset- self.pos)

        #update center dependingon whether circle center is inside the polygon
        if polygon.isIn_poly(self.pos):
            center = self.pos + self.radius * 10 * (-dir)
        else:
            center = self.pos

        #find intersection point on polygon edge and the centroid, center line
        for i in range(len(polygon.verts)-1):
            v1 = np.array(polygon.verts[i]) + polygon.offset
            v2 = np.array(polygon.verts[i+1]) +polygon.offset
            if if_intersect(center,polygon.centroid+polygon.offset,v1,v2):
                p_star = find_intersect(center,polygon.centroid+polygon.offset,v1,v2)
                break
        
        new_pos = p_star + self.radius *(-dir)
        dist = np.linalg.norm(new_pos-self.pos)

        if opt == 'self':
            self.pos = self.pos + np.abs(dist) *(-dir)
        elif opt == 'other':
            polygon.offset += np.abs(dist) * (dir)
        else:
            self.pos = self.pos + np.abs(dist) *(-dir)/2
            polygon.offset += np.abs(dist) * (dir)/2
        
                                  

    ''' check is the polygon intersects or is out of the template'''
    def isIn_template(self,template):
        #template is given in a numpy array in form of [xmin,xmax,ymin,ymax]
        xmin = self.pos[0]-self.radius;  xmax = self.pos[0]+self.radius; 
        ymin = self.pos[1]-self.radius;  ymax = self.pos[1]+self.radius; 
        eps=1e-8
        return (xmin > template[0]-eps and xmax < template[1]+eps and ymin > template[2]-eps and ymax < template[3]+eps)

    ''' compute distance to a circle'''
    def dist_to_circle(self,circle):
        assert(isinstance(circle,obj_circle))
        return distance(self.pos,circle.pos) - (self.radius + circle.radius) 

    ''' compute distance to center of an object'''
    def centroid_dist(self,obj):
        if isinstance(obj,obj_circle):
            return np.linalg.norm(self.pos - obj.pos)
        else:
            return np.linalg.norm(self.pos - (obj.centroid + obj.offset))

    '''compute displacement direction to an object'''
    def centroid_dir(self,obj):
        if isinstance(obj,obj_circle):
            return (self.pos -obj.pos)/np.linalg.norm(self.pos - obj.pos)
        else:
            return (obj.centroid+obj.offset - self.pos)/np.linalg.norm(self.pos - (obj.centroid + obj.offset))

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
        self.centroid = self.polygon_centroid(self.verts,self.area)
        self.bounding_box = self.find_bounding_box(self.verts)

    '''set polygon offset'''
    def set_pos(self,new_offset):
        self.offset = new_offset
    
    '''increment the position of polygon offset'''
    def increment_pos(self,disp):
        self.offset += disp

    ''' return the height of the bounding box'''
    def get_height(self):
        return self.bounding_box[3] - self.bounding_box[2]

    '''return the width of the bounding box'''
    def get_width(self):
        return self.bounding_box[1] - self.bounding_box[0]

    '''return the current position of the polygon'''
    def get_position(self):
        return self.centroid + self.offset

    ''' compute distance to center of an object'''
    def centroid_dist(self,obj):
        if isinstance(obj,obj_circle):
            return np.linalg.norm(self.centroid + self.offset - obj.pos)
        else:
            return np.linalg.norm(self.centroid + self.offset - (obj.centroid + obj.offset))

    '''compute displacement direction to an object'''
    def centroid_dir(self,obj):
        if isinstance(obj,obj_circle):
            return (self.centroid + self.offset - obj.pos)/ \
                   np.linalg.norm(self.centroid + self.offset - obj.pos)
        else:
            return (self.centroid + self.offset - (obj.centroid + obj.offset))/ \
                   np.linalg.norm(self.centroid + self.offset - (obj.centroid + obj.offset))

    '''untangle collision with an object '''
    def untangle_collision(self,obj,opt='self'):
        if isinstance(obj,obj_circle):
            untangle_collision_circle(self,obj,opt)
        elif isinstance(ob,obj_polygon):
            untangle_collision_polygon(self,obj,opt)
        else:
            assert(0)
    

    '''untangle collision with a cirlce'''
    def untangle_collision_circle(self,circle,opt ='self'):
        assert(isinstance(circle,obj_circle))
        if opt == 'self':
            circle.untangle_collision_polygon(self,'other')
        elif opt == 'other':
            circle.untangle_collision_polygon(self,'self')
        else:
            circle.untangle_collison_polygon(self,'equal')


    '''untangle collision with another polygon'''
    def untangle_collision_polygon(self,polygon,opt = 'self'):
        assert(isinstance(polygon,obj_polygon))
        assert(self.ifCollide_polygon(polygon))
        
        #resolve collision in direction of centroid
        dir = (polygon.centroid + polygon.offset - (self.centroid + self.offset))
        dir = dir / np.linalg.norm(dir)

        print 'polygon direction is:', dir

        #update center dependingon whether circle center is inside the polygon
        if polygon.isIn_poly(self.offset+self.centroid):
            center1 = self.offset + self.centroid + 10 * (-dir)
        else:
            center1 = self.offset + self.centroid

        #update center dependingon whether circle center is inside the polygon
        if self.isIn_poly(polygon.offset):
            center2 = polygon.offset + polygon.centroid + 10 * (dir)
        else:
            center2 = polygon.offset + polygon.centroid 

        print 'centers are', center1,center2
        #find intersection point on polygon edge and the centroid, center line
        for i in range(len(polygon.verts)-1):
            v11 = np.array(polygon.verts[i]) + polygon.offset
            v12 = np.array(polygon.verts[i+1]) +polygon.offset
            if if_intersect(center1,center2,v11,v12):
                break

        for i in range(len(self.verts)-1):
            v21 = np.array(self.verts[i]) + self.offset
            v22 = np.array(self.verts[i+1]) +self.offset
            if if_intersect(center1,center2,v21,v22):
                break

        print 'valid intersecting indices are:', v11,v12,v21,v22
        
        dist11=0.0; dist12=0.0; dist21=0.0; dist22=0.0
        if self.isIn_poly(v11):
            print 'v11 inside'
            dist11 = distPointToSegment(v21,v22,v11)
        if self.isIn_poly(v12):
            print 'v12 inside'
            dist12 = distPointToSegment(v21,v22,v12)
        if polygon.isIn_poly(v21):
            print 'v21 inside'
            dist21 = distPointToSegment(v11,v12,v21)
        if polygon.isIn_poly(v22):
            print 'v22 inside'
            dist22 = distPointToSegment(v11,v12,v22)
        
        dist = np.max(np.array([dist11,dist12,dist21,dist22]))
        print 'distance to shift =', dist
        if opt == 'self':
            self.offset +=  np.abs(dist) *(-dir)
        elif opt == 'other':
            polygon.offset += np.abs(dist) * (dir)
        else:
            self.offset += np.abs(dist) *(-dir)/2
            polygon.offset += np.abs(dist) * (dir)/2

        

    '''print information of the current object'''
    def print_info(self):
        print 'polygon: ,', 'verts = ', self.verts,\
              ',shifted centroid =', self.centroid+ self.offset, 'offset =', self.offset

    ''' check is the polygon intersects or is out of the template'''
    def isIn_template(self,template):
        #template is given in a numpy array in form of [xmin,xmax,ymin,ymax]
        xmin = self.bounding_box[0]+self.offset[0];  xmax = self.bounding_box[1]+self.offset[0]; 
        ymin = self.bounding_box[2]+self.offset[1];  ymax = self.bounding_box[3]+self.offset[1]; 

        #print 'template is', template
        #print 'bounding box is', xmin,ymin,xmax,ymax
        eps = 1e-8
        return (xmin > template[0]-eps and xmax < template[1]+eps and ymin > template[2]-eps and ymax < template[3]+eps)

    ''' find the bounding box [xmin,xmax,ymin,ymax] of the polygon with first vertex at (0,0)'''
    def find_bounding_box(self,verts):
        temp = np.array(verts)
        min_vec= np.min(temp,axis=0)
        max_vec= np.max(temp,axis=0)
        #print temp,min_vec,max_vec
        #print np.array([min_vec[0],max_vec[0],min_vec[1],max_vec[1]])
        return np.array([min_vec[0],max_vec[0],min_vec[1],max_vec[1]])

    ''' find the polygon centroid on ccw input vertices'''
    def polygon_centroid(self,verts,area):
        cx=0; cy=0
        for i in range(len(verts)-1):
            x_i = self.verts[i][0]; y_i = self.verts[i][1]
            x_ip1 = self.verts[i+1][0]; y_ip1 = self.verts[i+1][1]
            cx += (x_i+x_ip1) * (x_i*y_ip1 - x_ip1*y_i)
            cy += (y_i+y_ip1) * (x_i*y_ip1 - x_ip1*y_i)

        return np.array([cx/(6.0*area),cy/(6.0*area)])

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
        pt = pt - self.offset

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
    def __init__ (self,circle_list,polygon_list,template):
        # save circles and polygons
        self.circles = circle_list
        self.polygons = polygon_list
        self.num_circles = len(circle_list)
        self.num_polygons = len(polygon_list)
        self.template = template

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
                #use the iith (diagonal entry to test border crossing detection)
                if i==j and not all_objects[i].isIn_template(self.template): 
                    #print 'cutting template boundary', all_objects[i].isIn_template(self.template), 'index=', i, '\n'
                    collision_mat[i,i] = 1
                elif  i!=j and all_objects[i].ifCollide(all_objects[j]):
                    collision_mat[i,j] = 1

        temp = np.sum(collision_mat,axis=0)
        for i in range(num_obj):
            if temp[i] == 0:
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

    print 'debugging centroid distance.........'
    print 'distance from circle 1 to circle 2 is: ', circ1.centroid_dist(circ2)
    print 'distance from circle 1 to poly 1 is: ', circ1.centroid_dist(poly1)
    print 'distance from poly 1 to circle 1 is: ', poly1.centroid_dist(circ1)
    print 'distance from poly 1 to poly 2 is: ', poly1.centroid_dist(poly2)

    #formulate a list and do all pre-processing (finding distance + collision, etc)
    template = np.array([-10.0,10.0,-10.0,10.0])
    item_lists = object_lists(circ_list,poly_list,template)
    item_lists.print_items_info()
    #print item_lists.poly_verts

    delta_x = np.array([[1,0],[1,0],[1,0],[1,0],[1,0],[1,0]],dtype='float')

    item_lists.update_delta_pos(delta_x)
    item_lists.print_items_info()

    #print item_lists.poly_verts

    #debug function
    #print if_intersect([0,0],[2,0],[1,1],[1,-1])
    #print if_intersect([0,0],[2,0],[3,1],[3,-1])
    print 'debugging collision resolving function...'
    c3 = copy.deepcopy(circ1)
    c4 = copy.deepcopy(circ2)
    poly3 = copy.deepcopy(poly1)
    c3.print_info()
    c4.print_info()
    poly3.print_info()

    print 'testing untangle circle collision..................'
    c3.untangle_collision_circle(c4,opt='self')
    c3.print_info()
    c4.print_info()

    print 'testing untangle circle polygon collision.............'
    c4.untangle_collision_polygon(poly3)
    c4.print_info()
    poly3.print_info()

    print 'testing untangle polygon collision....................'
    poly_p1 = obj_polygon([(0,0),(3,0),(3,2),(0,2),(0,0.0)],np.array([0.0,0.0]))
    poly_p2 = obj_polygon([(0,-0.1),(3,-0.1),(3,2.1),(0,2.1),(0,-0.1)],np.array([2.0,0.0]))
    poly_p1.print_info()
    poly_p2.print_info()
    print poly_p2.isIn_poly([3,0])
    poly_p1.untangle_collision_polygon(poly_p2)
    poly_p1.print_info()

    
    
    


