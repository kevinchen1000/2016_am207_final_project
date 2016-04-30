#stochastic gradient descent solver

from objects import *
from animation import *
import numpy as np
import random

#Stochastic Gradient Descent to sequentially tile a template algorithm
#circ_list : array of circle objects
#poly_list : array of polygon objects
#template : numpy array of [xmin,xmax,ymin,ymax]
def SDG_tiling(circ_list,poly_list,template,animator = None,item_lists =None):
    obj_list = circ_list + poly_list
    num_objects = len(circ_list)+len(poly_list)

    xmin = template[0]; xmax = template[1]
    ymin = template[2]; ymax = template[3]

    #initialize arrangement at T = 0  
    #-->enforce no collision but may cross boundary of template (infeasible)

    num_iter =100
    iterate = 0
    tol = 1e-2
    converge = False
    
    #initialize centroid positions
    centroid_pos = np.zeros((num_objects,2) ,dtype = 'float')
    initialize_tiling_positions(obj_list,xmin,xmax,ymin,ymax)

    #initialize to screen
    if animator is not None:
        init_pos =  obtain_centroid_pos(obj_list)
        item_lists.update_delta_pos(init_pos)

    #optimization step
    while iterate < num_iter and not converge:

        #sequential movement of each object
        for i in range(len(obj_list)):
            #plot update for debugging
            delta_x = SDG_update(obj_list,i,xmin,xmax,ymin,ymax)
            if animator is not None:
                delta_vec = np.zeros((num_objects,2),dtype = 'float')
                delta_vec[i,:] = delta_x
                item_lists.update_delta_pos(delta_vec)
                animator.update_circular_objects(item_lists.circ_position,item_lists.circ_collision)
                animator.update_polygon_objects(item_lists.poly_verts,item_lists.poly_collision)
                #for debug
                print 'current center positions =' ,  obtain_centroid_pos(obj_list)
                variable = raw_input('input something!: ')

               
                #animator.remove_circles()
                #animator.remove_polygons()
                
                

        #update counter
        iterate +=1

        #comparsion for convergence
        update_centroid_pos = obtain_centroid_pos(obj_list)

        print 'update_centroid_pos= ',update_centroid_pos
        converge = np.sum(np.sum((centroid_pos - update_centroid_pos)**2)) /(num_objects+0.0) <tol

        print 'local converge =', np.sum(np.sum((centroid_pos - update_centroid_pos)**2)) /(num_objects+0.0)
        centroid_pos = copy.deepcopy(update_centroid_pos)


    print centroid_pos
    print update_centroid_pos
    print 'error =',  np.sum(np.sum((centroid_pos - update_centroid_pos)**2)) /(num_objects+0.0)
    print 'return solution after ', iterate, 'iterations'
    #return the updated object list and convergence flag
    return circ_list, poly_list, converge

''' function that copies all centroid positions of the object into a numpy array through DEEPCOPY'''
def obtain_centroid_pos(obj_list):
    num_objects = len(obj_list)
    centroid_pos =  np.zeros((num_objects,2) ,dtype = 'float')
    for i in range(len(obj_list)):
        if isinstance(obj_list[i],obj_circle):
            centroid_pos[i,:] = copy.deepcopy(obj_list[i].pos)
        else: 
            centroid_pos[i,:] = copy.deepcopy(obj_list[i].centroid + obj_list[i].offset)
    return centroid_pos

''' function that restricts displacement amplitude to avoid collision'''
def restrict_disp(disp,obj_list,ind):
    disp_mag =  np.linalg.norm(disp)
    disp_dir = disp / disp_mag

    temp_obj = copy.deepcopy(obj_list[ind])
    
    #detect if there is collision if incremented by full disp
    temp_obj.increment_pos(disp)

    if if_collision_free(obj_list,ind,temp_obj):
        print 'ok to move entire distance, no collision'
        return disp

    #otherwise, need to find the closest point
    print 'not collision free if moved by total distance '
    num_iter = 20
    iterate = 1
    collision_free = False
    temp_obj.increment_pos(-disp)
    total_disp_mag =0

    #continuously half the displacement magnitude until no collision
    while iterate < num_iter:

        #update object position
        delta_mag = (0.5**iterate)*disp_mag
        temp_obj.increment_pos(delta_mag*disp_dir)

        if if_collision_free(obj_list,ind,temp_obj):
            total_disp_mag += delta_mag
        else:
            temp_obj.increment_pos(-delta_mag*disp_dir)
    
        iterate +=1

    return total_disp_mag * disp_dir

'''detect if the ind th object is collision free with the rest of obstacles'''
def if_collision_free(obj_list,ind,temp_obj):
    collision_free = True
    for i in range(len(obj_list)):
        if i != ind and temp_obj.ifCollide(obj_list[i]):
            collision_free = False
            return collision_free
    return collision_free

'''function that updates the center (or centroid) position of an object based on stochastic gradient descent'''
#ind indicates the ind_th object in the list (under the potential influence of other)
def SDG_update(obj_list,ind,xmin,xmax,ymin,ymax):

    #select the object to be updated
    obj = obj_list[ind]

    #randomly select a subset of items
    num_to_select = len(obj_list)-1
    select_array = np.concatenate([np.array(range(0,ind)),np.array(range(ind+1,len(obj_list)))])
    #print 'select array is: ', select_array
    #print 'num_to_select is: ', num_to_select
    neighbor_indices = random.sample(select_array,num_to_select)

    #define constants
    K = 0.1
    G = 0.1*0.1
    origin = np.array([xmin,ymin])

    #compute displacement due to global potenial
    disp = -K * (obj.get_position() - origin)
   
    #compute influence due to other selected objects
    for neighbor_ind in neighbor_indices:
        neighbor_ind = int(neighbor_ind)
        local_dir = obj.centroid_dir(obj_list[neighbor_ind])
        disp += G * obj_list[neighbor_ind].area / obj.centroid_dist(obj_list[neighbor_ind]) * local_dir

    print 'initial disp =', disp
    #restrict the magnitude of disp to avoid collision
    disp = restrict_disp(disp,obj_list,ind)

    print 'object [',ind, '] increments by', disp

    #update position of object
    #obj_list[ind].increment_pos(disp)
    
    return disp

'''initialize the tiling positions of objects that are collision free but may intersect template
   borders
'''
def initialize_tiling_positions(obj_list,xmin,xmax,ymin,ymax):
    mean = np.array([0.0,0.0])
    scale = 0.1
    cov = np.array([[(xmax-xmin)**2*scale,0.0],[0.0,scale*(ymax-ymin)**2]])

    #initialize the positions sequentially
    for i in range(len(obj_list)):
        init = np.random.multivariate_normal(mean,cov,1)
        #debug
        if i==0:
            init= np.array([[0.0,0.0]])+10
        elif i == 1:
            init= np.array([[1.0,-1.0]])+10

        init_x = init[0][0]; init_y = init[0][1]
        obj_list[i].set_pos(np.abs(np.array([init_x,init_y])) +np.array([xmin,ymin]))

        print 'initial positon for ', i, 'th object is at', np.abs(np.array([init_x,init_y])) +np.array([xmin,ymin])
        #detect and 
        for j in range(i):
            if obj_list[i].ifCollide(obj_list[j]):
                obj_list[i].untangle_collision(obj_list[j])
                   
    return obj_list



#driver script
if __name__ == '__main__':
    print 'stochastic gradient descent starts...'


   

    #initialize a few circles and polygon objects
    circ1 = obj_circle(np.array([0.0,0.0]),0.8)
    circ2 = obj_circle(np.array([0.0,0.0]),0.9)
    circ3 = obj_circle(np.array([0.0,0.0]),1.1)
    circ4 = obj_circle(np.array([0.0,0.0]),0.4)
    circ5 = obj_circle(np.array([0.0,0.0]),0.6)
    circ6 = obj_circle(np.array([0.0,0.0]),0.8)
    circ7 = obj_circle(np.array([0.0,0.0]),1.2)
    circ8 = obj_circle(np.array([0.0,0.0]),1.3)
    circ9 = obj_circle(np.array([0.0,0.0]),1.7)
    circ10 = obj_circle(np.array([0.0,0.0]),1.9)
    circ11 = obj_circle(np.array([0.0,0.0]),2.1)
    circ12 = obj_circle(np.array([0.0,0.0]),2.2)
    circ13 = obj_circle(np.array([0.0,0.0]),3.3)
    circ14 = obj_circle(np.array([0.0,0.0]),2.4)
    circ15 = obj_circle(np.array([0.0,0.0]),3.5)

    poly1 = obj_polygon([(0,0),(3,0),(3,2),(0,2),(0,0)],np.array([0.0,0.0]))
    poly2 = obj_polygon([(0,0),(3,6),(0,4),(0,0)],np.array([0.0,0.0]))
    poly3 = obj_polygon([(0,0),(2,0),(1,1),(0,0)],np.array([0.0,0.0]))

    poly4 = obj_polygon([(0,0),(4,0),(4,2),(0,2),(0,0)],np.array([0.0,0.0]))
    poly5 = obj_polygon([(0,0),(0,1),(-2,-1),(0,0)],np.array([0.0,0.0]))
    poly6 = obj_polygon([(0,0),(2,0),(3,2),(0,0)],np.array([0.0,0.0]))

    poly7 = obj_polygon([(0,0),(1,0),(1,6),(0,6),(0,0)],np.array([0.0,0.0]))
    poly8 = obj_polygon([(0,0),(3,0),(0,0.5),(0,0)],np.array([0.0,0.0]))
    poly9 = obj_polygon([(0,0),(2,0),(0,2),(0,0)],np.array([0.0,0.0]))

    poly10 = obj_polygon([(0,0),(3,0),(5,2),(4,3),(0,0)],np.array([0.0,0.0]))
    poly11 = obj_polygon([(0,0),(2,0),(3,1),(1,3),(-1,1),(0,0)],np.array([0.0,0.0]))
    poly12 = obj_polygon([(0,0),(2,0),(3,1),(1,3),(-1,1),(0,0)],np.array([0.0,0.0]))

    poly13 = obj_polygon([(0,0),(3,0),(3,2),(0,2),(0,0)],np.array([0.0,0.0]))
    poly14 = obj_polygon([(0,0),(2,0),(3,1.5),(2,3),(0,3),(-1,1.5),(0,0)],np.array([0.0,0.0]))
    poly15 = obj_polygon([(0,0),(2,0),(3,1.5),(2,3),(0,3),(-1,1.5),(0,0)],np.array([0.0,0.0]))

    #form arrays of objects
    circ_list= [circ1,circ2,circ3,circ4,circ5,circ6,circ7,circ8,\
                circ9,circ10,circ11,circ12,circ13,circ14,circ15]
    poly_list= [poly1,poly2,poly3,poly4,poly5,poly6,poly7,poly8,\
                poly9,poly10,poly11,poly12,poly13,poly14,poly15]

    #formulate a list and do all pre-processing (finding distance + collision, etc)
    template = np.array([-10.0,10.0,-10.0,10.0])

    #debug
    circ_list_short = [circ1]
    poly_list_short = [poly1]
    #use stochastic methods to tile all items and return the convergene flag
    #incl_circ, incl_poly, converge =  SDG_tiling(circ_list_short,poly_list_short,template)
    #print 'method has converged = ', converge

    #plot solution
    #item_lists = object_lists(incl_circ,incl_poly,template)
    #animator = Animator(20,20)
    #animator.add_circular_objects(item_lists.circ_diameter,\
    #                              item_lists.circ_position,\
    #                              item_lists.circ_collision,0.1)

    #for i in range(len(item_lists.poly_collision)):
    #    item_lists.poly_collision[i] = False
    #animator.add_polygon_objects(item_lists.poly_verts,\
    #                              item_lists.poly_collision,100)


    #use plotting version
    animator = Animator(20,20)
    item_lists = object_lists(circ_list_short,poly_list_short,template)
    animator.add_circular_objects(item_lists.circ_diameter,\
                                  item_lists.circ_position,\
                                  item_lists.circ_collision,0.1)

                
    animator.add_polygon_objects(item_lists.poly_verts,\
                                 item_lists.poly_collision,0.5)

    circ_list, poly_list, converge = SDG_tiling(circ_list_short,poly_list_short,template,animator,item_lists)

