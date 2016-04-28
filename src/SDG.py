#stochastic gradient descent solver

from objects import *
from animation import *
import numpy as np

#Stochastic Gradient Descent to sequentially tile a template algorithm
#circ_list : array of circle objects
#poly_list : array of polygon objects
#template : numpy array of [xmin,xmax,ymin,ymax]
def SDG_tiling(circ_list,poly_list,template):
    obj_list = circ_list + poly_list
    num_objects = len(circ_list)+len(poly_list)

    xmin = template[0]; xmax = template[1]
    ymin = template[2]; ymax = template[3]

    #initialize arrangement at T = 0  
    #-->enforce no collision but may cross boundary of template (infeasible)

    num_iter =100
    iter = 0
    tol = 1
    converge = False
    
    #initialize centroid positions
    centroid_pos = np.zeros((num_objects,2) ,dtype = 'float')

    #optimization step
    while iter < num_iter and not converge:

        #sequential movement of each object
        for i in range(len(obj_list)):
            SDG_update(obj_list,i,xmin,xmax,ymin,ymax)

        #update counter
        num_iter +=1

        #comparsion for convergence
        update_centroid_pos = obtain_centroid_pos(obj_list)
        converge = np.sum(np.sum((centroid_pos - update_centroid_pos)**2)) /(num_objects+0.0) <tol
        centroid_pos = update_centroid_pos

    #return the updated object list and convergence flag
    return obj_list, converge

''' function that copies all centroid positions of the object into a numpy array through DEEPCOPY'''
def obtain_centroid_pos(obj_list):
    centroid_pos =  np.zeros((num_objects,2) ,dtype = 'float')
    for i in range(len(obj_list)):
        if isinstance(obj_list[i],obj_circle):
            centroid_pos[i,:] = copy.deepcopy(obj_list[i].pos)
        else: 
            centroid_pos[i,:] = copy.deepcopy(obj_list[i].centroid)
    return centroid_pos

'''function that updates the center (or centroid) position of an object based on stochastic gradient descent'''
#ind indicates the ind_th object in the list (under the potential influence of other)
def SDG_update(obj_list,ind,xmin,xmax,ymin,ymax):
    return obj_list
