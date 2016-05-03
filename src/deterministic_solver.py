#deterministic solver 
#implemented based on method discussed in http://cgi.csc.liv.ac.uk/~epa/surveyhtml.html
from objects import *
from animation import *
import numpy as np
from operator import add
#Strip packing First-Fit Decreasing Height algorithm
#circ_list : array of circle objects
#poly_list : array of polygon objects
#template : numpy array of [xmin,xmax,ymin,ymax]
def deterministic_FFDH(circ_list,poly_list,template):
    obj_list = circ_list + poly_list
    num_objects = len(circ_list)+len(poly_list)

    xmin = template[0]; xmax = template[1]
    ymin = template[2]; ymax = template[3]

    #measure height in descending order
    height_vec = np.zeros(num_objects,dtype ='float')
    width_vec = np.zeros(num_objects,dtype ='float')

    for i in range(num_objects):
        height_vec[i] = obj_list[i].get_height()
        width_vec[i] = obj_list[i].get_width()

    height_ind_vec = np.argsort(height_vec)[::-1]

    #initialize a list of list of tiling origins
    tile_list = [[]]
    #tile_width = []
    #tile in the objects in sequence
    num_include_obj = 0
    ind_to_include =[]
    for i in range(num_objects):
        #if can be put into the template
        if fit_tile(obj_list[height_ind_vec[i]],tile_list,obj_list):
            update_tile_list(tile_list,obj_list[height_ind_vec[i]],height_ind_vec[i],obj_list)
            ind_to_include.append(height_ind_vec[i])
            num_include_obj +=1

    #compute total area
    tot_area = 0
    for i in range(len(ind_to_include)):
        tot_area +=obj_list[i].area

    print 'total area= ',tot_area , 'with indices', ind_to_include

    #assemble the list of tiled objects
    tile_obj = []
    tile_circ = []
    tile_poly = []
    for i in range(len(tile_list)):
        for j in range(len(tile_list[i])):
            tile_obj.append(obj_list[tile_list[i][j][1]])
            if isinstance(tile_obj[-1],obj_polygon):
                tile_obj[-1].offset = tile_list[i][j][0] + np.array([-tile_obj[-1].bounding_box[0],-tile_obj[-1].bounding_box[2]])
                tile_poly.append(tile_obj[-1])
            else:
                tile_obj[-1].pos = np.array(map(add,tile_list[i][j][0],[tile_obj[-1].radius, tile_obj[-1].radius]))
                tile_circ.append(tile_obj[-1])


    #return the list of tiled objects
    return tile_circ,tile_poly
        
'''update tile_list locally '''
def update_local_tile_list(tile_list,row_ind,xpos,ypos,obj_ind,obj):
    if len(tile_list) == row_ind:
        tile_list.append([])

    print 'current row ind = ,',row_ind

    if isinstance(obj,obj_polygon):
        tile_list[row_ind].append([[xpos,ypos],obj_ind])
    else:
        tile_list[row_ind].append([[xpos+0*obj.radius,ypos+0*obj.radius],obj_ind])
    

'''function that update tile_list'''
#assume that the object can be added
def update_tile_list(tile_list,obj,obj_ind,obj_list,xmin=-10.0,xmax=10.0,ymin=-10.0,ymax=10.0):
    #define a small offset
    eps = 1e-3

    #empty tile_list
    if len(tile_list[0]) == 0:
        update_local_tile_list(tile_list,0,xmin+eps,ymin+eps,obj_ind,obj)
        print 'to be tiled at the position,', xmin,ymin
        return

    #non-empty list
    for i in range(len(tile_list)):
        ypos =  tile_list[i][-1][0][1] + 0*obj_list[tile_list[i][-1][1]].get_height()
        xpos =  tile_list[i][-1][0][0] + obj_list[tile_list[i][-1][1]].get_width() 
        if (ymax -ypos)> obj.get_height() and (xmax - xpos) > obj.get_width():
            update_local_tile_list(tile_list,i,xpos+eps,tile_list[i][-1][0][1]+eps,obj_ind,obj)
            print 'to be tiled at the position,', xpos,ypos
            return

    #add a new row
    if (ymax - (tile_list[-1][0][0][1] + obj_list[tile_list[-1][0][1]].get_height())> obj.get_height()) and \
       (xmax - xmin) > obj.get_width():
        ypos = (tile_list[-1][0][0][1] + obj_list[tile_list[-1][0][1]].get_height())
        update_local_tile_list(tile_list,len(tile_list),xmin+eps,ypos+eps,obj_ind,obj)
        print 'to be tiled at the position,', xmin,ypos
        return 

    #else return an error
    assert(0)
                 

    

'''function that update tile_width'''
#assume that the object can be added
#def update_tile_width(tile_width,obj,obj_ind,xmin=template[0],xmax=template[1],ymin=template[2],ymax=template[3]):

'''function that checks whether an object can fit the current template'''
def fit_tile(obj,tile_list,obj_list,xmin=-10.0,xmax=10.0,ymin=-10.0,ymax=10.0):
    fit = False

    if len(tile_list[0]) == 0:
        print 'first tiled object'
        return (ymax - ymin) > obj.get_height() and  (xmax - xmin) > obj.get_width()

    #check each row
    for i in range(len(tile_list)):
        ypos =  tile_list[i][-1][0][1] + 0*obj_list[tile_list[i][-1][1]].get_height()
        xpos =  tile_list[i][-1][0][0] + obj_list[tile_list[i][-1][1]].get_width() 
        if (ymax - ypos) > obj.get_height() and \
           (xmax - xpos) > obj.get_width() :
            print 'tiled after row', i
            return True

    #whether can create a new row
    if (ymax - (tile_list[-1][0][0][1] + obj_list[tile_list[-1][0][1]].get_height())> obj.get_height()) and \
       (xmax - xmin) > obj.get_width():
        print 'can fit into a new row'
        return True

    return fit

if __name__  == '__main__' :
    print 'start testing deterministic solver...'
        

    #initialize animator
    animator = Animator(20,20)

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

    poly16 = obj_polygon([(0,0),(2,0),(2,15),(0,15),(0,0)],np.array([0.0,0.0]))

    #form arrays of objects
    circ_list= [circ1,circ2,circ3,circ4,circ5,circ6,circ7,circ8,\
                circ9,circ10,circ11,circ12,circ13,circ14,circ15]
    poly_list= [poly1,poly2,poly3,poly4,poly5,poly6,poly7,poly8,\
                poly9,poly10,poly11,poly12,poly13,poly14,poly15,\
                poly16]

    #formulate a list and do all pre-processing (finding distance + collision, etc)
    template = np.array([-10.0,10.0,-10.0,10.0])

    #use deterministic methods to choose the set of items to include 
    incl_circ, incl_poly = deterministic_FFDH(circ_list,poly_list,template)

    #detect collision etc
    print len(incl_circ)
    for i in range(len(incl_circ)):
        incl_circ[i].print_info()

    print len(incl_poly)
    for i in range(len(incl_poly)):
        incl_poly[i].print_info()

    item_lists = object_lists(incl_circ,incl_poly,template)

    #add to animator
    #for i in range(len(item_lists.circ_collision)):
    #    item_lists.circ_collision[i] = False
    animator.add_circular_objects(item_lists.circ_diameter,\
                                  item_lists.circ_position,\
                                  item_lists.circ_collision,0.1)

    #for i in range(len(item_lists.poly_collision)):
    #    item_lists.poly_collision[i] = False
    #animator.show_title(42.2544,100)
    animator.add_polygon_objects(item_lists.poly_verts,\
                                  item_lists.poly_collision,50)
    
