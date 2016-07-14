from objects import *  
import scipy.io as sio

def from_file(filename):
    #load in file
    shape_contents = sio.loadmat(filename)
    #print shape_contents
    wings_location = shape_contents['wings_location']
    
    poly_list= []
    for i in range(len(wings_location)):
        xpos = wings_location[i][0][0]
        ypos = wings_location[i][0][1]
        offset = np.array([xpos[0],ypos[0]])
        verts =[]
        for j in range(len(xpos)):
            verts.append((xpos[j]-xpos[0],ypos[j]-ypos[0]))

        poly_list.append(obj_polygon(verts,offset))

    #form arrays of objects
    circ_list= []

    #load in template information
    #template_info = shape_contents['template_info']
    template_info =  shape_contents['template_info']
    temp_dim = float(template_info[0])
    #print temp_dim
    template = np.array([-temp_dim, temp_dim, -temp_dim, temp_dim])

    template_obs_list = [] 
    radius = 1.0
    for i in range(1,len(template_info)):
        #print template_info[i][0][0][0]
        template_obs_list.append(obj_circle(np.array([template_info[i][0][0][0], \
                                                      template_info[i][0][1][0]],'float'),radius))
    

    item_lists = object_lists(circ_list,poly_list,template)


    return circ_list, poly_list, item_lists,template,template_obs_list

def save_to_file(filename,poly_lists,indices):
    offsets =[]
    for i in range(len(poly_lists)):
        offsets.append(poly_lists[i].offset)
    sio.savemat(filename,{'offsets':offsets,'indices':indices})
        


def median_input_set():
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
    template = np.array([-10.0,10.0,-10.0,10.0])
    item_lists = object_lists(circ_list,poly_list,template)

    print 'total item number= ', item_lists.num_circles + item_lists.num_polygons, 'total area = ',  item_lists.total_area()
    return circ_list, poly_list, item_lists,template

def large_input_set():
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
    circ16 = obj_circle(np.array([0.0,0.0]),1.8)
    circ17 = obj_circle(np.array([0.0,0.0]),1.6)
    circ18 = obj_circle(np.array([0.0,0.0]),1.4)
    circ19 = obj_circle(np.array([0.0,0.0]),1.2)
    circ20 = obj_circle(np.array([0.0,0.0]),1.0)
    circ21 = obj_circle(np.array([0.0,0.0]),0.8)
    circ22 = obj_circle(np.array([0.0,0.0]),3.0)
    circ23 = obj_circle(np.array([0.0,0.0]),3.5)
    circ24 = obj_circle(np.array([0.0,0.0]),0.4)
    circ25 = obj_circle(np.array([0.0,0.0]),0.3)

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

    poly16 = obj_polygon([(0,0),(1,0),(1,11),(0,11),(0,0)],np.array([0.0,0.0]))
    poly17 = obj_polygon([(0,0),(11,0),(11,1),(0,1),(0,0)],np.array([0.0,0.0]))
    poly18 = obj_polygon([(0,0),(1,0),(1,1),(0,1),(0,0)],np.array([0.0,0.0]))
    poly19 = obj_polygon([(0,0),(2,0),(2,2),(0,2),(0,0)],np.array([0.0,0.0]))
    poly20 = obj_polygon([(0,0),(3,0),(3,3),(0,3),(0,0)],np.array([0.0,0.0]))
    poly21 = obj_polygon([(0,0),(4,0),(4,4),(0,4),(0,0)],np.array([0.0,0.0]))
    poly22 = obj_polygon([(0,0),(10,0),(15,1),(0,0)],np.array([0.0,0.0]))
    poly23 = obj_polygon([(0,0),(0,7),(-3,10),(0,0)],np.array([0.0,0.0]))
    poly24 = obj_polygon([(0,0),(0,4),(-1,6),(0,0)],np.array([0.0,0.0]))
    poly25 = obj_polygon([(0,0),(1,0),(1,6),(0,0)],np.array([0.0,0.0]))

    #form arrays of objects
    circ_list= [circ1,circ2,circ3,circ4,circ5,circ6,circ7,circ8,\
                circ9,circ10,circ11,circ12,circ13,circ14,circ15,\
                circ16,circ17,circ18,circ19,circ20,circ21,circ22,\
                circ23,circ24,circ25]
    poly_list= [poly1,poly2,poly3,poly4,poly5,poly6,poly7,poly8,\
                poly9,poly10,poly11,poly12,poly13,poly14,poly15,\
                poly16,poly17,poly18,poly19,poly20,poly21,poly22,\
                poly23,poly24,poly25]
    template = np.array([-10.0,10.0,-10.0,10.0])
    item_lists = object_lists(circ_list,poly_list,template)
   

    print 'total item number= ', item_lists.num_circles + item_lists.num_polygons, 'total area = ',  item_lists.total_area()

    return circ_list, poly_list, item_lists,template

if __name__ == '__main__':
    #circ_list,poly_list, item_list,template = median_input_set()
    #circ_list,poly_list, item_list,template = large_input_set()
    circ_list,poly_list, item_list,template,template_obs_list = from_file('tiling_test.mat')
