#optimization that combines SA and SGD

from sa_sgd import *
import pickle
from generate_input_r1 import *


#driver function
if __name__ == '__main__':
    
    #generate circle and polygon set
    circ_list,poly_list,item_lists,template = large_input_set()

    #prepare graphics
    animator = Animator(20,20)

    animator.add_circular_objects(item_lists.circ_diameter,\
                                  item_lists.circ_position,\
                                  item_lists.circ_collision,0.1)

                
    animator.add_polygon_objects(item_lists.poly_verts,\
                                 item_lists.poly_collision,0.1)

    animator.show_title(50,0.4)


    #solve num_tries times
    num_tries = 100
    index_mat = np.zeros((num_tries,len(circ_list)+len(poly_list)), dtype ='int')
    area_vec = np.zeros(num_tries,dtype='float')
    
    
    
    for i in range(num_tries):
        incl_circ_list, incl_poly_list, area, converge,potential_vec,index_vec = SDG_tiling(circ_list,poly_list,template,100,animator,item_lists)
        index_mat[i,:] = index_vec
        area_vec[i] = area
        print 'area=', area, '...', index_vec
    
    f = open('data.p','wb')
    pickle.dump(index_mat,f)
    pickle.dump(area_vec,f)
    f.close()

    f = open('data.p','rb')
    test_index = pickle.load(f)
    test_area = pickle.load(f)
    f.close()

    print test_index, test_area
