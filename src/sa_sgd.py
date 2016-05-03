from input import *
#from SDG import *
from fast_SGD import *
from animation import *
import copy

# Simulated Annealing with incorporating Stochastic Gradient Descent

# step 1. Propose subset with area smaller than template
# step 2. run SGD to check whether it is feasible
#           if False, propose another subset randomly
# step 3. compare the area and replace

def smaller_area_subset(object_list, target_area):
    # save index only
    whole_index = range(len(object_list)) #index 0 ~ len-1
    subset_index = []
    subset_area = 0.

    while True:
        index = np.random.choice(whole_index)
        if index not in subset_index:
            tmp_area = subset_area + object_list[index].area
            if tmp_area >= target_area:
                # print subset_area
                return subset_index
            # else, it's not yet over the target area
            # then add the object
            subset_index.append(index)
            subset_area+=object_list[index].area

def propose_subset(target_area, input_object_list):
    converge = False
    while not converge:
        # -------------------
        # propose subset
        # -------------------
        subset_index = smaller_area_subset(input_object_list, target_area)
        subset_circle = []
        subset_poly = []
        for index in subset_index:
            object = input_object_list[index]
            if isinstance(object, obj_circle):
                subset_circle.append(object)
            elif isinstance(object, obj_polygon):
                subset_poly.append(object)
        # print subset_index
        # print len(subset_index)

        # ---------------------------
        # run SGD to see if feasible
        # ---------------------------
        item_lists = object_lists(subset_circle,subset_poly,template)
        animator = Animator(20,20)
        animator.add_circular_objects(item_lists.circ_diameter,\
                                  item_lists.circ_position,\
                                  item_lists.circ_collision,0.1)

                
        animator.add_polygon_objects(item_lists.poly_verts,\
                                 item_lists.poly_collision,0.1)

        propose_circle, propose_poly, area, converge, potential_vec = \
        SDG_tiling(subset_circle, subset_poly, template,30,animator,item_lists)

        matplotlib.pyplot.close("all")

    return propose_circle, propose_poly

def sum_total_area(circle_list, poly_list):
    l = circle_list + poly_list
    area = 0.
    for figure in l:
        area += figure.area
    return area

def simulated_anealing_sgd(template):
    # -------------------
    # create input
    # -------------------
    input_circle, input_poly = input_generator(template)
    input_object_list = input_circle + input_poly

    # -----------------------------------------
    # Propose new subset that AREA & FEASIBLE
    # -----------------------------------------
    target_area = (template[1]-template[0])*(template[3]-template[2])
    current_circle, current_poly = propose_subset(target_area, input_object_list)
    current_area = sum_total_area(current_circle, current_poly)

    # -----------------------
    # keep Best Result
    # -----------------------
    best_area = current_area
    best_circle = copy.deepcopy(current_circle)
    best_poly = copy.deepcopy(current_poly)

    # ------------------------------------------
    # now we have propose objects, let's run SA
    # ------------------------------------------
    T = np.arange(3, 0.1, -0.5)
    sa_iteration = 20

    for iteration in range(sa_iteration):
        for temperature in T:
            # -------------------
            # compare the area
            # -------------------
            propose_circle, propose_poly = propose_subset(target_area, input_object_list)
            propose_area = sum_total_area(propose_circle, propose_poly)

            if current_area < propose_area: # always accept!
                current_area = propose_area
                current_circle = copy.deepcopy(propose_circle)
                current_poly = copy.deepcopy(propose_poly)

            # if proposed area is bad, accept bases on probability
            else:
                delta_area = propose_area - current_area # < 0
                if np.random.rand() < np.exp(delta_area/temperature):
                    current_area = propose_area
                    current_circle = copy.deepcopy(propose_circle)
                    current_poly = copy.deepcopy(propose_poly)

            # to keep the best solution
            if propose_area > best_area:
                best_area = propose_area
                best_circle = copy.deepcopy(propose_circle)
                best_poly = copy.deepcopy(propose_poly)

    # return best_area, best_circle, best_poly
    print "# ============================"
    print best_area, 'circle:', len(best_circle), 'poly:', len(best_poly)
    print "# ============================"

    item_lists = object_lists(best_circle,best_poly,template)

    animator = Animator(20,20)
    animator.add_circular_objects(item_lists.circ_diameter,\
                                  item_lists.circ_position,\
                                  item_lists.circ_collision,0.1)

    #for i in range(len(item_lists.poly_collision)):
    #    item_lists.poly_collision[i] = False
    animator.add_polygon_objects(item_lists.poly_verts,\
                                  item_lists.poly_collision,100)




if __name__ == '__main__':
    template = np.array([-10.0,10.0,-10.0,10.0])
    simulated_anealing_sgd(template)
