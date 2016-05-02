# Simulated Annealing Solver
from objects import *
from animation import *
import numpy as np
import random
from operator import add
from matplotlib import pyplot
import copy

# Initialized 15 circle objects 
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

# Initialize 15 polygon objects
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

circ_list= [circ1,circ2,circ3,circ4,circ5,circ6,circ7,circ8,\
            circ9,circ10,circ11,circ12,circ13,circ14,circ15]
poly_list= [poly1,poly2,poly3,poly4,poly5,poly6,poly7,poly8,\
            poly9,poly10,poly11,poly12,poly13,poly14,poly15]

def index_to_list(circ_list, poly_list, circ_idx, poly_idx):
    result_circ_list = [copy.deepcopy(circ_list[i])  for i in range(len(circ_idx)) if circ_idx[i] == 1]
    result_poly_list = [copy.deepcopy(poly_list[i])  for i in range(len(poly_idx)) if poly_idx[i] == 1]
    return result_circ_list, result_poly_list

def get_value(circ_list, poly_list, circ_soln, poly_soln):
    score = 0
    count = 0
    for i in range(len(circ_soln)):
        if circ_soln[i] == 1:
            score = score + circ_list[i].area
            count += 1
    for j in range(len(poly_soln)):
        if poly_soln[j] == 1:
            count += 1
            score = score + poly_list[j].area
    print count, "items"
    return score

def tile_item(circ_list, poly_list, circ_soln, poly_soln, index, circle, template, max_iter):
    success = True
    if circle:
        for k in range(max_iter):
            temp_circle = copy.deepcopy(circ_list[index])
            new_pos = np.random.uniform(-10., 10., size = 2)
            temp_circle.set_pos(new_pos)

            for i in range(len(circ_soln)):
                if(circ_soln[i] == 1):
                    if temp_circle.ifCollide_circle(circ_list[i]):
                        success = False
                        break
            for i in range(len(poly_soln)):
                if(poly_soln[i] == 1):
                    if temp_circle.ifCollide_polygon(poly_list[i]):
                        success = False
                        break
            if not temp_circle.isIn_template(template):
                success = False
            if success:
                return success, new_pos
    if not circle:
        for k in range(max_iter):
            temp_poly = copy.deepcopy(poly_list[index])
            new_pos = np.random.uniform(-10., 10., size = 2)
            temp_poly.set_pos(new_pos)

            for i in range(len(poly_soln)):
                if(poly_soln[i] == 1):
                    if temp_poly.ifCollide_polygon(poly_list[i]):
                        success = False
                        break
            for i in range(len(circ_soln)):
                if(circ_soln[i] == 1):
                    if temp_poly.ifCollide_circle(circ_list[i]):
                        success = False
                        break
            if not temp_poly.isIn_template(template):
                success = False
            if success:
                return success, new_pos
    return success, new_pos

def random_item(circ_list, poly_list, init_circ, init_poly,template, step):

        
    len_circ = len(init_circ)
    len_poly = len(init_poly)
    max_len = len_circ + len_poly
    
    new_circ = copy.deepcopy(init_circ)
    new_poly = copy.deepcopy(init_poly)
    
    rand_idx = np.random.choice(max_len, size = step, replace = False)
    rand_circ = []
    rand_poly = []
    for i in rand_idx:
        if i < len_circ:
            rand_circ.append(i)
        else:
            rand_poly.append(i - len_circ)
    
    for i in rand_circ:
        if new_circ[i] == 1:
            new_circ[i] = 0
        if new_circ[i] == 0:
            success, pos = tile_item(circ_list, poly_list, new_circ, new_poly, i, True,template, 1000)
            if success == False:
                pass
            else:
                circ_list[i].set_pos(pos)
                new_circ[i] = 1
                
                
    for i in rand_poly:
        if new_poly[i] == 1:
            new_poly[i] = 0
        if new_poly[i] == 0:
            success, pos = tile_item(circ_list, poly_list, new_circ, new_poly, i, False,template, 500)
            if success == False:
                pass
            else:
                poly_list[i].set_pos(pos)
                new_poly[i] = 1
                
    return new_circ, new_poly

def simulated_annealing(circ_list, poly_list, templete, init_circ, init_poly,init_temp, thermostat, reannealing, itol):
    # number of accepted solution
    step = 1
    accepted = 0 
    hist=[]
    num_iter = 0
    T = init_temp
    best_circ, best_poly = init_circ, init_poly
    best_value = get_value(circ_list, poly_list, init_circ, init_poly)
    
    current_circ = init_circ
    current_poly = init_poly
    current_value = get_value(circ_list, poly_list, current_circ, current_poly)
    
    while True:
        num_iter += 1
        hist.append(current_value)
        # Generate a new proposed solution
        propose_circ, propose_poly  = random_item(circ_list, poly_list, current_circ, current_poly, template, step)
        propose_value = get_value(circ_list, poly_list,propose_circ, propose_poly)
        
        delta = propose_value - current_value
        
        # if not accept
        if np.random.random() > np.exp(delta / T):
            continue
        else:
            accepted += 1
            current_circ, current_poly = propose_circ, propose_poly
            current_value = get_value(circ_list, poly_list, current_circ, current_poly)
        if  current_value > best_value:
            best_circ, best_poly, best_value = current_circ, current_poly, current_value  
        print num_iter, current_value
        # check if it is time to cool down
        if num_iter % reannealing == 0:
            T = thermostat * T;
            #temperature =  temperature/np.log(it)
            
            #if we get too cold, reheat
            if T < 0.01:
                T = 20
        
        #stopping criteria
        if num_iter > itol:
            print 'itol'
            break   
    #plt.plot(hist)
    #plt.show()
    return current_circ, current_poly, best_value, num_iter, accepted, hist

x = np.zeros([15])
y = np.zeros([15])
template = np.array([-10.0,10.0,-10.0,10.0])
c, p, v, n, accpet, hist = simulated_annealing(circ_list, poly_list, template, x, y, 30, 0.8, 50, 2000)

plt.plot(hist)
plt.savefig("hist.png")

a, b = index_to_list(circ_list, poly_list, c, p)
animator = Animator(20,20)
template = np.array([-10.0,10.0,-10.0,10.0])

item_lists = object_lists(a,b,template)
#print item_lists.collision_vec


animator.add_circular_objects(item_lists.circ_diameter,\
                                  item_lists.circ_position,\
                                  item_lists.circ_collision,0.1)


animator.add_polygon_objects(item_lists.poly_verts,\
                                  item_lists.poly_collision,100)

