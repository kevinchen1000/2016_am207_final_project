#script to do post analysis
import pickle
from matplotlib import pyplot as plt
from generate_input_r1 import *

#driver post-analysis function
if __name__ == "__main__":

    #generate circle and polygon set
    circ_list,poly_list,item_lists,template = large_input_set()
    obj_list = circ_list + poly_list

    #compute shape information for post_analysis
    obj_area_vec = np.zeros(len(obj_list),dtype = 'float')
    obj_area_ratio = np.zeros(len(obj_list),dtype = 'float')
    l_w_ratio = np.zeros(len(obj_list),dtype = 'float')
    for i in range(len(obj_list)):
        print obj_list[i].area
        obj_area_vec[i] = obj_list[i].area
        obj_area_ratio[i]= obj_list[i].area / obj_list[i].get_height() / obj_list[i].get_width()
        if obj_list[i].get_height() > obj_list[i].get_width():
            l_w_ratio[i] = obj_list[i].get_height() / obj_list[i].get_width()
        else:
            l_w_ratio[i] = obj_list[i].get_width() / obj_list[i].get_height()
        
        


    #load in data
    f = open('data_kc.p','rb')
    test_index = pickle.load(f)
    test_area = pickle.load(f)
    f.close()

    print test_index, test_area
    print type(test_index)

    #1. compute probability of item i being chosen
    num_objects = len(circ_list) + len(poly_list)
    num_tries = test_area.shape[0]
    num_chosen =  np.sum(test_index,axis=0)
    prob_chosen = num_chosen/(num_tries+0.0)
    num_items_incl = np.sum(test_index,axis=1)

    #2. compute expectation area contributed by an item
    E_area_chosen =  np.zeros(num_objects,dtype='float')
    E_area_not_chosen = np.zeros(num_objects,dtype='float')

    for i in range(num_tries):
        for j in range(num_objects):
            E_area_chosen[j] += test_index[i,j] * test_area[i]
            E_area_not_chosen[j] += (1-test_index[i,j]) * test_area[i]

            
    #print E_area_chosen.shape, num_chosen.shape
    E_area_chosen = E_area_chosen / (num_chosen+0.0)
    E_area_not_chosen = E_area_not_chosen / ((num_tries - num_chosen)+0.0)

    #plotting###################################################
    #1. bar graph of P(chosen)
    fig1 = plt.figure()
    plt.subplot(1,3,1)
    bins = np.arange(1,51)
    plt.title('chosen probability')
    plt.bar(bins,prob_chosen,label='chosen probability',color='green',alpha=0.5)
    plt.xlabel('item index')
    plt.ylabel('probability')
    plt.legend()

    plt.subplot(1,3,2)
    bins = np.arange(1,51)
    plt.title('expected area given included')
    plt.bar(bins,E_area_chosen,label='E_area_given_selected',color='green',alpha=0.5)
    plt.xlabel('item index')
    plt.ylabel('area')
    plt.legend()

    plt.subplot(1,3,3)
    bins = np.arange(1,51)
    plt.title('expected area given not included')
    plt.bar(bins,E_area_not_chosen,label='E_area_given_notselected',color='green',alpha=0.5)
    plt.xlabel('item index')
    plt.ylabel('area')
    plt.legend()
    #plt.show()

    fig2 = plt.figure()
    bins = np.arange(220,280,5)
    plt.title('distribution of area')
    plt.hist(test_area,bins,label='area distribution',color='green',alpha=0.5,normed="True")
    plt.xlabel('area')
    plt.ylabel('probability')
    plt.legend()

    fig3 = plt.figure()
    plt.subplot(3,3,1)
    plt.plot(obj_area_vec,prob_chosen,'go',label='object area vs prob chosen')
    plt.xlabel('object area')
    plt.ylabel('prob chosen')
    plt.legend()

    plt.subplot(3,3,2)
    plt.plot(obj_area_ratio,prob_chosen,'go',label='object area ratio vs prob chosen')
    plt.xlabel('object area ratio')
    plt.ylabel('prob chosen')
    plt.legend()

    plt.subplot(3,3,3)
    plt.plot(l_w_ratio,prob_chosen,'go',label='aspect ratio vs prob chosen')
    plt.xlabel('object aspect ratio')
    plt.ylabel('prob chosen')
    plt.legend()

    plt.subplot(3,3,4)
    plt.plot(obj_area_vec,E_area_chosen,'go',label='object area vs E_area_chosen')
    plt.xlabel('object area')
    plt.ylabel('E_area_chosen')
    plt.legend()

    plt.subplot(3,3,5)
    plt.plot(obj_area_ratio,E_area_chosen,'go',label='object area ratio vs E_area_chosen')
    plt.xlabel('object area raio')
    plt.ylabel('E_area_chosen')
    plt.legend()

    plt.subplot(3,3,6)
    plt.plot(l_w_ratio,E_area_chosen,'go',label='aspect ratio vs E_area_chosen')
    plt.xlabel('object aspect ratio')
    plt.ylabel('E_area_chosen')
    plt.legend()

    plt.subplot(3,3,7)
    plt.plot(obj_area_vec,E_area_not_chosen,'go',label='object area vs E_area_not_chosen')
    plt.xlabel('object area')
    plt.ylabel('E_area_not_chosen')
    plt.legend()

    plt.subplot(3,3,8)
    plt.plot(obj_area_ratio,E_area_not_chosen,'go',label='object area ratio vs E_area_not_chosen')
    plt.xlabel('object aspect ratio')
    plt.ylabel('E_area_not_chosen')
    plt.legend()

    plt.subplot(3,3,9)
    plt.plot(l_w_ratio,E_area_not_chosen,'go',label='aspect ratio vs E_area_not_chosen')
    plt.xlabel('object aspect ratio')
    plt.ylabel('E_area_not_chosen')
    plt.legend()

    plt.show()

    
    #3. Baysian inference
    
    



    



    

    
