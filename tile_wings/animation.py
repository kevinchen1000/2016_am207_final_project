#script that contains plotting functions

import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
from matplotlib.collections import PolyCollection
import numpy as np
#from timer import Timer
import time
import matplotlib

from objects import *

#define class Animator

class Animator(object):
    def __init__(self, template_w, template_h):
        
        #set up figure
        plt.ion()
        #fig = plt.figure(figsize = (10,10))
        fig = plt.figure(figsize = (template_w,template_h))
        #print template_w, template_h
        #raw_input()
        ax = fig.gca()

        #set axis dimension 
        ax.axis( [-2*template_w/2.0, 2*template_w/2.0,-2*template_h/2.0,2*template_h/2.0]) 
        #ax.axis([-template_w/2.0, template_w/2.0,-template_h/2.0,template_h/2.0]) 
        ax.get_xaxis().set_visible(True)
        ax.get_yaxis().set_visible(True)

        #draw template
        self.plot_template(template_w,template_h)

        #draw to screen
        plt.draw()
        plt.pause(0.1)

        #save axes to object
        self.ax = ax

        #save animator canvas area
        self.width = template_w
        self.height = template_h
        self.title = plt.text(0,20,'start title',color = 'k')

    #plot template
    def plot_template(self,l = 20.0, w = 20.0):
        x_pos = [-l/2, l/2, l/2, -l/2, -l/2 ]
        y_pos = [-w/2, -w/2, w/2, w/2, -w/2]
        plt.plot(x_pos,y_pos,'k-',linewidth = 3.0)
        plt.xlabel('x')
        plt.ylabel('y')

    def generate_color_array(self,collision):
        #define colors (red = collision, blue = no collision)
        colors = [None] * len(collision)
        for i in range(len(collision)):
            if collision[i] == True:
                colors[i] = (1.0,0.0,0.0)
            else:
                colors[i] = (0.0,0.0,1.0)
        return colors


    #add circular objects
    def add_circular_objects(self,diameters,positions,collision,pause_time=0.1):
        
        circ_colors = self.generate_color_array(collision)
        self.circles = EllipseCollection(widths=diameters,
                                         heights=diameters,
                                         angles=np.zeros_like(diameters),
                                         units='xy',
                                         offsets=positions,
                                         transOffset=self.ax.transData,
                                         edgecolor='k', facecolor=circ_colors)


        self.ax.add_collection(self.circles)

        #add text label
        text_labels = [None] * len(collision)
        for i in range(len(collision)):
            text_labels[i]= plt.text(positions[i,0],positions[i,1],str(i),color = 'k')
        self.circle_labels = text_labels

        #draw to screen
        #plt.draw()
        #plt.pause(pause_time)

        #remove text labels
        #for i in range(len(collision)):
        #    text_labels[i].remove()

    #add polygon objects
    def add_polygon_objects(self,verts,collision,pause_time=0.1):

        poly_colors = self.generate_color_array(collision)
        self.polygons = PolyCollection(verts, facecolors=poly_colors)

        self.ax.add_collection(self.polygons)

        #add text label
        num_circles = self.circles.get_offsets().shape[0]
        #print 'number of circles =', num_circles
      
        text_labels = [None] * len(collision)
        for i in range(len(collision)):
            temp = np.array(verts[i])
            x = np.mean(temp[:,0])
            y = np.mean(temp[:,1])
        
            text_labels[i]= plt.text(x,y,str(i+num_circles),color = 'k')

        self.polygon_labels = text_labels


        plt.draw()
        plt.pause(pause_time)

        #remove text labels
        #for i in range(len(collision)):
        #    text_labels[i].remove()


    #remove circular objects:
    def remove_circles(self):
        self.circles.remove()
        for label in self.circle_labels:
            label.remove()

    #remove polygon objects:
    def remove_polygons(self):
        self.polygons.remove()
        for label in self.polygon_labels:
            label.remove()

    #update circular objects
    def update_circular_objects(self,positions,collision,pause_time = 0.1):
        
        #set circle colors
        circ_colors = self.generate_color_array(collision)
        self.circles.set_facecolors(circ_colors)

        #set circle positions
        self.circles.set_offsets(positions)

        #remove labels
        for label in self.circle_labels:
            label.remove()

        #add labels
        text_labels = [None] * len(collision)
        for i in range(len(collision)):
            text_labels[i]= plt.text(positions[i,0],positions[i,1],str(i),color = 'k')
        self.circle_labels = text_labels


        plt.draw()
        plt.pause(pause_time)


    #update polygon objects
    def update_polygon_objects(self,positions,collision,pause_time=0.1):

        #set polygon colors
        poly_colors = self.generate_color_array(collision)
        self.polygons.set_facecolors(poly_colors)

        #set polygon positions
        #print 'new verts positions=' , positions
        self.polygons.set_verts(positions)

        #remove labels
        for label in self.polygon_labels:
            label.remove()
        #self.title.remove()

        #add new labels
        num_circles = self.circles.get_offsets().shape[0]
        #print self.polygons.get_offsets()

        #assert(0)
        text_labels = [None] * len(collision)
        for i in range(len(collision)):
            temp = np.array(positions[i])
            x = np.mean(temp[:,0])
            y = np.mean(temp[:,1])
            text_labels[i]= plt.text(x,y,str(i+num_circles),color = 'k')

        self.polygon_labels = text_labels

        plt.draw()
        plt.pause(pause_time)

    #display the total area covered
    def show_title(self,area,pause_time=0.1):
        title = plt.text(-self.width/4,self.height,'total covered area =  '+ str(area),color = 'k',fontsize = 20)
        self.title.remove()
        self.title = title


#testing 
if __name__ == '__main__':
    print 'graphics testing here....\n'

    #initialize animator
    animator = Animator(20,20)

    #initialize a few circles and polygon objects
    circ1 = obj_circle(np.array([0.0,0.0]),1.0)
    circ2 = obj_circle(np.array([1.5,0]),1.0)
    circ3 = obj_circle(np.array([8.0,8.0]),2.0)

    poly1 = obj_polygon([(0,0),(3,0),(3,2),(0,2),(0,0)],np.array([2.0,0.0]))
    poly2 = obj_polygon([(0,0),(3,6),(0,4),(0,0)],np.array([3.0,0.0]))
    poly3 = obj_polygon([(0,0),(2,0),(1,1),(0,0)],np.array([-6.0,0.0]))
   
    #form arrays of objects
    circ_list= [circ1,circ2,circ3]
    poly_list= [poly1,poly2,poly3]

    #formulate a list and do all pre-processing (finding distance + collision, etc)
    template = np.array([-10.0,10.0,-10.0,10.0])
    item_lists = object_lists(circ_list,poly_list,template)

    #add to animator
    animator.add_circular_objects(item_lists.circ_diameter,\
                                  item_lists.circ_position,\
                                  item_lists.circ_collision)

    animator.add_polygon_objects(item_lists.poly_verts,\
                                  item_lists.poly_collision)


    #time.sleep(100)
    #number of optimization iterations
    num_items_renew = 100
    num_position_update = 1000
    
    num_obj = item_lists.num_circles+ item_lists.num_polygons
    while True:
        dx = 2* (np.random.random_sample((num_obj,2))- 0.5)
        #dx = np.array([[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0],[0.0,0.0]])
        #print dx
        item_lists.update_delta_pos(dx)
        #item_lists.print_items_info()

        animator.update_circular_objects(item_lists.circ_position,item_lists.circ_collision)
        #print 'next vertices positions are:', item_lists.poly_verts
        animator.update_polygon_objects(item_lists.poly_verts,item_lists.poly_collision)
    '''
        #switch in and out a new subset of items
        for i in range(num_items_renew):

            #reposition the items
            for j in range(num_position_update):
    '''
