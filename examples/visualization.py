#script to visualize the template and filled with circles or polygons

import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
from matplotlib.collections import PolyCollection
import numpy as np
from timer import Timer
import time
import matplotlib

class Animator(object):
    def __init__(self, positions, diameter):
        self.circ_count = positions.shape[0]

        plt.ion()
        fig = plt.figure(figsize=(10, 10))
        ax = fig.gca()
        self.ax = ax

        diameters = np.ones(self.circ_count) * diameter
        circ_colors = [(0.0,0.0,1.0) for _ in range(self.circ_count)]

        #add circles
        self.circles = EllipseCollection(widths=diameters,
                                         heights=diameters,
                                         angles=np.zeros_like(diameters),
                                         units='xy',
                                         offsets=positions,
                                         transOffset=ax.transData,
                                         edgecolor='face', facecolor=circ_colors)
        ax.add_collection(self.circles)

        #add polygons
        self.poly_count = 3
        verts =[[(0,1),(1,0),(2,2)],[(6,5),(3,7),(7,6)],[(0,-1),(-1,0),(4,5)]]
        poly_colors = [(0.0,0.0,1.0) for _ in range(self.poly_count)]
        self.polygons = PolyCollection(verts, facecolors=poly_colors)

        ax.add_collection(self.polygons)

        ax.axis([-20, 20, -20, 20])
        ax.get_xaxis().set_visible(True)
        ax.get_yaxis().set_visible(True)
        #ax.set_axis_bgcolor('black')
        plt.draw()
        plt.pause(0.1)

    def update(self, positions):
        
        #update number of balls

        self.circles.set_offsets(positions)
        colors = [(1.0,0.0,0.0) for _ in range(self.circ_count)]
        self.circles.set_facecolors(colors)
        
        #redefine circles 
        diameter= 2
        diameters = np.ones(self.circ_count+100) * diameter
        circ_colors = [(0.0,0.0,1.0) for _ in range(self.circ_count)]

        #add circles
        self.circles.remove()
        self.circles = EllipseCollection(widths=diameters,
                                         heights=diameters,
                                         angles=np.zeros_like(diameters),
                                         units='xy',
                                         offsets=positions,
                                         transOffset=self.ax.transData,
                                         edgecolor='face', facecolor=circ_colors)
        self.ax.add_collection(self.circles)

        #label
        text_labels = [None] * positions.shape[0]
        for i in range(positions.shape[0]):
            text_labels[i] =  plt.text(positions[i,0],positions[i,1],str(i),color='white')
            

        
        

        #remove polygon
        #self.polygons.remove()
        poly_color = np.random.uniform(0.0, 0.89, (3,)) + 0.1
        verts =[[(10,1),(9,0),(8,2),(10,2)],[(8,8),(9,9),(6,7)]]
        self.polygons.set_verts(verts)
        
        plt.draw()
        plt.pause(0.1)

        #remove labels
        for i in range(positions.shape[0]):
            text_labels[i].remove()
 

if __name__ == '__main__':
    num_balls = 10
    radius = 0.4
    positions = np.random.uniform(-20 + radius, 20 - radius,
                                  (num_balls, 2)).astype(np.float32)
    velocities = np.random.uniform(-0.5, 0.5,
                                   (num_balls, 2)).astype(np.float32)

    animator = Animator(positions, radius * 2)


    anim_step = 1.0 / 30  # FPS
    frame_count = 0
    total_time = 0

    while True:
        with Timer() as t:

            positions = np.random.uniform(-20 + radius, 20 - radius,
                                  (num_balls, 2)).astype(np.float32)
            
        total_time += t.interval
        frame_count += 1
        if total_time > anim_step:
            animator.update(positions)
            #plt.show()
            #time.sleep(1)
            print("{} simulation frames per second".format(frame_count / total_time))
            frame_count = 0
            total_time = 0


            #update(positions, velocities, grid,
            #       radius, grid_size, physics_step, locks_ptr)
