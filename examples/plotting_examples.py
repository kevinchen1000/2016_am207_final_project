#example script to plot a moving colored circle

import matplotlib.pyplot as plt
from matplotlib.collections import EllipseCollection
import numpy as np
from timer import Timer
import time

def plot_template(l = 20.0, w = 20.0):
  x_pos = [-l/2, l/2, l/2, -l/2, -l/2 ]
  y_pos = [-w/2, -w/2, w/2, w/2, -w/2]
  plt.plot(x_pos,y_pos,'k-',linewidth = 3.0)
  plt.xlabel('x')
  plt.ylabel('y')
  #plt.show()

def randcolor():
    return np.random.uniform(0.0, 0.89, (3,)) + 0.1

class Animator(object):
    def __init__(self, positions, diameter):
        self.count = positions.shape[0]

        plt.ion()
        fig = plt.figure(figsize=(10, 10))
        ax = fig.gca()
        self.ax = ax

        diameters = np.ones(self.count) * diameter
        colors = [randcolor() for _ in range(self.count)]
        self.circles = EllipseCollection(widths=diameters,
                                         heights=diameters,
                                         angles=np.zeros_like(diameters),
                                         units='xy',
                                         offsets=positions,
                                         transOffset=ax.transData,
                                         edgecolor='face', facecolor=colors)
        ax.add_collection(self.circles)

        ax.axis([0, 1, 0, 1])
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.set_axis_bgcolor('black')
        plt.draw()
        plt.show()

    def update(self, positions):
        self.circles.set_offsets(positions)
        plt.draw()
        plt.show()

'''
if __name__ == '__main__':
    print 'start'
    #plot_template()
    circle1=plt.Circle((0,0),2,color='r')
    fig = plt.gcf()
    fig.gca().add_artist(circle1)
    plt.draw()
    plt.show()
'''

if __name__ == '__main__':
    num_balls = 10000
    radius = 0.004
    positions = np.random.uniform(0 + radius, 1 - radius,
                                  (num_balls, 2)).astype(np.float32)
    velocities = np.random.uniform(-0.5, 0.5,
                                   (num_balls, 2)).astype(np.float32)

    animator = Animator(positions, radius * 2)


    anim_step = 1.0 / 30  # FPS
    frame_count = 0
    total_time = 0

    while True:
        with Timer() as t:

            positions = np.random.uniform(0 + radius, 1 - radius,
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


