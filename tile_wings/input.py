import numpy as np
from objects import *
from animation import *
import scipy.io as sio

def input_generator(template):
	'''
	Parameters
	----------
	template : numpy array for representing the template
				we set it to be [-10, 10, -10, 10] so total area = 20*20

	Returns
	-------
	class object_lists
		contains several circles and polygons in random positions.
		Total sum of all obejcts area is greater than Animator canvas.
	'''

	input_area = 0.0
	animator_width = template[1]-template[0] # 20
	animator_height = template[3]-template[2] # 20
	animator_area = animator_width * animator_height
	circle_list = []
	polygon_list = []
	
	while input_area < animator_area:
		# a random [0,1] number whether to make a new circle or a new polygon
		is_circle = np.random.randint(low=0, high=2)

		if is_circle:
			# randomly pick radius size
			radius = np.random.uniform(low=0.1, high=animator_width*1./4.)
			# maybe later when we want to put object starting from random positions...
			center = np.random.uniform(low=-animator_width*1./2.+radius, high=animator_width*1./2.-radius)
			circle = obj_circle(np.array([0., 0.]), radius)
			circle = obj_circle(np.array([0.0,0.0]), radius)
			circle_list.append(circle)
			input_area += circle.area
		else:
			# randomly pick either triangle(3), square(4), pentagon(5), hexagon(6)
			num_vertex = np.random.randint(low=3, high=7)
			vertex = [(0,0)]

			if num_vertex == 3: #triangle = [(0,0), (x1, y1), (x2, y2), (0,0)]
				# (x1, y1)
				x1 = np.random.uniform(low=-6., high=6.)
				y1 = np.random.uniform(low=-6., high=6.)
				# bottom_point = np.random.uniform(low=0.1, high=animator_width*1./2.0)
				vertex.append((x1, y1))
				# (x2, y2)
				# top_point = np.random.uniform(low=0.1, high=bottom_point)
				x2 = np.random.uniform(low=-6., high=6.)
				y2 = np.random.uniform(low=-6., high=6.)
				vertex.append((x2, y2))
				# add (0,0)
				vertex.append((0,0))
				# offset
				# offset_x = np.random.uniform(low=(-animator.width+bottom_point)*1./2., high=(animator.width-bottom_point)*1./2.)
				# offset_y = np.random.uniform(low=(-animator.width+top_point)*1./2., high=(animator.width-top_point)*1./2.)

				triangle = obj_polygon(vertex, np.array([0.0,0.0]))

				polygon_list.append(triangle)

				input_area += triangle.area

			elif num_vertex == 4: #rectangle
				# [(0,0),(3,0),(3,2),(0,2),(0,0)]
				# [(0,0),(4,0),(4,2),(0,2),(0,0)]
				# [(0,0),(1,0),(1,6),(0,6),(0,0)]
				# [(0,0),(3,0),(5,2),(4,3),(0,0)]

				# point 1
				x1 = np.random.uniform(low=0.1, high=6.)
				# bottom_point = np.random.uniform(low=0.1, high=animator.width*1./2.0)
				vertex.append((x1, 0))
				# point 2
				# top_point = np.random.uniform(low=0.1, high=bottom_point)
				x2 = x1 + np.random.uniform(low=0, high=x1)
				y2 = np.random.uniform(low=0.1, high=6.)
				vertex.append((x2, y2))
				# point 3
				x3 = x2 - np.random.uniform(low=1., high=6.)
				y3 = y2 + np.random.uniform(low=-y2+1., high=2.)
				vertex.append((x3,y3))
				# add (0,0)
				vertex.append((0,0))
				# offset
				# offset_x = np.random.uniform(low=(-animator.width+bottom_point)*1./2., high=(animator.width-bottom_point)*1./2.)
				# offset_y = np.random.uniform(low=(-animator.width+top_point)*1./2., high=(animator.width-top_point)*1./2.)

				rectangle = obj_polygon(vertex, np.array([0., 0.]))
				polygon_list.append(rectangle)

				input_area += rectangle.area

			elif num_vertex == 5: #pentagon
				# [(0,0),(2,0),(3,1),(1,3),(-1,1),(0,0)]
				# [(0,0),(2,0),(3,1),(1,3),(-1,1),(0,0)]

				# point 1
				x1 = np.random.uniform(low=0.1, high=3.)
				vertex.append((x1, 0))

				# point 2
				x2 = x1 + np.random.uniform(low=0, high=3.)
				y2 = np.random.uniform(low=1., high=5.)
				vertex.append((x2, y2))

				# point 3
				x3 = x2 - np.random.uniform(low=1., high=3.)
				y3 = y2 + np.random.uniform(low=1., high=3.)
				vertex.append((x3,y3))

				# point 4
				x4 = x3 - np.random.uniform(low=1., high=3.)
				y4 = y3 - np.random.uniform(low=1., high=y3-1.)
				vertex.append((x4, y4))
				# add (0,0)
				vertex.append((0,0))
				# offset
				# offset_x = np.random.uniform(low=(-animator.width+bottom_point)*1./2., high=(animator.width-bottom_point)*1./2.)
				# offset_y = np.random.uniform(low=(-animator.width+top_point)*1./2., high=(animator.width-top_point)*1./2.)

				pentagon = obj_polygon(vertex, np.array([0., 0.]))
				polygon_list.append(pentagon)
				input_area += pentagon.area

			elif num_vertex == 6: #hexagone
				# [(0,0),(2,0),(3,1.5),(2,3),(0,3),(-1,1.5),(0,0)]
				# [(0,0),(2,0),(3,1.5),(2,3),(0,3),(-1,1.5),(0,0)]

				# point 1
				x1 = np.random.uniform(low=0.1, high=3.)
				vertex.append((x1, 0))

				# point 2
				x2 = x1 + np.random.uniform(low=0, high=3.)
				y2 = np.random.uniform(low=1., high=5.)
				vertex.append((x2, y2))

				# point 3
				x3 = x2 - np.random.uniform(low=1., high=3.)
				y3 = y2 + np.random.uniform(low=1., high=3.)
				vertex.append((x3,y3))

				# point 4
				x4 = x3 - np.random.uniform(low=1., high=3.)
				y4 = y3 - np.random.uniform(low=1., high=y3-1.)
				vertex.append((x4, y4))

				# point 5
				x5 = x4 - np.random.uniform(low=1., high=3.)
				y5 = y4 - np.random.uniform(low=1., high=y4-1.)
				# add (0,0)
				vertex.append((0,0))
				# offset
				# offset_x = np.random.uniform(low=(-animator.width+bottom_point)*1./2., high=(animator.width-bottom_point)*1./2.)
				# offset_y = np.random.uniform(low=(-animator.width+top_point)*1./2., high=(animator.width-top_point)*1./2.)

				hexagon = obj_polygon(vertex, np.array([0., 0.]))
				polygon_list.append(hexagon)
				input_area += hexagon.area

	# when input area is greater than Animator canvas area
	# input_object = object_lists(circle_list, polygon_list, template)
	return circle_list, polygon_list




if __name__ == '__main__': 
    template = np.array([-10.0,10.0,-10.0,10.0])
    print input_generator(template)

    circ_list, poly_list = input_generator(template)

    item_lists = object_lists(circ_list,poly_list,template)

    #add to animator
    animator = Animator(20,20)
    animator.add_circular_objects(item_lists.circ_diameter,\
                                  item_lists.circ_position,\
                                  item_lists.circ_collision,0.1)

    animator.add_polygon_objects(item_lists.poly_verts,\
                                  item_lists.poly_collision,100)

