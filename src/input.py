import numpy as np
import objects

def input_generator(animator):
	'''
	Parameters
    ----------
    animator : class Animator
        the animator object we are working on


    Returns
    -------
    class object_lists
    	contains several circles and polygons in random positions.
    	Total sum of all obejcts area is greater than Animator canvas.
	'''
	input_area = 0.0
	animator_area = animator.width * animator.height
	circle_list = []
	polygon_list = []
	
	while input_area < animator_area:
		# a random [0,1] number whether to make a new circle or a new polygon
		is_circle = np.random.randint(low=0, high=2)

		if is_circle:
			# randomly pick radius size
			raidus = np.random.uniform(low=0.1, high=animator.width*1./4.0)
			center = np.random.uniform(low=-animator.width*1./2.+radius, high=animator.width*1./2.-radius)
			circle = objects.obj_circle(np.array([center, center], radius))
			circle_list.append(circle)
			input_area += circle.area
		else:
			# randomly pick either triangle, square
			num_vertex = np.random.randint(low=3, high=5)
			vertex = [(0,0)]
			
			if num_vertex == 3: #triangle
				# point 1
				bottom_point = np.random.uniform(low=0.1, high=animator.width*1./2.0)
				vertex.append((0, bottom_point))
				# point 2
				top_point = np.random.uniform(low=0.1, high=bottom_point)
				vertex.append((top_point, top_point))
				# add (0,0)
				vertex.append((0,0))
				# offset
				offset_x = np.random.uniform(low=(-animator.width+bottom_point)*1./2., high=(animator.width-bottom_point)*1./2.)
				offset_y = np.random.uniform(low=(-animator.width+top_point)*1./2., high=(animator.width-top_point)*1./2.)

				triangle = objects.obj_polygon(vertex, np.array([offset_x, offset_y]))
				polygon_list.append(triangle)

				input_area += triangle.area

			elif num_vertex == 4: #rectangle
				# point 1
				bottom_point = np.random.uniform(low=0.1, high=animator.width*1./2.0)
				vertex.append((0, bottom_point))
				# point 2
				top_point = np.random.uniform(low=0.1, high=bottom_point)
				vertex.append((0, top_point))
				# point 3
				vertex.append((bottom_point, top_point))
				# add (0,0)
				vertex.append((0,0))
				# offset
				offset_x = np.random.uniform(low=(-animator.width+bottom_point)*1./2., high=(animator.width-bottom_point)*1./2.)
				offset_y = np.random.uniform(low=(-animator.width+top_point)*1./2., high=(animator.width-top_point)*1./2.)

				rectangle = objects.obj_polygon(vertex, np.array([offset_x, offset_y]))
				polygon_list.append(rectangle)

				input_area += rectangle.area

	# when input area is greater than Animator canvas area
	input_object = objects.object_lists(circle_list, polygon_list)
	return input_object

			
