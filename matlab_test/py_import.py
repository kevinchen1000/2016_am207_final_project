#test script for importing matlab file
import scipy.io as sio
import numpy as np
shape_contents = sio.loadmat('tiling_test.mat')
print shape_contents['wings_location'][0][0][0]


offsets = np.array([[0,1],[1,0],[0,0]])
indices = np.array([1,0,0,1,0,0,1])
#sio.savemat('tiled_result.mat',{'offsets':offsets})
sio.savemat('tiled_result.mat',{'offsets':offsets,'indices':indices})
