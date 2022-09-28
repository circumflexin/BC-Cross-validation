#%% shebangs
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 17 10:35:14 2021
@author: Alec Burslem; acb35@st-andrews.ac.uk; github.com/circumflexin
"""

#to do
"""
make sure to use
3x3 or 4x4 ArUco library as the marker squares are larger.
"""

#%% load packages & set directories.
import cv2
import numpy as np
import os
import math
import pickle
in_dir = "D:\\Analysis\\BC Cross validation\\UAV photogrammetry\\calibration\\scaletest"
out_dir = "D:\\Analysis\\BC Cross validation\\UAV photogrammetry\\calibration"
os.chdir(in_dir)

#%% help the user design their calibration session
#work out marker size for printing
#minimum marker resolution = 20*20 px
#optimum marker resolution = 30*30 px

height = height = 10*10**3 #(mm),

#p4
sensor_width = 6.3 #mm
focal_length = 3.61 #mm
imwidth = 3840 #mm
fov = 94 #degrees

# ##p4pro
# sensor_width = 13.2 #mm
# focal_length = 8.8 #mm
# imwidth = 5472 #mm
# fov = 84 #degrees



approx_gsd = ((height*sensor_width)/(focal_length*imwidth))/10 #cm per pixel
print(approx_gsd) # cm / px

targ_sz = round(30*approx_gsd) #target size
fov = math.radians(fov)
samp_area =  (math.tan(fov/2) * height)*2*10**-3 # ground area imaged by
# the drone at the height specified (m)
grid_area = samp_area*0.6 # its good to cover at least 60% of the image.

#%% set up grid parameters
aruco_dict = cv2.aruco.getPredefinedDictionary(cv2.aruco.DICT_4X4_50)
board = cv2.aruco.GridBoard_create(6,7,0.24,0.001,aruco_dict)
# 	squaresX = 6,
# 	squaresY = 7,
# 	squareLength = .025,
# 	sep = .001,
# 	dictionary = aruco_dict
img = cv2.aruco.drawMarker(aruco_dict,1, 700)
imboard = board.draw((2000,2000))
#cv2.imshow('aruco.png',imboard)
#cv2.waitKey(0)

# create markers for printing
dpi = 300
dpcm = dpi/2.54
sidePixels = round(dpcm * targ_sz)
for i in range(6*7):
	marker = cv2.aruco.drawMarker(aruco_dict,i,sidePixels)
	name = ''.join(['marker',str(i),'.png'])
	cv2.imwrite(os.path.join(in_dir,'markers',name),marker)
cv2.imwrite(os.path.join(in_dir,'markers','aruco_board.png'),imboard)

#%%  load the input video, show the marker detection and perform the calibration

all_counts = []
cap = cv2.VideoCapture('test_cal_Trim.mp4')
n = 1
while cap.isOpened():
	ret, frame = cap.read()
	if ret == True:
		grey = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
		parameters = cv2.aruco.DetectorParameters_create()
		parameters.cornerRefinementMethod = cv2.aruco.CORNER_REFINE_SUBPIX
		#detect markers
		corners,ids, rejectedImgPoints = cv2.aruco.detectMarkers(grey, aruco_dict, parameters=parameters)
		# make sure there are markers
		if ids is not None and len(ids) > 11:
			if n == 1:
				all_corners = corners
				all_ids = ids
			else:
				all_ids = np.concatenate([all_ids,ids])
				all_corners = np.concatenate([all_corners,corners])
			all_counts.append(len(ids))
			print(n)
			n+=1
			if cv2.waitKey(1) & 0xFF == ord('q'):
				break
	else:
		break
cap.release()
cv2.destroyAllWindows()

all_counts = np.array(all_counts)

ret2,cameraMatrix, distCoeffs,_,_ = cv2.aruco.calibrateCameraAruco(
 				corners = all_corners,
 				ids = all_ids,
 				counter = all_counts,
 				board = board,
 				imageSize = grey.shape,
 				cameraMatrix=None,
 				distCoeffs = None
 				)

print(cameraMatrix)
print(distCoeffs)
log = open('cal_log.txt','wb')
pickle.dump((cameraMatrix, distCoeffs), log)
log.close()
print('calibration complete')

#%% perform and record calibration diagnostics
cap = cv2.VideoCapture('scale_test2.mp4')
# Define the codec and create VideoWriter object
fourcc = cv2.VideoWriter_fourcc(*'MJPG')
out = cv2.VideoWriter('scale_test_out.avi',fourcc,20, (1920,1080),True)
#for each frame in the video obeject, apply the calibration object and project the orientation.
while(cap.isOpened()):
	ret, frame = cap.read()
	if ret == True:
		corners, ids, rejectedImgPoints = cv2.aruco.detectMarkers(frame, aruco_dict, parameters=parameters)
		if ids is not None:
			for i in range(0, len(ids)):
				rvec, tvec, markerPoints = cv2.aruco.estimatePoseSingleMarkers(corners[i], 0.02, cameraMatrix,distCoeffs)
				(rvec - tvec).any()
				cv2.aruco.drawDetectedMarkers(frame,corners,ids)
				cv2.aruco.drawAxis(frame, cameraMatrix, distCoeffs, rvec, tvec, 0.03)
			cv2.imshow('frame', frame)
		out.write(frame)
	else:
		break
	if cv2.waitKey(1) & 0xFF == ord('q'):
		break
cap.release()
out.release()
cv2.destroyAllWindows()

