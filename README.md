# Circles Detector
**Circle Detection using the Circle Hough Transform and Phase Coding**

The circle Hough Transform is a technique used in image processing for detecting circles in imperfect images. The circle candidates are produced by “voting” in the Hough parameter space and then selecting local maxima in an accumulator array.

The problem with the classical Hough Transform is the requirement of a 3-D array for storing votes for multiple radii (2-D for the coordinates of the center and the 3rd Dimension for the radii), which results in large storage requirements and long processing times. 

Phase Coding method solve this problem by using a single 2-D accumulator array for all the radii. The key idea is the use of complex values in the accumulator array with the radius information encoded in the phase of the array entries.


# A Glimpse Into the Algorithm
![מסמך ללא שם-1](https://user-images.githubusercontent.com/82455000/120086676-a26e9080-c0e9-11eb-952b-878eab33c9fc.png)
