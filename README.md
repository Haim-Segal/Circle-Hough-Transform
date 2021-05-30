# Circles Detector
## **Circle Detection using the Circle Hough Transform and Phase Coding**

The circle Hough Transform is a technique used in image processing for detecting circles in imperfect images. The circle candidates are produced by “voting” in the Hough parameter space and then selecting local maxima in an accumulator array.

The problem with the classical Hough Transform is the requirement of a 3-D array for storing votes for multiple radii (2-D for the coordinates of the center and the 3rd Dimension for the radii), which results in large storage requirements and long processing times. 

Phase Coding method solve this problem by using a single 2-D accumulator array for all the radii. The key idea is the use of complex values in the accumulator array with the radius information encoded in the phase of the array entries.


## A Glimpse Into the Algorithm

imperfect image vs detected circles

![dotted and found](https://user-images.githubusercontent.com/82455000/120094977-fac68200-c12b-11eb-97f7-41884f14ddfe.png)

absolute accumulator array through the process

![first](https://user-images.githubusercontent.com/82455000/120094854-78d65900-c12b-11eb-8180-c933029b5c56.png)
![second](https://user-images.githubusercontent.com/82455000/120094865-84c21b00-c12b-11eb-817f-fc00bfdba170.png)
![third](https://user-images.githubusercontent.com/82455000/120094878-9277a080-c12b-11eb-96e3-1107c55c06ea.png)


