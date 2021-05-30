# Circles Detector
## **Circle Detection using the Circle Hough Transform and Phase Coding**

The circle Hough Transform is a technique used in image processing for detecting circles in imperfect images. The circle candidates are produced by “voting” in the Hough parameter space and then selecting local maxima in an accumulator array.

The problem with the classical Hough Transform is the requirement of a 3-D array for storing votes for multiple radii (2-D for the coordinates of the center and the 3rd Dimension for the radii), which results in large storage requirements and long processing times. 

Phase Coding method solve this problem by using a single 2-D accumulator array for all the radii. The key idea is the use of complex values in the accumulator array with the radius information encoded in the phase of the array entries.


## A Glimpse Into the Algorithm

imperfect image

![imperfect image](https://user-images.githubusercontent.com/82455000/120093199-4f183480-c121-11eb-9344-a760d5e5d0e2.png)

absolute accumulator array

![1](https://user-images.githubusercontent.com/82455000/120093519-81c32c80-c123-11eb-91a0-38c90328bf7b.png)
![2](https://user-images.githubusercontent.com/82455000/120093520-82f45980-c123-11eb-966f-ae37b451c69f.png)
![3](https://user-images.githubusercontent.com/82455000/120093521-84258680-c123-11eb-8f1c-261939728352.png)

detected circles

![found](https://user-images.githubusercontent.com/82455000/120093653-5c82ee00-c124-11eb-8eba-6d6339b3e88f.png)
