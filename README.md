# Circle Hough Transform 
## Circle Detection using the Circle Hough Transform and Phase Coding

The Circle Hough Transform is a technique used in image processing for detecting imperfect circles in  images. The circle candidates are produced by “voting” in the Hough parameter space and then selecting local maxima in an accumulator array.

The problem with the classical Hough Transform is the requirement of a 3-D array for storing votes for multiple radii (2-D for the coordinates of the center and the 3rd Dimension for the radii), which results in large storage requirements and long processing times. 

Phase Coding method solve this problem by using a single 2-D accumulator array for all the radii. The key idea is the use of complex values in the accumulator array with the radii information encoded in the phase of the array entries.

But, a new problem arises using the Phase Coding method; that is the accuracy of the radius estimated by decoding the phase information from the estimated center location in the accumulator array.

To solve this problem I use Phase Coding twice for each circle, once to get a global accumulator array and than one more time to form a local accumulator array of a thin ring of radii around the estimated center and radius I've found out of the global accumulator. from the local accumulator array I get a more accurate center candidate, than I use a kind of voting process to get even more accurate center and radius, avoiding the complexity of the classical algorithm by considering only the near vicinity pixels of the estimated center as candidates for the center and than calculating the distances between each of the center candidates and the pixels found in the thin ring. More over, I keep the Phase Coding approach and use a complex 2-D array instead of 3-D array, storing the estimated radius in the real part and its votes in the imaginary part of each array entries.
## A Glimpse Into the Algorithm

imperfect circles (left) vs detected circles (right)

![dotted and found](https://user-images.githubusercontent.com/82455000/120094977-fac68200-c12b-11eb-97f7-41884f14ddfe.png)

global and local absolute accumulator arrays through the process

![first](https://user-images.githubusercontent.com/82455000/120094854-78d65900-c12b-11eb-8180-c933029b5c56.png)
![second](https://user-images.githubusercontent.com/82455000/120094865-84c21b00-c12b-11eb-817f-fc00bfdba170.png)
![third](https://user-images.githubusercontent.com/82455000/120094878-9277a080-c12b-11eb-96e3-1107c55c06ea.png)
## Installation and Usage
Save the scripts in the same folder.

Run DottedCirclesMaker.m to create and save an image of impefect circles. In addition to the imperfect circles image, the script plots a figure of the original whole circles which the imperfect circles created from. this can be used for comparison with the detcted circles.<br/>
as the fnction input arguments you can insert the number of dotted circles you wish to create, the number of noise dots and the minimum number of dots for each dotted circle. also, you can call the funtion without any arguments to get a random number of dotted circles betwwen 3 to 10, without noise dots at all and a minimum number of 30 dots for each dotted circle. you can also insert only one argument to choose the number of circles or only two arguments to choose the number of noise points too.

After creating the imperfect circles image, run CirclsDetector.m to detect the circles. the script will plot the detected circles and the rest of the unused pixels. in the function call you can insert the minimun number of dots needed to form a circle and also the minimum and maximun radius to set the radii range of the circles you wish to find, or you can leave it empty to set the defult minimum of 3 dots and the radii range between quarter and half of the image side. you can also insert only the first argument, leaving the radii range to the default set.

Enjoy!

## References
[T.J Atherton, D.J. Kerbyson. "Size invariant circle detection." Image and Vision Computing. Volume 17, Number 11, 1999, pp. 795-803.](https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.24.3310&rep=rep1&type=pdf)
