clear, clc, close all
tic
dottedCircles(9)
toc



function [radii, xCen, yCen, pointsSize] = setRadiiAndCenters( ...
        numOfcircles, scaleConst)
    radii = randi(round(scaleConst * [0.25, 0.5]), numOfcircles, 1);
    xCen = randi((scaleConst - 20) * [- 1, 1], numOfcircles, 1);
    yCen = randi((scaleConst - 20) * [- 1, 1], numOfcircles, 1);
    pointsSize = rand(numOfcircles, 1) * 2 + 1;
end

function [xNoise, yNoise, markerSizeNoise] = setNoise(numOfNoisePoints, scaleConst)
    xNoise = randi((scaleConst - 2) * [- 1, 1], 1, numOfNoisePoints);
    yNoise = randi((scaleConst - 2) * [- 1, 1], 1, numOfNoisePoints);
    markerSizeNoise = rand * 2 + 1;
end

function plotCorners(scaleConst)
    figure('visible', 'off')
    hold on
    plot(scaleConst * [- 1, 1, 1, - 1], scaleConst * [- 1, - 1, 1, 1], ...
        '.k', 'MarkerSize', 0.1)
    axis equal off
end

function outBordersPoints = rotateOutBordersPoints(rotMat, x, y, ...
        xCen, yCen, outBordersPoints, circle, SCALE_CONST)
    z = rotMat * [x(outBordersPoints) - xCen(circle); y(outBordersPoints) - yCen(circle)];
    x(outBordersPoints) = z(1, :) + xCen(circle);
    y(outBordersPoints) = z(2, :) + yCen(circle);
    outBordersPoints = find(abs(x) > SCALE_CONST | abs(y) > SCALE_CONST);
end

function [x, y] = KeepInBordersPointsOnly(SCALE_CONST, minPoints, x, y, ...
            xCen, yCen, circle)
    outBordersPoints = find(abs(x) > SCALE_CONST | abs(y) > SCALE_CONST);
    if nargin == 3
        while length(x) - length(outBordersPoints) < minPoints
            outBordersPoints = rotateOutBordersPoints([0, - 1; 1, 0], x, y, ...
            xCen, yCen, outBordersPoints, circle, SCALE_CONST);
        end
    end
    
    x(outBordersPoints) = [];
    y(outBordersPoints) = [];
end

function plotAndSaveDottedCircels(SCALE_CONST, minPoints, radii, xCen, yCen, ...
    pointsSize, numOfcircles, scaleConst, xNoise, yNoise, noisePointsSize)
    plotCorners(scaleConst)
    for circle = 1:numOfcircles
        t = rand(1, randi([minPoints minPoints + 10])) * 2 * pi;
        x = radii(circle) .* cos(t) + xCen(circle);
        y = radii(circle) .* sin(t) + yCen(circle);
        [x,y] = KeepInBordersPointsOnly(SCALE_CONST, minPoints, x, y, ...
            xCen, yCen, circle);
        plot(x, y, '.k', "MarkerSize", pointsSize(circle))
    end
    plot(xNoise, yNoise, '.k', "MarkerSize", noisePointsSize)
    saveas(gcf, 'DottedCircles.png')
end

function plotAndSaveWholeCircles(pointsSize, radii, xCen, yCen, minPoints, numOfcircles, SCALE_CONST, xNoise, yNoise, noisePointsSize)
    plotCorners(SCALE_CONST)
    for circle = 1:numOfcircles
        t = linspace(1 / radii(circle), 2 * pi, 2 * pi * radii(circle));
        x = radii(circle) .* cos(t) + xCen(circle);
        y = radii(circle) .* sin(t) + yCen(circle);
        [x, y] = KeepInBordersPointsOnly(SCALE_CONST, minPoints, x, y, ...
            xCen, yCen, circle);
        plot(x, y, '.k', "MarkerSize", pointsSize(circle))
    end
    plot(xNoise, yNoise, '.k', "MarkerSize", noisePointsSize)
    saveas(gcf, 'WholeCircles.png')
    close all
end

function imageShow(imageName)
    figure
    image = imbinarize(im2gray(imread([imageName,'.png'])));
    image = removeBorders(image);
    imshow(image)
end

function mat = removeBorders(mat)
    for side = 1:4
        while all(mat(:, 1))
            mat(:, 1) = [];
        end
        mat = rot90(mat);
    end
    mat([1:2, end - 1:end], [1:2, end - 1:end]) = 1;
end

function dottedCircles(numOfcircles, numOfNoisePoints, minPoints)
setInputArg(nargin);
SCALE_CONST = 200;
[radii, xCen, yCen, pointsSize] = setRadiiAndCenters(numOfcircles, SCALE_CONST);
[xNoise, yNoise, noisePointsSize] = setNoise(numOfNoisePoints, SCALE_CONST);
plotAndSaveDottedCircels(SCALE_CONST, minPoints, radii, xCen, yCen, pointsSize, numOfcircles, SCALE_CONST, xNoise, yNoise, noisePointsSize)
plotAndSaveWholeCircles(pointsSize, radii, xCen, yCen, minPoints, numOfcircles, SCALE_CONST, xNoise, yNoise, noisePointsSize)
imageShow('DottedCircles')
imageShow('WholeCircles')

    function setInputArg(numOfInputArguments)
        switch numOfInputArguments
            case 0
                numOfcircles = randi([3 10]);
                numOfNoisePoints = 0;
                minPoints = 30;
            case 1
                numOfNoisePoints = 0;
                minPoints = 30;
            case 2
                minPoints = max(30, round(3 * sqrt(numOfNoisePoints)));
        end
    end
end