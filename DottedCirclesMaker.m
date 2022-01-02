clear, clc, close all
tic
dottedCircles()
toc



function [radii, xCen, yCen, pointsSize] = setRadiiAndCenters( ...
        numOfCircles, scaleConst)
    radii = randi(round(scaleConst * [0.25, 0.5]), numOfCircles);
    xCen = randi((scaleConst - 20) * [- 1, 1], numOfCircles);
    yCen = randi((scaleConst - 20) * [- 1, 1], numOfCircles);
    pointsSize = rand(numOfCircles, 1) * 2 + 1;
end

function [xNoise, yNoise, noiseSize] = setNoise( ...
        numOfNoisePoints, scaleConst)
    xNoise = randi((scaleConst - 2) * [- 1, 1], 1, numOfNoisePoints);
    yNoise = randi((scaleConst - 2) * [- 1, 1], 1, numOfNoisePoints);
    noiseSize = rand * 2 + 1;
end

function plotCorners(scaleConst)
    figure('visible', 'off')
    hold on
    plot(scaleConst * [- 1, 1, 1, - 1], scaleConst * [- 1, - 1, 1, 1], ...
        '.k', 'MarkerSize', 0.1)
    axis equal off
end

function outsideBordersPoints = rotateOutsideBordersPoints( ...
        rotMat, x, y, xCen, yCen, outsideBordersPoints, circle, SCALE_CONST)
    rotatedPoints = rotMat * [x(outsideBordersPoints) - xCen(circle);
                                y(outsideBordersPoints) - yCen(circle)];
    x(outsideBordersPoints) = rotatedPoints(1, :) + xCen(circle);
    y(outsideBordersPoints) = rotatedPoints(2, :) + yCen(circle);
    outsideBordersPoints = find( ...
        abs(x) > SCALE_CONST | abs(y) > SCALE_CONST);
end

function [x, y] = KeepInsideBordersPointsOnly(wholeCircle, SCALE_CONST, minPoints, ...
        x, y, xCen, yCen, circle)
    outsideBordersPoints = find(abs(x) > SCALE_CONST | abs(y) > SCALE_CONST);
    if wholeCircle
        while length(x) - length(outsideBordersPoints) < minPoints
            outsideBordersPoints = rotateOutsideBordersPoints( ...
                [0, - 1; 1, 0], x, y,xCen, yCen, outsideBordersPoints, ...
                circle, SCALE_CONST);
        end
    end
    x(outsideBordersPoints) = [];
    y(outsideBordersPoints) = [];
end

function plotAndSaveCircles(wholeCircle ,imgName, pointsSize, radii, ...
        xCen, yCen, minPoints, numOfCircles, SCALE_CONST, xNoise, yNoise, noisePointsSize)
    plotCorners(SCALE_CONST)
    for circle = 1:numOfCircles
        if wholeCircle
            theta = linspace(1 / radii(circle), 2 * pi, 2 * pi * radii(circle));
        else
            theta = rand(1, randi([minPoints minPoints + 10])) * 2 * pi;
        end
        x = radii(circle) .* cos(theta) + xCen(circle);
        y = radii(circle) .* sin(theta) + yCen(circle);
        [x, y] = KeepInsideBordersPointsOnly(wholeCircle, SCALE_CONST, ...
            minPoints, x, y, xCen, yCen, circle);
        plot(x, y, '.k', "MarkerSize", pointsSize(circle))
    end
    plot(xNoise, yNoise, '.k', "MarkerSize", noisePointsSize)
    saveas(gcf, imgName)
end

function mat = removeWhiteBorders(mat)
    for side = 1:4
        while all(mat(:, 1))
            mat(:, 1) = [];
        end
        mat = rot90(mat);
    end
    mat([1:2, end - 1:end], [1:2, end - 1:end]) = 1;
end

function imageShow(imageName)
    figure
    image = imbinarize(im2gray(imread([imageName,'.png'])));
    image = removeWhiteBorders(image);
    imshow(image)
end

function dottedCircles(numOfCircles, numOfNoisePoints, minPoints, SCALE_CONST)

    function setInputArg(numOfInputArguments)
        if numOfInputArguments < 4
            SCALE_CONST = 200;
            if numOfInputArguments < 3
                minPoints = 30;
                if numOfInputArguments < 2
                    numOfNoisePoints = 0;
                    if numOfInputArguments < 1
                        numOfCircles = randi([3 10]);
                    end
                end
            end
        end
    end
        
    setInputArg(nargin);
    [radii, xCen, yCen, pointsSize] = setRadiiAndCenters( ...
        numOfCircles, SCALE_CONST);
    [xNoise, yNoise, noisePointsSize] = setNoise( ...
        numOfNoisePoints, SCALE_CONST);
    plotAndSaveCircles(false, 'DottedCircles.png', pointsSize, radii, ...
        xCen, yCen, minPoints, numOfCircles, SCALE_CONST, xNoise, ...
        yNoise, noisePointsSize)
    plotAndSaveCircles(true, 'WholeCircles.png', pointsSize, radii, ...
        xCen, yCen, minPoints, numOfCircles, SCALE_CONST, xNoise, ...
        yNoise, noisePointsSize)
    close all
    imageShow('DottedCircles')
    imageShow('WholeCircles')
end