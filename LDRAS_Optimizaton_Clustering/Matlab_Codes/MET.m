function [ Level ] = MET( IMG )
%Maximum Error Thresholding By Kittler
%   Finding the Min of a cost function J in any possible thresholding. The
%   function output is the Optimal Thresholding.

for t = 0:255 % Assuming 8 bit image
    I1 = IMG;
    I1 = I1(I1 <= t);
    q1 = sum(hist(I1, 256));

    I2 = IMG;
    I2 = I2(I2 > t);
    q2 = sum(hist(I2, 256));

    % J is proportional to the Overlapping Area of the 2 assumed Gaussians
    J(t + 1) = 1 + 2 * (q1 * log(std(I1, 1)) + q2 * log(std(I2, 1)))...
        -2 * (q1 * log(q1) + q2 * log(q2));
end

[~, Level] = min(J);

%Level = (IMG <= Level);

end