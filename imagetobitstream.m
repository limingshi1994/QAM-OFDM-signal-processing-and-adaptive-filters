% Function converts imput (bmp) image to bitstream
% Input arguments
%   filename: file name of BMP image in string format
% Output arguments
%   imageBitStream: image bit stream in vector format
%   colorMap: image color map
%   imageSize: dimensions of image
%   bitsPerPixel: nb of bits per pixel

function [imageBitStream, imageData, colorMap,imageSize,bitsPerPixel] = imagetobitstream(filename)

imageObj = importdata(filename); % read in image
colorMap = imageObj.colormap;
bitsPerPixel=log2(size(colorMap,1));  % number of bits to assign color to pixel
imageData = double(imageObj.cdata);   % conversion to doubles
imageSize = size(imageData);  % dimensions of the image
imageString = imageData(:);  % put all pixels in a datavector instead of matrix

imageBitStream = zeros(length(imageString)*bitsPerPixel,1); % imagebits
for m = 1:length(imageString)
    pixel = zeros(bitsPerPixel,1);
    val = imageString(m);
    for k = 1:bitsPerPixel
        pixel(k) = rem(val,2);
        val = floor(val/2);
    end
    imageBitStream(m*bitsPerPixel:-1:(m-1)*bitsPerPixel+1) = pixel; 
end

imageData = uint8(imageData);
