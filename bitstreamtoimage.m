% Construct BMP image from bitstream
% Input arguments
%   bitStream: bit stream in vector format
%   imageSize: dimensions of output image
%   bitsPerPixel: number of bits per pixel
% Output arguments
%   rxImage: image matrix

function rxImage = bitstreamtoimage(bitStream, imageSize, bitsPerPixel)

newPixelIndex=1; %index of the next received pixel in the vector rxDataValues
rxDataValues=zeros(imageSize(1)*imageSize(2),1);  %received data values should come in this vector
bitValues=2.^[bitsPerPixel-1:-1:0];    %the decimal value assigned to each bit position

for m=1:bitsPerPixel:length(bitStream)-bitsPerPixel+1
    if newPixelIndex<=length(rxDataValues)
        pixel=bitStream(m:m+bitsPerPixel-1);
        rxDataValues(newPixelIndex)=bitValues*pixel;
    end
    newPixelIndex=newPixelIndex+1;
end

% Construct image matrix out of image data stream
rxImage = reshape(rxDataValues,imageSize(1),imageSize(2));

% Convert to the correct datatype (necessary for the image command to
% use the correct value from the colormap)
if bitsPerPixel==1
    rxImage=logical(rxImage);
elseif bitsPerPixel==8
    rxImage=uint8(rxImage);
end
