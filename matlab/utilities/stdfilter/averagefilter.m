function image=averagefilter(image, varargin)
%AVERAGEFILTER 2-D mean filtering.
%   B = AVERAGEFILTER(A) performs mean filtering of two dimensional 
%   matrix A with integral image method. Each output pixel contains 
%   the mean value of the 3-by-3 neighborhood around the corresponding
%   pixel in the input image. 
%
%   B = AVERAGEFILTER(A, [M N]) filters matrix A with M-by-N neighborhood.
%   M defines vertical window size and N defines horizontal window size. 
%   
%   B = AVERAGEFILTER(A, [M N], PADDING) filters matrix A with the 
%   predefinned padding. By default the matrix is padded with zeros to 
%   be compatible with IMFILTER. But then the borders may appear distorted.
%   To deal with border distortion the PADDING parameter can be either
%   set to a scalar or a string: 
%       'circular'    Pads with circular repetition of elements.
%       'replicate'   Repeats border elements of matrix A.
%       'symmetric'   Pads array with mirror reflections of itself. 
%
%   Comparison
%   ----------
%   There are different ways how to perform mean filtering in MATLAB. 
%   An effective way for small neighborhoods is to use IMFILTER:
%
%       I = imread('eight.tif');
%       meanFilter = fspecial('average', [3 3]);
%       J = imfilter(I, meanFilter);
%       figure, imshow(I), figure, imshow(J)
%
%   However, IMFILTER slows down with the increasing size of the 
%   neighborhood while AVERAGEFILTER processing time remains constant.
%   And once one of the neighborhood dimensions is over 7 pixels,
%   AVERAGEFILTER is faster. Anyway, both IMFILTER and AVERAGEFILTER give
%   the same results.
%
%   Remark
%   -------
%   The output class type is the same as the class type of input matrix A.
%
%   Example
%   -------
%       I = imread('eight.tif');
%       J = averagefilter(I, [5 5], 'replicate');
%       figure, imshow(I), figure, imshow(J)
%
%   See also IMFILTER, FSPECIAL, PADARRAY.

%   Contributed by Jan Motl (jan@motl.us)
%   $Revision: 1.3 $  $Date: 2013/05/04 16:58:01 $


% Parameter checking.
numvarargs = length(varargin);
if numvarargs > 2
    error('myfuns:somefun2Alt:TooManyInputs', ...
        'requires at most 2 optional inputs');
end
 
optargs = {[3 3] 0};            % set defaults for optional inputs
optargs(1:numvarargs) = varargin;
[window, padding] = optargs{:}; % use memorable variable names
m = window(1);
n = window(2);

if (ndims(image)~=2)            % check for color pictures
    display('The input image must be a two dimensional array.')
    display('Consider using rgb2gray or similar function.')
    return
end

% Initialization.
[rows columns] = size(image);   % size of the image

% Pad the image.
imageP  = padarray(image, [floor((m+1)/2) floor((n+1)/2)], padding, 'pre');
imagePP = padarray(imageP, [ceil((m-1)/2) ceil((n-1)/2)], padding, 'post');

% Always use double because uint8 would be too small.
imageD = double(imagePP);

% Calculate the integral image - the sum of numbers above and left.
t = cumsum(cumsum(imageD),2);

% Calculate the mean values from the look up table 't'.
imageI = t(1+m:rows+m, 1+n:columns+n) + t(1:rows, 1:columns)...
    - t(1+m:rows+m, 1:columns) - t(1:rows, 1+n:columns+n);

% Now each pixel contains sum of the window. But we want the average value.
imageI = imageI/(m*n);

% Return matrix in the original type class.
image = cast(imageI, class(image));
