function deviation = stdfilter(image, varargin)
%STDFILTER 2-D standard deviation filtering.
%   B = STDFILTER(A) performs standard deviation filtering of two dimensional 
%   matrix A with integral image method. Each output pixel contains 
%   the standard deviation value of the 3-by-3 neighborhood around the corresponding
%   pixel in the input image.  
%
%   B = STDFILTER(A, [M N]) filters matrix A with M-by-N neighborhood.
%   M defines vertical window size and N defines horizontal window size. 
%
%   B = STDFILTER(A, [M N], CORRECTION) boolean value CORRECTION decides
%   whether to use corrected variance estimate (the degree of freedom is
%   decremented by one) or uncorected estimate (the degree of freedom is
%   equal to the sample size). 
%   The default setting is set to FALSE (use uncorrected estimate). 
%   
%   B = STDFILTER(A, [M N], CORRECTION, PADDING) filters matrix A with the 
%   predefinned padding. By default the matrix is padded with zeros to 
%   be compatible with STDFILT. But then the borders may appear distorted.
%   To deal with border distortion the PADDING parameter can be either
%   set to a scalar or a string: 
%       'circular'    Pads with circular repetition of elements.
%       'replicate'   Repeats border elements of matrix A.
%       'symmetric'   Pads array with mirror reflections of itself. 
%
%   Comparison
%   ----------
%   There are different ways how to perform standard deviation filtering in MATLAB. 
%   An effective way for small neighborhoods is to use STDFILT:
%
%       I = imread('eight.tif');
%       J = stdfilt(I);
%       figure, imshow(I), figure, imshow(J)
%
%   However, STDFILT slows down with the increasing size of the 
%   neighborhood while STDFILTER processing time remains constant.
%   And once one of the neighborhood dimensions is over 23 pixels,
%   STDFILT is faster. Anyway, both STDFILT and STDFILTER give
%   the same results.
%
%   Remark
%   -------
%   The output class type is the same as the class type of input matrix A.
%
%   Example
%   -------
%       I = imread('eight.tif');
%       J = stdfilter(I, [5 5], true, 'replicate');
%       figure, imshow(I), figure, imshow(J)
%
%   See also STDFILT, PADARRAY.

%   Contributed by Jan Motl (jan@motl.us)
%   $Revision: 1.0 $  $Date: 2013/08/05 16:58:01 $


% Parameter checking.
numvarargs = length(varargin);
if numvarargs > 4
    error('myfuns:somefun2Alt:TooManyInputs', ...
          'requires at most 4 optional inputs');
end
 
optargs = {[3 3] 'false' 0};             % set defaults for optional inputs
optargs(1:numvarargs) = varargin;
[window, corrected, padding] = optargs{:}; % use memorable variable names
m = window(1);
n = window(2);

if (ndims(image)~=2)                     % check for color pictures
    display('The input image must be a two dimensional array.')
    display('Consider using rgb2gray or similar function.')
    return
end

if ~isscalar(m)                          % check for window specification
    [m n] = size(m);
    display('STDFILT kernels are unsupported.')
    display('Continuing with assumption the kernel contains only ones...')
end

% Decide whether to use corrected variance estimate
if corrected
    normalisation = m*n/(m*n-1);    % denominator is 'n-1'
else
    normalisation = 1;              % denominator is 'n'
end

% Type conversion (integer sqrt is not supported on all MATLAB versions)
image = double(image);

% Mean value
mean = averagefilter(image, [m n], padding);

% Standard deviation
meanSquare = averagefilter(image.^2, [m n], padding);
deviation = (normalisation*(meanSquare - mean.^2)).^0.5;
