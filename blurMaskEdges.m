function mask = blurMaskEdges(mask0, kernSize)

% Smooth the mask to remove distracting edges
%
% From MAS
% Modified by Rachel Denision
% 2011 Nov 15

if nargin==1
    kernSize = 3;
end

mask = mask0;

% Pad on each side to start to deal with blurring at image edges
padWidth = kernSize*2;
for i=1:padWidth
    mask = growIm(mask);
end

skern0 = repmat(1, kernSize, kernSize) ./ kernSize^2;
mask = conv2(mask, skern0);
mask = conv2(mask, skern0);
mask = conv2(mask, skern0);

% Return mask to original size. Get rid of extra pixels due to initial
% padding and convolutions
for iDim = 1:2
    pad = size(mask,iDim) - size(mask0,iDim);
    if rem(pad,2)==1
        pad = pad+1;
    end
    inds{iDim} = pad/2+1:pad/2+size(mask0,iDim);
end
mask = mask(inds{1},:);
mask = mask(:,inds{2});

end

function imgrown = growIm(im)

im1 = [im(:,1) im im(:,end)]; % grow the left and right sides
im2 = [im1(1,:); im1; im1(end,:)]; % grow the top and bottom sides

imgrown = im2;

end