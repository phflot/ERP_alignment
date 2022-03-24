function SqInteg = IntegralImage(v,xShift,yShift)

t = ImShift(v,xShift,yShift);
% Sqaured difference image
SqDiff = (v-t).^2;
% Cumulative sum  along rows
SqInteg = cumsum(SqDiff,1);
% Cumulative sum along columns
SqInteg = cumsum(SqInteg,2);
end