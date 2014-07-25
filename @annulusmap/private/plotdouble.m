function plotdouble(pOuter,pInner)

if (nargin==1)
    pInner = pOuter{2};
    pOuter = pOuter{1};
end

turn_off_hold = ~ishold;

plot(pOuter)
hold on
plot(pInner)

if turn_off_hold, hold off, end;

end