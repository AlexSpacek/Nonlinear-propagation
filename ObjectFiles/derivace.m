function f = derivace(aq,bq,x)
  a2 = abs(aq - x);
  [m,ind] = min(a2);
  if ind <= 3
      ind = ind + 3;
  elseif ind >= length(aq)-3
      ind = ind -3;
  end    

f = ( bq(ind-2)-8*bq(ind-1)+8*bq(ind+1)-bq(ind+2) ) / (12*(aq(ind)-aq(ind-1)));