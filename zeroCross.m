function zCross = zeroCross(row)

zCross=size(find(row==0),1);

for i=1:size(row,1)-1
    
    if row(i,1)<0 && row(i+1,1)>0
        zCross = zCross + 1;
    elseif row(i,1)>0 && row(i,1)<0
        zCross = zCross + 1;
    end
end
end
        

