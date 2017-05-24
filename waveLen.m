function f = waveLen(row)

a=row(2:end,1)-row(1:end-1,1);

f = sum(abs(a))/size(row,1);

end