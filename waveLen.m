function f = waveLen(row)

f = sum(abs(row(2:end,1)-row(1:end-1,1)))/size(row,1);

end