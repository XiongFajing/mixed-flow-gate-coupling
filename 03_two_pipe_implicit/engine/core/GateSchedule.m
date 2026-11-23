function G=GateSchedule(t,G0,operation)

G=G0(:).*interp1(operation.times,operation.factors.', ...
    min(max(t,operation.times(1)),operation.times(end)),'linear').';
G=G(:);
end
