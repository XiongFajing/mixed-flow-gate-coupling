function s=MethodSettings(method)
s.method=char(method); s.timeCoupling='endpoint';
switch s.method
 case 'E1'
  s.number=1;s.gateLaw='classical';s.timeCoupling='explicit_old_new_opening';s.relaxation='none';s.label='E1 old heads';
 case 'R1'
  s.number=2;s.gateLaw='classical';s.relaxation='fixed';s.label='R1 raw/fixed';
 case 'R2'
  s.number=3;s.gateLaw='classical';s.relaxation='aitken';s.label='R2 raw/Aitken';
 case 'S1'
  s.number=4;s.gateLaw='smooth';s.relaxation='fixed';s.label='S1 smooth/fixed';
 case 'S2'
  s.number=5;s.gateLaw='smooth';s.relaxation='aitken';s.label='S2 smooth/Aitken';
 otherwise
  error('IFVM:Method','Use E1, R1, R2, S1 or S2.');
end
end
