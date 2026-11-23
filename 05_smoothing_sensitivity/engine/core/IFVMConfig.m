function cfg=IFVMConfig(caseID,profile)
if nargin<1,caseID='two_strong';end
if nargin<2,profile='screen';end
cfg=MakeCase(caseID,profile,USER_CONFIG());
end
