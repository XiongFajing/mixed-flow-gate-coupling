function s=MethodSettings(method)
assert(strcmp(method,'E1'),'本包仅包含A2单管显式方案。');
s=struct('method','E1','number',1,'gateLaw','classical', ...
 'timeCoupling','explicit_old_new_opening','relaxation','none');
end
