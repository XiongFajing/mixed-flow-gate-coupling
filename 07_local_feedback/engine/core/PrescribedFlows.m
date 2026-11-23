function q=PrescribedFlows(cfg,t,nGate)

q=nan(nGate,1);
if isfield(cfg,'prescribedFlowConstant') && ~isempty(cfg.prescribedFlowConstant)
    q=cfg.prescribedFlowConstant(:);
end
if isfield(cfg,'prescribedFlowTable') && ~isempty(cfg.prescribedFlowTable)
    tab=cfg.prescribedFlowTable;
    q=interp1(tab(:,1),tab(:,2:end),min(max(t,tab(1,1)),tab(end,1)),'linear').';
end
if isfield(cfg,'prescribedFlowFunction') && ~isempty(cfg.prescribedFlowFunction)
    q=cfg.prescribedFlowFunction(t); q=q(:);
end
if numel(q)~=nGate || any(isinf(q)), error('IFVM:Prescribed','Prescribed flow size/value invalid.'); end
if nGate>2 && any(isfinite(q(2:end-1)))
    error('IFVM:InternalFlow','Internal gate flows must be solved, not prescribed.');
end
end
