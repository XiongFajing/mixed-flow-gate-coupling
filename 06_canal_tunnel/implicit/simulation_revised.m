clc; clear; close all;
runTimer=tic;
scriptFolder=fileparts(mfilename('fullpath'));
addpath(scriptFolder);
outputFolder=fullfile(scriptFolder,'results');
if ~exist(outputFolder,'dir'), mkdir(outputFolder); end
makeFigures=false;
methodName='隐式边界（Aitken-Picard及Newton后备）';
g = 9.8;

L = [9000, 8000, 23000];
n_canal = length(L);
n = [0.015, 0.015, 0.015];
cs = [CityGate(8.3,5.96,0.01), CityGate(8.3,5.96,0.01), CityGate(8.4,5.88,0.01)];
S0 = [1/4300, 1/4300, 1/4200];

height = zeros(n_canal,1);
for i = 1:n_canal
    height(i) = cs(i).H + 0.5*cs(i).D;
end

Dl = 100;
n_grid = L./Dl;
x = cell(n_canal,1); x_real = cell(n_canal,1);
z = cell(n_canal,1); z_real = cell(n_canal,1);
for i = n_canal:-1:1
    x{i}  = (0.5 : 1 : n_grid(i)-0.5) * Dl;
    z{i}  = (L(i) - x{i}) * S0(i);
    if i == n_canal
        z_real{i} = z{i};
    else
        z_real{i} = z{i} + z_real{i+1}(1) + 0.5*Dl*S0(i);
    end
end
for i = 1:n_canal
    if i==1
        x_real{i} = x{i};
    else
        x_real{i} = x{i} + x_real{i-1}(end) + 0.5*Dl;
    end
end

y_design = [7.71, 7.71, 7.528];
Q_design = [135, 135, 135];
y = cell(n_canal,1); Q = cell(n_canal,1);
u = cell(n_canal,1); A = cell(n_canal,1);
c = cell(n_canal,1); a = zeros(n_canal,1);
for i = 1:n_canal
    y{i} = linspace(y_design(i)-0.5, y_design(i), n_grid(i));
    Q{i} = Q_design(i)*ones(1,n_grid(i));
    A{i} = cs(i).getAofY(y{i});
    u{i} = Q{i}./A{i};
    c{i} = cs(i).getCofY(y{i});
    a(i) = cs(i).a;
end

yall = cell(n_canal,1); Qall = cell(n_canal,1); Aall = cell(n_canal,1);
for i = 1:n_canal
    yall{i}  = [y{i}(1), y{i}, y{i}(end)];
    Aall{i} = cs(i).getAofY(yall{i});
    Qall{i}  = [Q{i}(1), Q{i}, Q{i}(end)];
end

yuu = 9.5;
ydd = 6.0;
Cd  = 0.83;

headEps = 0.02;
allowReverse = true;

solverOptions.tolAbs = 0.05;
solverOptions.tolRel = 1e-5;
solverOptions.flowScale = max(abs(Q_design));
solverOptions.closedOpeningTol = 1e-10;
solverOptions.picardMaxIter = 20;
solverOptions.picardOmegaInitial = 0.02;
solverOptions.picardOmegaMin = 1e-4;
solverOptions.picardOmegaMax = 0.5;
solverOptions.picardMaxFlowStep = 30;
solverOptions.newtonMaxIter = 8;
solverOptions.newtonFdAbs = 0.02;
solverOptions.newtonFdRel = 1e-4;
solverOptions.newtonMaxFlowStep = 50;
solverOptions.lineSearchMaxIter = 10;

u_max = max(cellfun(@max, u));
a_max = max(a);
Cr = 0.5;
Dt = Cr * Dl / max(u_max + a_max);
tf = 12 * 3600;
step_all = ceil(tf/Dt);
t_plot = Dt * (0:step_all);

gateStartTime = 2 * 3600;
delta_tu = 7200;
delta_td = 600;
delta_mid = [8, 8] * 3600;

[Hu0, Hd0] = GetGateHeads(y, yuu, ydd, S0, Dl);
G0 = zeros(n_canal+1,1);
Q_bc_previous = zeros(n_canal+1,1);
for i = 1:n_canal+1
    if i == 1
        sectionIndex = 1;
    elseif i == n_canal+1
        sectionIndex = n_canal;
    else
        sectionIndex = i-1;
    end

    targetQ = Q_design(sectionIndex);
    dH = Hu0(i) - Hd0(i);
    headFunction = dH/(dH^2 + headEps^2)^(1/4);
    requiredArea = abs(targetQ)/(Cd*sqrt(2*g)*abs(headFunction));

    if targetQ*dH <= 0 || requiredArea > cs(sectionIndex).Af
        error('闸门%d初始水头差与设计流量不相容，请检查初始条件。', i);
    end

    G0(i) = cs(sectionIndex).getYofA(requiredArea);
    Q_bc_previous(i) = targetQ;
end

Qout_cur = zeros(n_canal, 1);
Qout_plot = zeros(n_canal, step_all+1);

Q_plot = zeros(n_canal+1, step_all+1);
G_plot = zeros(n_canal+1, step_all+1);
hd_plot = zeros(n_canal, step_all+1);
GIFname = '水位过程线.gif';
step = 1;

stateTime=nan(1,step_all+1); stateTime(1)=0;
huHistory=nan(n_canal,step_all+1); hdHistory=huHistory;
quHistory=huHistory; qdHistory=huHistory;
for pool=1:n_canal
    huHistory(pool,1)=y{pool}(1); hdHistory(pool,1)=y{pool}(end);
    quHistory(pool,1)=Q{pool}(1); qdHistory(pool,1)=Q{pool}(end);
end
stepStart=nan(1,step_all); stepEnd=stepStart;
solveSeconds=stepStart; evalHistory=stepStart; iterHistory=stepStart;
newtonHistory=stepStart; residualHistory=stepStart;
runStatus='正常完成'; failureMessage=''; failedSolveSeconds=0;
loopTimer=tic;

t = 0;
while t < tf
    DtStep=min(Dt,tf-t);

    gateElapsed = max(0, t + DtStep - gateStartTime);
    G_u  = t2G(G0(1), 0, delta_tu,    gateElapsed);
    G_m1 = t2G(G0(2), 0, delta_mid(1),gateElapsed);
    G_m2 = t2G(G0(3), 0, delta_mid(2),gateElapsed);
    G_d  = t2G(G0(4), 0, delta_td,    gateElapsed);
    G_bc = [G_u; G_m1; G_m2; G_d];

    solveTimer=tic;
    try
    [coupledState, solverStats] = SolveGateCoupling( ...
        Q_bc_previous, G_bc, A, Q, y, ...
        n, S0, Dl, DtStep, cs, Qout_cur, Cd, ...
        yuu, ydd, headEps, allowReverse, solverOptions);
    catch ME
        failedSolveSeconds=toc(solveTimer);
        runStatus='求解失败，已保存此前完成的时间步';
        failureMessage=ME.message;
        warning('%s',failureMessage);
        break;
    end
    solveSeconds(step)=toc(solveTimer);
    evalHistory(step)=solverStats.evaluations;
    iterHistory(step)=solverStats.iterations;
    newtonHistory(step)=solverStats.newtonIterations;
    residualHistory(step)=solverStats.residualNorm;
    stepStart(step)=t; stepEnd(step)=t+DtStep;

    A = coupledState.A;
    Q = coupledState.Q;
    y = coupledState.y;
    Q_bc = coupledState.qGate;
    Q_bc_previous = Q_bc;

    G_plot(:,step) = G_bc;
    Q_plot(:,step) = Q_bc;
    Qout_plot(:,step) = Qout_cur;

    for pool = 1:n_canal
        yall{pool} = [y{pool}(1), y{pool}, y{pool}(end)];
        Aall{pool} = [A{pool}(1), A{pool}, A{pool}(end)];
        Qall{pool} = [ ...
            2*Q_bc(pool) - Q{pool}(1), ...
            Q{pool}, ...
            2*(Q_bc(pool+1) + Qout_cur(pool)) - Q{pool}(end)];
    end

    for pool = 1:n_canal
        hd_plot(pool, step) = yall{pool}(end-1);
    end

    stateTime(step+1)=t+DtStep;
    for pool=1:n_canal
        huHistory(pool,step+1)=y{pool}(1);
        hdHistory(pool,step+1)=y{pool}(end);
        quHistory(pool,step+1)=Q{pool}(1);
        qdHistory(pool,step+1)=Q{pool}(end);
    end
    step = step + 1;
    t = t + DtStep;
    if mod(step, 500) == 0
        fprintf('Progress: %.2f / %.2f h\n', t/3600, tf/3600);
    end
end

loopSeconds=toc(loopTimer);
completedSteps=step-1;
initializationSeconds=toc(runTimer)-loopSeconds;
exportTimer=tic;
ExportEngineering(outputFolder,methodName,makeFigures, ...
    stateTime(1:step),huHistory(:,1:step),hdHistory(:,1:step), ...
    quHistory(:,1:step),qdHistory(:,1:step), ...
    stepStart(1:completedSteps),stepEnd(1:completedSteps), ...
    G_plot(:,1:completedSteps),Q_plot(:,1:completedSteps), ...
    Qout_plot(:,1:completedSteps),solveSeconds(1:completedSteps), ...
    evalHistory(1:completedSteps),iterHistory(1:completedSteps), ...
    newtonHistory(1:completedSteps),residualHistory(1:completedSteps));
exportSeconds=toc(exportTimer);
WriteTiming(outputFolder,methodName,runStatus,failureMessage, ...
    Cr,Dl,Dt,tf,t,completedSteps,initializationSeconds,loopSeconds, ...
    sum(solveSeconds(1:completedSteps))+failedSolveSeconds, ...
    exportSeconds,sum(evalHistory(1:completedSteps)), ...
    sum(iterHistory(1:completedSteps)), ...
    sum(newtonHistory(1:completedSteps)>0));
fprintf('\n%s：完成 %d 步，终点 %.6f s\n',methodName,completedSteps,t);
fprintf('纯求解耗时 %.3f s；时间循环耗时 %.3f s；输出耗时 %.3f s\n', ...
    sum(solveSeconds(1:completedSteps))+failedSolveSeconds,loopSeconds,exportSeconds);
fprintf('输出目录：%s\n',outputFolder);
if ~isempty(failureMessage), fprintf('停止原因：%s\n',failureMessage); end
