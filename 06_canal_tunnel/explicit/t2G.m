function G = t2G(G_ini,G_end,delta_t,t)
n = length(t);
G = ones(n,1);
for i = 1:n
    if t(i) < delta_t
        G(i) = G_ini - (G_ini-G_end)*t(i)/delta_t;
    else
        G(i) = G_end;
    end
end
