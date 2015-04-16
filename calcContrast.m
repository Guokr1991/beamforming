function CNR = calcContrast(env_les, env_bg)

u_les = mean(env_les);
v_les = var(env_les);

u_bg = mean(env_bg);
v_bg = var(env_bg);

CNR = abs(u_bg-u_les)/sqrt(v_bg+v_les);