function [CNR CR] = calcContrast(env_les, env_bg)
% Contrast calculations - Will Long. Latest revision: 4/20/15
% calculate contrast ratio and contrast to noise ratio given envelope
% detected echo signals from inside target and outside target

u_les = mean(env_les);
v_les = var(env_les);

u_bg = mean(env_bg);
v_bg = var(env_bg);

CNR = abs(u_bg-u_les)/sqrt(v_bg+v_les); % contrast to noise ratio
CR = 20.*log10(u_bg/u_les); % contrast ratio in dB