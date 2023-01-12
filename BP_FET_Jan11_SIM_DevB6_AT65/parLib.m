function pa=parLib()
% bvp4c acc
pa.acc_a=1e-6;
pa.acc_r=1e-6;
pa.nmax=1e3;

% physical constants
pa.ec = 1.602176634e-19; % C, elementary_charge
pa.bc = 1.380649e-23; % J.K^-1, Boltzmann_constant
pa.pc = 6.62607015e-34; % J.s, Planck_constant
pa.es = 3.9*8.8541878128e-12;% F.m^-1, epsilon_sio2
pa.fe = 9.1093837015e-31;% kg, free_electron_mass

% environment parameters
pa.TT = 300;% K

% geometrical parameters
pa.ch = 3.*1.6405017e-6;% m
lab=[6,5];
lab_tag=sprintf('%d-%d',lab);
switch lab_tag
case '3-1'
    wa=5.976757.*(2.36+2.6+2.77).^-1;
    wb=9.7091707.*(2.36+2.6+2.77).^-1;
    wl=(wb-wa).*(log(wb)-log(wa)).^-1;
case '3-2'
    wa=8.4840286.*(2.77).^-1;
    wb=9.7091707.*(2.77).^-1;
    wl=(wb-wa).*(log(wb)-log(wa)).^-1;
case '3-4'
    wa=10.992238.*(2.7406842).^-1;
    wb=10.641416.*(2.7406842).^-1;
    wl=(wb-wa).*(log(wb)-log(wa)).^-1;
case '5-3'
    wa=5.1948369.*(2.404594).^-1;
    wb=5.0005911.*(2.404594).^-1;
    wl=(wb-wa).*(log(wb)-log(wa)).^-1;
case '6-5'
    wa=5.9043575.*(1.2354405).^-1;
    wb=7.2295442.*(1.2354405).^-1;
    wl=(wb-wa).*(log(wb)-log(wa)).^-1;
end

pa.wl = wl;% W/L ratio
pa.ox = 20e-9/(9.1/3.9);% m

% material parameters
% pa.eg = 0.39*1; % eV
pa.eg = 2.20; % eV
pa.pn = 0.1*pa.eg; %eV
pa.pp = pa.eg-pa.pn; % eV

pa.me = 0.15; % m0
pa.mh = 0.14; % m0
pa.ue = 100e-4;% m^2.V^-1.s^-1
pa.uh = 100e-4;% m^2.V^-1.s^-1

% SRH equivalent mobility
pa.t0 = 25e-6;% s
% pa.t0 = 25e-15;% s

% threshold
% pa.threshold=0;

end