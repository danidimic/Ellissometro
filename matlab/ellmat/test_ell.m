%% test ellipsometry simulation scripts

% define instrument parameters
parset(1240/632.8,1,pi/4,0.01)

% define a model for the interface
ellmod(1,1.5+i*0.02)

% plot the computed quantities
ellplot([0:0.01:pi/2]);

% build fake experimental data
theta_exp = [20:10:70]*pi/180;
[psi, delta] = ell(theta_exp);
psi_exp = psi+0.002*randn(size(psi))
delta_exp = delta+0.002*randn(size(delta))

% use the scripts for data fitting
ellfit([1.1 0],[1 0 NaN NaN],psi_exp,delta_exp,theta_exp)
