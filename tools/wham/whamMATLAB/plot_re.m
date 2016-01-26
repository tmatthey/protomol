% Load the results of the RE simulations
dih200K = makePositiveAngle(load('../200k/200k.dihedral.out'));
dih239K = makePositiveAngle(load('../239k/239k.dihedral.out'));
dih286K = makePositiveAngle(load('../286k/286k.dihedral.out'));
dih342K = makePositiveAngle(load('../342k/342k.dihedral.out'));
dih409K = makePositiveAngle(load('../409k/409k.dihedral.out'));
dih489K = makePositiveAngle(load('../489k/489k.dihedral.out'));
dih585K = makePositiveAngle(load('../585k/585k.dihedral.out'));
dih700K = makePositiveAngle(load('../700k/700k.dihedral.out'));
dih = [dih200K dih239K dih286K dih342K dih409K dih489K dih585K dih700K];
%dih = [dih342K dih409K dih489K dih585K dih700K];
dihE200K = load('../200k/200k.dihedralE.out');
dihE239K = load('../239k/239k.dihedralE.out');
dihE286K = load('../286k/286k.dihedralE.out');
dihE342K = load('../342k/342k.dihedralE.out');
dihE409K = load('../409k/409k.dihedralE.out');
dihE489K = load('../489k/489k.dihedralE.out');
dihE585K = load('../585k/585k.dihedralE.out');
dihE700K = load('../700k/700k.dihedralE.out');
dihE = [dihE200K dihE239K dihE286K dihE342K dihE409K dihE489K dihE585K dihE700K];
%dihE = [dihE342K dihE409K dihE489K dihE585K dihE700K];
pot200K = load('../200k/200k.potential.out');
pot239K = load('../239k/239k.potential.out');
pot286K = load('../286k/286k.potential.out');
pot342K = load('../342k/342k.potential.out');
pot409K = load('../409k/409k.potential.out');
pot489K = load('../489k/489k.potential.out');
pot585K = load('../585k/585k.potential.out');
pot700K = load('../700k/700k.potential.out');
pot = [pot200K pot239K pot286K pot342K pot409K pot489K pot585K pot700K];
%pot = [pot342K pot409K pot489K pot585K pot700K];

% Use the WHAM function to plot the results
T = [200; 239; 286; 342; 409; 489; 585; 700]; % Simulation temperatures
%T = [342; 409; 489; 585; 700]; % Simulation temperatures
refTemp = 300; % Reference temperature
V = zeros(length(dih700K(1:end,1)),1); % Biasing potentials (none)
wham_result = WHAM(T,dih,dihE,pot,V,refTemp);
hold on;

% Determine the analytical solution
x = 0:0.01:2*pi;
u = 1.6*(1 + cos(3*x)) + 0.6*(1+cos(x));
kb = 0.00198719;
T = 300;
B = 1/(kb * T);
y = exp(-B*u);
y = y / trapezoidRule(x,y);
x = 180 * x / pi;
plot(x,y,':k');

% Format the plot
axis([0 360 0 1.5]);
legend1 = legend(...
  {'REM Dist','Analytical Dist'},...
  'FontSize',12,...
  'Position',[0.6176 0.7862 0.2707 0.1174]);
hold off;

print -deps trial.eps;

% Determine the probability of each configuration
whamlen = length(wham_result(1:end,1));
ind0 = 1;
ind1 = ind0;
while wham_result(ind1,1) < 100
    ind1 = ind1 + 1;
end
ind2 = ind1;
while wham_result(ind2,1) < 250
    ind2 = ind2 + 1;
end
ind3 = whamlen;

wham1 = trapezoidRule(pi / 180 * wham_result(ind0:ind1,1),wham_result(ind0:ind1,2));
wham2 = trapezoidRule(pi / 180 * wham_result(ind1:ind2,1),wham_result(ind1:ind2,2));
wham3 = trapezoidRule(pi / 180 * wham_result(ind2:ind3,1),wham_result(ind2:ind3,2));
whamsum = wham1 + wham2 + wham3;
wham1 = wham1 / whamsum
wham2 = wham2 / whamsum
wham3 = wham3 / whamsum

% Determine the probability of each configuration
ylen = length(y);
ind0 = 1;
ind1 = ind0;
while x(ind1) < 120
    ind1 = ind1 + 1;
end
ind2 = ind1;
while x(ind2) < 240
    ind2 = ind2 + 1;
end
ind3 = ylen;

y1 = trapezoidRule(pi / 180 * x(ind0:ind1),y(ind0:ind1))
y2 = trapezoidRule(pi / 180 * x(ind1:ind2),y(ind1:ind2))
y3 = trapezoidRule(pi / 180 * x(ind2:ind3),y(ind2:ind3))
