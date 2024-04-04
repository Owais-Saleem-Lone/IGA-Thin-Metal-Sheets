function testLinearKirchoffLoveShellMultipatchAnalysis(testCase)

%
%% Function definition
%
% Test the linear Kirchoff-Love shell isogeometric analysis over the 
% 4-patch semispherical shell with internal pressure solved with the
% penalty and the Lagrange Multipliers methods
%
% Function layout :
%
% 0. Read input
%
% 1. Define the multipatch geometry
%
% 2. Define the material constants
%
% 3. GUI
%
% 4. Refinement
%
% 5. Define boundary conditions
%
% 6. Create the patches and the Lagrange Multiplier fields
%
% 7. Solve the linear coupled system using the penalty method
%
% 8. Solve the linear coupled system using the Lagrange Multipliers method
%
% 9. Verify the results
%
%% Function main body

%% 0. Read input
                       
%% 1. Define the multipatch geometry

% Global variables
Radius = 0.075000000000000;

% Patch 1 :
% _________

% Polynomial degrees
p1 = 2;
q1 = 2;

% Knot vectors
Xi1 = [0 0 0 1 1 1];
Eta1 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP1(:, :, 1) = [0      0      0
                Radius Radius Radius
                Radius Radius Radius];
         
% y-coordinates
CP1(:, :, 2) = [-Radius -Radius 0
                -Radius -Radius 0
                0       0       0];
         
% z-coordinates
CP1(:, :, 3) = [0 Radius Radius
                0 Radius Radius
                0 0      0];
       
% Weights
weight = sqrt(2)/2;
CP1(:, :, 4) = [1      weight  1
                weight weight^2 weight 
                1      weight  1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS1 = false;
nxi1 = length(CP1(:, 1, 1));
neta1 = length(CP1(1, :, 1));
for i = 1:nxi1
    for j = 1:neta1
        if CP1(i, j, 4) ~= 1
            isNURBS1 = true;
            break;
        end
    end
    if isNURBS1
        break;
    end
end

% Patch 2 :
% _________

% Polynomial degrees
p2 = 2;
q2 = 2;

% Knot vectors
Xi2 = [0 0 0 1 1 1];
Eta2 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP2(:, :, 1) = [0      0      0
                Radius Radius Radius
                Radius Radius Radius];
         
% y-coordinates
CP2(:, :, 2) = [Radius Radius 0
                Radius Radius 0
                0      0      0];
         
% z-coordinates
CP2(:, :, 3) = [0 Radius Radius
                0 Radius Radius
                0 0      0];
       
% Weights
weight = sqrt(2)/2;
CP2(:, :, 4) = [1      weight  1
                weight weight^2 weight 
                1      weight  1];
          
% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS2 = false;
nxi2 = length(CP2(:, 1, 1));
neta2 = length(CP2(1, :, 1));
for i = 1:nxi2
    for j = 1:neta2
        if CP2(i, j, 4) ~= 1
            isNURBS2 = true;
            break;
        end
    end
    if isNURBS2
        break;
    end
end

% Patch 3 :
% _________

% Polynomial degrees
p3 = 2;
q3 = 2;

% Knot vectors
Xi3 = [0 0 0 1 1 1];
Eta3 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP3(:, :, 1) = [0       0       0
                -Radius -Radius -Radius
                -Radius -Radius -Radius];
         
% y-coordinates
CP3(:, :, 2) = [Radius Radius 0
                Radius Radius 0
                0      0      0];
         
% z-coordinates
CP3(:, :, 3) = [0 Radius Radius
                0 Radius Radius
                0 0      0];
       
% Weights
weight = sqrt(2)/2;
CP3(:, :, 4) = [1      weight  1
                weight weight^2 weight 
                1      weight  1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS3 = false;
nxi3 = length(CP3(:, 1, 1));
neta3 = length(CP3(1, :, 1));
for i = 1:nxi3
    for j = 1:neta3
        if CP3(i, j, 4) ~= 1
            isNURBS3 = 1;
            break;
        end
    end
    if isNURBS3
        break;
    end
end

% Patch 4 :
% _________

% Polynomial degrees
p4 = 2;
q4 = 2;

% Knot vectors
Xi4 = [0 0 0 1 1 1];
Eta4 = [0 0 0 1 1 1];

% Control Point coordinates

% x-coordinates
CP4(:, :, 1) = [0       0       0
                -Radius -Radius -Radius
                -Radius -Radius -Radius];
         
% y-coordinates
CP4(:, :, 2) = [-Radius -Radius 0
                -Radius -Radius 0
                0       0       0];
         
% z-coordinates
CP4(:, :, 3) = [0 Radius Radius
                0 Radius Radius
                0 0      0];
       
% Weights
weight = sqrt(2)/2;
CP4(:, :, 4) = [1      weight  1
                weight weight^2 weight 
                1      weight  1];

% Find whether the geometrical basis is a NURBS or a B-Spline
isNURBS4 = false;
nxi4 = length(CP4(:, 1, 1));
neta4 = length(CP4(1, :, 1));
for i = 1:nxi4
    for j = 1:neta4
        if CP4(i, j, 4) ~= 1
            isNURBS4 = true;
            break;
        end
    end
    if isNURBS4
        break;
    end
end



%% 2. Define the material constants

% general parameters
EYoung = 2.1e6;
nue = 0.0;
thickness = .03;

% Patch 1 :
% _________

% Young's modulus
parameters1.E = EYoung;

% Poisson ratio
parameters1.nue = nue;

% Thickness of the plate
parameters1.t = thickness;

% Patch 2 :
% _________

% Young's modulus
parameters2.E = EYoung;

% Poisson ratio
parameters2.nue = nue;

% Thickness of the plate
parameters2.t = thickness;

% Patch 3 :
% _________

% Young's modulus
parameters3.E = EYoung;

% Poisson ratio
parameters3.nue = nue;

% Thickness of the plate
parameters3.t = .03;

% Patch 4 :
% _________

% Young's modulus
parameters4.E = EYoung;

% Poisson ratio
parameters4.nue = nue;

% Thickness of the plate
parameters4.t = thickness;

%% 3. UI

% On the analysis
analysis.type = 'isogeometricKirchhoffLoveShellAnalysis';

% Define linear equation system solver
solve_LinearSystem = @solve_LinearSystemMatlabBackslashSolver;

% Integration scheme
% type = 'default' : default FGI integration element-wise
% type = 'user' : manual choice of the number of Gauss points

% Patch 1 :
% _________

int1.type = 'default';
if strcmp(int1.type, 'user')
    int1.xiNGP = 6;
    int1.etaNGP = 6;
    int1.xiNGPForLoad = 6;
    int1.etaNGPForLoad = 6;
    int1.nGPForLoad = 6;
    int1.nGPError = 12;
end

% Patch 2 :
% _________

int2.type = 'default';
if strcmp(int2.type, 'user')
    int2.xiNGP = 6;
    int2.etaNGP = 6;
    int2.xiNGPForLoad = 6;
    int2.etaNGPForLoad = 6;
    int2.nGPForLoad = 6;
    int2.nGPError = 12;
end

% Patch 3 :
% _________

int3.type = 'default';
if strcmp(int3.type, 'user')
    int3.xiNGP = 6;
    int3.etaNGP = 6;
    int3.xiNGPForLoad = 6;
    int3.etaNGPForLoad = 6;
    int3.nGPForLoad = 6;
    int3.nGPError = 12;
end

% Patch 4 :
% _________

int4.type = 'default';
if strcmp(int4.type, 'user')
    int4.xiNGP = 6;
    int4.etaNGP = 6;
    int4.xiNGPForLoad = 6;
    int4.etaNGPForLoad = 6;
    int4.nGPForLoad = 6;
    int4.nGPError = 12;
end

% Interface integration :
% _______________________

intC.type = 'default';
intC.method = 'Nitsche';
if strcmp(intC.type, 'user')
    if strcmp(intC.method, 'lagrangeMultipliers')
        intC.nGP1 = 12;
        intC.nGP2 = 12;
    else
        intC.nGP = 12;
    end
    intC.nGPError = 12;
end

% On the coupling

% Compute the material matrices for the membrane and the bending part
Dm = EYoung*thickness/(1 - nue^2)*...
      [1   nue  0
       nue 1    0
       0   0    (1 - nue)/2];
Db = thickness^3/12*Dm;

%% 4. Refinement

%%%%%%%%%%%%%%%%%%%%
% Degree elevation %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

tp1 = 0;
tq1 = 0;
[Xi1, Eta1, CP1, p1, q1] = degreeElevateBSplineSurface ...
    (p1, q1, Xi1, Eta1, CP1, tp1, tq1, '');

% Patch 2 :
% _________

tp2 = 0;
tq2 = 0;
[Xi2, Eta2, CP2, p2, q2] = degreeElevateBSplineSurface ...
    (p2, q2, Xi2, Eta2, CP2, tp2, tq2, '');

% Patch 3 :
% _________

tp3 = 0;
tq3 = 0;
[Xi3, Eta3, CP3, p3, q3] = degreeElevateBSplineSurface ...
    (p3, q3, Xi3, Eta3, CP3, tp3, tq3, '');

% Patch 4 :
% _________

tp4 = 0;
tq4 = 0;
[Xi4, Eta4, CP4, p4, q4] = degreeElevateBSplineSurface ...
    (p4, q4, Xi4, Eta4, CP4, tp4, tq4, '');



%%%%%%%%%%%%%%%%%%%%
% Knot insertion   %
%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

refXi1 = 3;
refEta1 = 3;
[Xi1, Eta1, CP1] = knotRefineUniformlyBSplineSurface ...
    (p1, Xi1, q1, Eta1, CP1, refXi1, refEta1, '');

% Patch 2 :
% _________

refXi2 = 4;
refEta2 = 4;
[Xi2, Eta2, CP2] = knotRefineUniformlyBSplineSurface ...
    (p2, Xi2, q2, Eta2, CP2, refXi2, refEta2, '');

% Patch 3 :
% _________

refXi3 = 3;
refEta3 = 3;
[Xi3, Eta3, CP3] = knotRefineUniformlyBSplineSurface ...
    (p3, Xi3, q3, Eta3, CP3, refXi3, refEta3, '');

% Patch 4 :
% _________

refXi4 = 4;
refEta4 = 4;
[Xi4, Eta4, CP4] = knotRefineUniformlyBSplineSurface ...
    (p4, Xi4, q4, Eta4, CP4, refXi4, refEta4, '');


%% 5. Define boundary conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dirichlet boundary conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs1 = [];
xisup1 = [0 1];
etasup1 = [0 0];
for dir = 1:3
    homDOFs1 = findDofs3D ...
        (homDOFs1, xisup1, etasup1, dir, CP1);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs1 = [];
valuesInhomDOFs1 = [];

% Patch 2 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs2 = [];
xisup2 = [0 1];
etasup2 = [0 0];
for dir = 1:3
    homDOFs2 = findDofs3D ...
        (homDOFs2, xisup2, etasup2, dir, CP2);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs2 = [];
valuesInhomDOFs2 = [];

% Patch 3 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs3 = [];
xisup3 = [0 1];
etasup3 = [0 0];
for dir = 1:3
    homDOFs3 = findDofs3D ...
        (homDOFs3, xisup3, etasup3, dir, CP3);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs3 = [];
valuesInhomDOFs3 = [];

% Patch 4 :
% _________

% Homogeneous Dirichlet boundary conditions
homDOFs4 = [];
xisup4 = [0 1];
etasup4 = [0 0];
for dir = 1:3
    homDOFs4 = findDofs3D ...
        (homDOFs4, xisup4, etasup4, dir, CP4);
end

% Inhomogeneous Dirichlet boundary conditions
inhomDOFs4 = [];
valuesInhomDOFs4 = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neumann boundary conditions   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% General parameter
loadAmplitude = 1e6;

% Patch 1 :
% _________

FAmp1 = loadAmplitude;
NBC1.noCnd = 1;
xib1 = [0 1];
etab1 = [0 1];
dirForce1 = 'normal';
NBC1.xiLoadExtension = {xib1};
NBC1.etaLoadExtension = {etab1};
NBC1.loadAmplitude = {FAmp1};
NBC1.loadDirection = {dirForce1};
NBC1.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC1.isFollower = true;
NBC1.isTimeDependent = false;

% Patch 2 :
% _________

FAmp2 = - loadAmplitude;
NBC2.noCnd = 1;
xib2 = [0 1];
etab2 = [0 1];
dirForce2 = 'normal';
NBC2.xiLoadExtension = {xib2};
NBC2.etaLoadExtension = {etab2};
NBC2.loadAmplitude = {FAmp2};
NBC2.loadDirection = {dirForce2};
NBC2.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC2.isFollower = true;
NBC2.isTimeDependent = false;

% Patch 3 :
% _________

FAmp3 = loadAmplitude;
NBC3.noCnd = 1;
xib3 = [0 1];
etab3 = [0 1];
dirForce3 = 'normal';
NBC3.xiLoadExtension = {xib3};
NBC3.etaLoadExtension = {etab3};
NBC3.loadAmplitude = {FAmp3};
NBC3.loadDirection = {dirForce3};
NBC3.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC3.isFollower = true;
NBC3.isTimeDependent = false;

% Patch 4 :
% _________

FAmp4 = - loadAmplitude;
NBC4.noCnd = 1;
xib4 = [0 1];
etab4 = [0 1];
dirForce4 = 'normal';
NBC4.xiLoadExtension = {xib4};
NBC4.etaLoadExtension = {etab4};
NBC4.loadAmplitude = {FAmp4};
NBC4.loadDirection = {dirForce4};
NBC4.computeLoadVct = {'computeLoadVctAreaIGAThinStructure'};
NBC4.isFollower = true;
NBC4.isTimeDependent = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interface parametrizations     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Patch 1 :
% _________

% Connection with patch 2:
xicoup12 = [0 1];   etacoup12 = [1 1];

% Connection with patch 4:
xicoup14 = [0 0];   etacoup14 = [0 1];

% Patch 2 :
% _________

% Connection with patch 1:
xicoup21 = [0 1];   etacoup21 = [1 1];

% Connection with patch 3:
xicoup23 = [0 0];   etacoup23 = [0 1];

% Patch 3 :
% _________

% Connection with patch 2:
xicoup32 = [0 0];   etacoup32 = [0 1];

% Connection with patch 4:
xicoup34 = [0 1];   etacoup34 = [1 1];

% Patch 4 :
% _________

% Connection with patch 3:
xicoup43 = [0 1];   etacoup43 = [1 1];

% Connection with patch 1:
xicoup41 = [0 0];   etacoup41 = [0 1];

% Define connections :
% ____________________

% Number of connections
noConnections = 4;

% Define connections by patch numbers
connections.No = noConnections;
connections.xiEtaCoup = zeros(noConnections, 10);
connections.xiEtaCoup(1, :) = [1 2 xicoup12 etacoup12 xicoup21 etacoup21];
connections.xiEtaCoup(2, :) = [2 3 xicoup23 etacoup23 xicoup32 etacoup32];
connections.xiEtaCoup(3, :) = [3 4 xicoup34 etacoup34 xicoup43 etacoup43];
connections.xiEtaCoup(4, :) = [1 4 xicoup14 etacoup14 xicoup41 etacoup41];

%% 6. Create the patches and the Lagrange Multiplier fields

% Patch 1 :
% _________

patch1 = fillUpPatch ...
    (analysis, p1, Xi1, q1, Eta1, CP1, isNURBS1, parameters1, homDOFs1, ....
    inhomDOFs1, valuesInhomDOFs1, [], [], NBC1, [], [], [], [], [], int1);

% Patch 2 :
% _________

patch2 = fillUpPatch ...
    (analysis, p2, Xi2, q2, Eta2, CP2, isNURBS2, parameters2, homDOFs2, ...
    inhomDOFs2, valuesInhomDOFs2, [], [], NBC2, [], [], [], [], [], int2);

% Patch 3 :
% _________

patch3 = fillUpPatch ...
    (analysis, p3, Xi3, q3, Eta3, CP3, isNURBS3, parameters3, homDOFs3, ...
    inhomDOFs3, valuesInhomDOFs3, [], [], NBC3, [], [], [], [], [], int3);

% Patch 4 :
% _________

patch4 = fillUpPatch ...
    (analysis, p4, Xi4, q4, Eta4, CP4, isNURBS4, parameters4, homDOFs4, ...
    inhomDOFs4, valuesInhomDOFs4, [], [], NBC4, [], [], [], [], [], int4);

% Collect all patches into an array :
% ___________________________________

BSplinePatches = {patch1 patch2 patch3 patch4};

%% 7. Solve the linear coupled system using the penalty method
propCouplingPenalty.alphaD = norm(Dm)*1e3*ones(4,1);
propCouplingPenalty.alphaR = norm(Db)*1e3*ones(4,1);
propCouplingPenalty.intC = intC;
[dHatPenaltyLinear, FCompletePenalty, ~, ~, ~, minElASizePenalty] = ...
    solve_DDMPenaltyIGAKirchhoffLoveShellMultipatchesLinear ...
    (BSplinePatches, connections, propCouplingPenalty, solve_LinearSystem, '');


end
