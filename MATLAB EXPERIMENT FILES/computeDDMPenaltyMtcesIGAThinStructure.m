function [KpDispI, KpRotI, CpDispI, CpRotI, KpDispJ, KpRotJ] = ...
    computeDDMPenaltyMtcesIGAThinStructure...
    (patchI, patchJ, alphaD, alphaR, haveSameOrientation, int)
%% Licensing
%
% License:         BSD License
%                  cane Multiphysics default license: cane/license.txt
%
% Main authors:    Andreas Apostolatos
%
%% Function documentation
%
%  Returns the stiffness and the coupling matrices for the penalty
%  decomposition method applied to the multipatch Kirchhoff-Love shell. The
%  coupling matrices are only return for patch I cause the respective
%  contributions in patch J are symmetric.
%
%                  Input : 
%          patchI,patchJ : The B-Spline patches which are sharing an
%                          interface
%                 alphaD : Penalty factor for the displacement coupling
%                 alphaR : Penalty factor for the rotation coupling
%    haveSameOrientation : Flag on whether the couling surfaces of the 
%                          shells are oriented in the same direction over 
%                          the coupling interface
%                   int  : On the interface quadrature :
%                           .type : 'default' or 'user'
%                          .noGPs : Number of Gauss Points
%
%                Output :
%        KpDispI,KpRotI : Displacement and rotation coupling contribution
%                         to the stiffness matrix for patch 1
%        KpDispJ,KpRotJ : Displacement and rotation coupling contribution
%                         to the stiffness matrix for patch 2
%        CpDispI,CpRotI : Displacement and rotation coupling contribution
%                         to the coupling matrix for patch 1
%
% Function layout :
%
% 0. Read input
%
% 1. Get the running and the fixed parameters on the patch interface and the coupling region
%
% 2. Compute the merged knot vector from both patches over the interface
%
% 3. Issue Gauss Point coordinates and weights
%
% 4. Loop over the elements on the coupling interface
% ->
%    4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
%
%  4iii. Loop over the Gauss points
%  ->
%        4iii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
%
%        4iii.2. Compute the NURBS basis functions
%
%        4iii.3. Create the element freedom tables for both patches
%
%        4iii.4. Compute the covariant base vectors
%
%        4iii.5. Compute the surface normal vectors
%
%        4iii.6. Compute the derivatives of the surface normal vectors
%
%        4iii.7. Compute the normal to the boundary vector
%
%        4iii.8. Compute the normal to the boundary vector
%
%        4iii.9. Compute the covariant metric coefficients
%
%        4iii.10. Compute the contravariant base vectors
%
%        4iii.11. Compute the basis functions matrix and their derivatives and the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
%
%        4iii.12. Transform the normal and the tangent vectors to the covariant bases
%
%        4iii.13. Compute the curvature coefficients
%
%        4iii.14. Compute the B-operator matrices for the rotations
%
%        4iii.15. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
%
%        4iii.16. Compute the element length at the GP
%
%        4iii.17. Compute the element stiffness matrix contributions on the GP and add them to the global matrices
%
%        4iii.18. Compute the element coupling matrices at the GP
%  <-
% <-
% 
%% Function main body

%% 0. Read input

% Initialize a tolerance for the geometrical map
tolDet = 1e-8;

% For patch I :
% _____________

% Reassign the analysis arrays
pI = patchI.p;
qI = patchI.q;
XiI = patchI.Xi;
EtaI = patchI.Eta;
CPI = patchI.CP;
isNURBSI = patchI.isNURBS;
xicoupI = patchI.xicoup;
etacoupI = patchI.etacoup;
numDOFsI = patchI.noDOFs;

% Get the DOF numbering
DOFNumberingI = patchI.DOFNumbering;

% Number of Control Points in xi-,eta- directions
numCPs_xiI = length(CPI(:, 1, 1));
numCPs_etaI = length(CPI(1, :, 1));

% Number of local DOFs
numCPsElI = (pI + 1)*(qI + 1);
numDOFsElI = 3*numCPsElI;

% For patch J :
% _____________

pJ = patchJ.p;
qJ = patchJ.q;
XiJ = patchJ.Xi;
EtaJ = patchJ.Eta;
CPJ = patchJ.CP;
isNURBSJ = patchJ.isNURBS;
xicoupJ = patchJ.xicoup;
etacoupJ = patchJ.etacoup;
numDOFsJ = patchJ.noDOFs;

% Get the DOF numbering
DOFNumberingJ = patchJ.DOFNumbering;

% Number of Control Points in xi-,eta- directions
numCPs_xiJ = length(CPJ(:, 1, 1));
numCPs_etaJ = length(CPJ(1, :, 1));

% Number of local DOFs
numCPsElJ = (pJ + 1)*(qJ + 1);
numDOFsElJ = 3*numCPsElJ;

% Initialize the element freedom tables
EFTI = zeros(1, numDOFsElI);
EFTJ = zeros(1, numDOFsElJ);

% Initialize auxiliary arrays
BDisplacementsGCI = zeros(3, numDOFsElI);
BDisplacementsGCJ = zeros(3, numDOFsElJ);


% Initialize the output arrays
KpDispI = zeros(numDOFsI, numDOFsI);
KpDispJ = zeros(numDOFsJ, numDOFsJ);
CpDispI = zeros(numDOFsI, numDOFsJ);

%% 1. Get the running and the fixed parameters on the patch interface and the coupling region



%% 3. Issue Gauss Point coordinates and weights
if isstruct(int)
    if isfield(int, 'type')
        if strcmp(int.type, 'default')
            if isOnXiI
                pDegreeI = pI + 1;
            else
                pDegreeI = qI + 1;
            end
            if isOnXiJ
                pDegreeJ = pJ + 1;
            else
                pDegreeJ = qJ + 1;
            end
            numGPs = ceil((pDegreeI + pDegreeJ + 1)/2);
        elseif strcmp(int.type,'user')
            numGPs = int.noGPs;
        else
            error('int must define the to type to be either "default" or "user"');
        end
    else
        error('int must define variable type');
    end
else
    error('int must be a structure');
end
[GP,GW] = getGaussPointsAndWeightsOverUnitDomain(numGPs);
GP = fliplr(GP);
GW = fliplr(GW);

%% 4. Loop over the elements on the coupling interface
for i = 1:length(couplingRegionOnKnotVector)-1
    if couplingRegionOnKnotVector(i) ~= couplingRegionOnKnotVector(i + 1)
        %% 4i. Compute the determinant of the Jacobian of the transformation from the parent to the integation domain
        detJxizeta = (couplingRegionOnKnotVector(i + 1) - couplingRegionOnKnotVector(i))/2;

        %% 4iii. Loop over the Gauss points
        for j = 1:numGPs
            %% 4iii.1. Transform the Gauss Point coordinate from the bi-unit interval to the knot span
            xiEta = ((1 - GP(j))*couplingRegionOnKnotVector(i) + (1 + GP(j))*couplingRegionOnKnotVector(i + 1))/2;

            %% 4iii.2. Compute the NURBS basis functions
            

            dRI = computeIGABasisFunctionsAndDerivativesForSurface ...
                (xiSpanI, pI, xiI, XiI, etaSpanI, qI, etaI, EtaI, CPI, isNURBSI, 2);
            

            dRJ = computeIGABasisFunctionsAndDerivativesForSurface...
                (xiSpanJ, pJ, xiJ, XiJ, etaSpanJ, qJ, etaJ, EtaJ, CPJ, isNURBSJ, 2);
            
            %% 4iii.3. Create the element freedom tables
            
            % For patch I :
            % _____________
            
            % Initialize of the counter
            rI = 1;

            % Relation global-local DoFs
            for cpj = etaSpanI - qI:etaSpanI
                for cpi = xiSpanI - pI:xiSpanI
                    EFTI(rI) = DOFNumberingI(cpi, cpj, 1);
                    EFTI(rI + 1) = DOFNumberingI(cpi, cpj, 2);
                    EFTI(rI + 2) = DOFNumberingI(cpi, cpj, 3);

                    % update counter
                    rI = rI + 3;
                end
            end
            
            % For patch J :
            % _____________
                        
            % Initialize of the counter
            rJ = 1;

            % Relation global-local DoFs
            for cpj = etaSpanJ - qJ:etaSpanJ
                for cpi = xiSpanJ - pJ:xiSpanJ
                    EFTJ(rJ) = DOFNumberingJ(cpi, cpj, 1);
                    EFTJ(rJ + 1) = DOFNumberingJ(cpi, cpj, 2);
                    EFTJ(rJ + 2) = DOFNumberingJ(cpi, cpj, 3);

                    % update counter
                    rJ = rJ + 3;
                end
            end
        


            
            %% 4iii.10. Compute the basis functions matrix and their derivatives and the determinant of the Jacobian to the transformation from the physical space (x-y) to the NURBS parameter space (xi-eta)
            
            % For patch I :
            % _____________

            % initialize counter
            kI = 0;
            
            % Loop over all the non-zero contributions at the span
            % under study
            for c = 0:qI
                for b = 0:pI
                    % Update counter
                    kI = kI + 1;
                    
                    % Matrix containing the basis functions
                    BDisplacementsGCI(1,3*kI - 2) = dRI(kI, 1);
                    BDisplacementsGCI(2,3*kI - 1) = dRI(kI, 1);
                    BDisplacementsGCI(3,3*kI) = dRI(kI, 1);
                     
                end
            end
            
            % For patch J :
            % _____________
            
            % initialize counter
            kJ = 0;
            
            % Loop over all the non-zero contributions at the span
            % under study
            for c = 0:qJ
                for b = 0:pJ
                    % Update counter
                    kJ = kJ + 1;
                    
                    % Matrix containing the basis functions
                    BDisplacementsGCJ(1, 3*kJ - 2) = dRJ(kJ, 1);
                    BDisplacementsGCJ(2, 3*kJ - 1) = dRJ(kJ, 1);
                    BDisplacementsGCJ(3, 3*kJ) = dRJ(kJ, 1);

                    % Matrix containing the derivatives of the basis functions
 
                end
            end

            %% 4iii.14. Compute the determinant of the Jacobian of the transformation from the physical to the parent domain on the GP
            if isOnXiI
                detJxxi = norm(dA1I(:, 1));
            else
                detJxxi = norm(dA2I(:, 1));
            end
            
            %% 4iii.16. Compute the element length at the GP
            elementLengthOnGP = detJxxi*detJxizeta*GW(j);
            
            %% 4iii.17. Compute the element stiffness matrix contributions on the GP and add them to the global matrices
            
            % For patch I :
            % _____________
            
            % Compute the displacement stiffness matrix
            if ~strcmp(alphaD, 'undefined')
                KpDispI(EFTI, EFTI) = KpDispI(EFTI, EFTI) + ...
                    alphaD*(BDisplacementsGCI'*BDisplacementsGCI)*elementLengthOnGP;
            end
            
            
            % For patch J :
            % _____________
            
            % Compute the displacement stiffness matrix
            if ~strcmp(alphaD, 'undefined')
                KpDispJ(EFTJ, EFTJ) = KpDispJ(EFTJ, EFTJ) + ...
                    alphaD*(BDisplacementsGCJ'*BDisplacementsGCJ)*elementLengthOnGP;
            end
            
     
            
            %% 4iii.18. Compute the element coupling matrices at the GP
        
            % Compute the displacement coupling matrix
            if ~strcmp(alphaD, 'undefined')
                CpDispI(EFTI, EFTJ) = CpDispI(EFTI, EFTJ) - ...
                    alphaD*(BDisplacementsGCI'*BDisplacementsGCJ)*elementLengthOnGP;
            end

           
        end
    end
end

end
