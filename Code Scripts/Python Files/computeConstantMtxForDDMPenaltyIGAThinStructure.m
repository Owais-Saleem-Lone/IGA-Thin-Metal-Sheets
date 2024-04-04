
% Initialize the master penalty matrix of the multi-patch system
KPenalty = zeros(numDOFs, numDOFs);

% Get the interface quadrature structure
intC = propCoupling.intC;

% Check input
isDisplacementCoupling = false;
if isfield(propCoupling,'alphaD')
    if length(propCoupling.alphaD) ~= connections.No
        error('The number of the penalty parameters in propCoupling.alphaD (%d) must be equal to the number of connections (%d)', ...
            length(propCoupling.alphaD), connections.No);
    end
    isDisplacementCoupling = true;
end
isRotationCoupling = false;
if isfield(propCoupling,'alphaR')
    if length(propCoupling.alphaR) ~= connections.No
        error('The number of the penalty parameters in propCoupling.alphaR (%d) must be equal to the number of connections (%d)', ...
            length(propCoupling.alphaD), connections.No);
    end
    isRotationCoupling = true;
end

%% 1. Loop over all the connections in the multipatches and assemble the global penalty matrix
for iConnections = 1:connections.No
    %% 1i. Get the patch IDs
    
    % Patch 1 :
    % _________
    
    ID1 = connections.xiEtaCoup(iConnections, 1);
    
    % Patch 2 :
    % _________
    
    ID2 = connections.xiEtaCoup(iConnections, 2);
    
    %% 1ii. Get the penalty factors for the current patch pair
    if isDisplacementCoupling
        alphaDIJ = propCoupling.alphaD(iConnections, 1);
        if alphaDIJ == 0
            alphaDIJ = 'undefined';
        end
    else
        alphaDIJ = 'undefined';
    end
    if isRotationCoupling
        alphaRIJ = propCoupling.alphaR(iConnections, 1);
        if alphaRIJ == 0
            alphaRIJ = 'undefined';
        end
    else
        alphaRIJ = 'undefined';
    end
    
    %% 1iii. Get the coupling boundaries for both patches
    
    % Patch 1 :
    % _________
    
    BSplinePatches{ID1}.xicoup = connections.xiEtaCoup(iConnections, 3:4);
    BSplinePatches{ID1}.etacoup = connections.xiEtaCoup(iConnections, 5:6);
    
	% Patch 2 :
    % _________
    
    BSplinePatches{ID2}.xicoup = connections.xiEtaCoup(iConnections, 7:8);
    BSplinePatches{ID2}.etacoup = connections.xiEtaCoup(iConnections, 9:10);
    
    %% 1iv. Determine the interface orientation
    isSameOrientation = findSubdomainInterfaceOrientation ...
        (BSplinePatches{ID1}.p, BSplinePatches{ID1}.Xi, BSplinePatches{ID1}.q, ...
        BSplinePatches{ID1}.Eta, BSplinePatches{ID1}.CP, BSplinePatches{ID1}.isNURBS, ...
        BSplinePatches{ID1}.xicoup, BSplinePatches{ID1}.etacoup, ...
        BSplinePatches{ID2}.p,BSplinePatches{ID2}.Xi,BSplinePatches{ID2}.q, ...
        BSplinePatches{ID2}.Eta, BSplinePatches{ID2}.CP, BSplinePatches{ID2}.isNURBS, ...
        BSplinePatches{ID2}.xicoup, BSplinePatches{ID2}.etacoup);
    
    %% 1v. Compute the penalty contributions to the coupled system 
    [K1PenaltyDisplacements, K1PenaltyRotations, C1PenaltyDisplacements, ...
        C1PenaltyRotations, K2PenaltyDisplacements, K2PenaltyRotations] = ....
        computeDDMPenaltyMtcesIGAThinStructure ...
        (BSplinePatches{ID1}, BSplinePatches{ID2}, alphaDIJ, alphaRIJ, ...
        isSameOrientation, intC);
    
    % We can use the transpose of C1PenaltyDisplacements and
    % C1PenaltyRotations see below the assembly to the global coupling 
    % matrix for both patches
    
    %% 1vi. Assemble the penalty matrices to the master penalty stiffness/coupling matrix
    
    % Assemble the stiffness matrix contribution to patch 1
    KPenalty(BSplinePatches{ID1}.EFTPatches,BSplinePatches{ID1}.EFTPatches) = ...
        KPenalty(BSplinePatches{ID1}.EFTPatches,BSplinePatches{ID1}.EFTPatches) + ...
        K1PenaltyDisplacements + K1PenaltyRotations;
    
    % Assemble the stiffness matrix contribution to patch 2
    KPenalty(BSplinePatches{ID2}.EFTPatches,BSplinePatches{ID2}.EFTPatches) = ...
        KPenalty(BSplinePatches{ID2}.EFTPatches,BSplinePatches{ID2}.EFTPatches) + ...
        K2PenaltyDisplacements + K2PenaltyRotations;
    
    % Assemble to the coupling matrix contribution of the variational
    % problem in patch 1
    KPenalty(BSplinePatches{ID1}.EFTPatches,BSplinePatches{ID2}.EFTPatches) = ...
        KPenalty(BSplinePatches{ID1}.EFTPatches,BSplinePatches{ID2}.EFTPatches) + ...
        C1PenaltyDisplacements + C1PenaltyRotations;
    
    % Assemble to the coupling matrix contribution of the variational
    % problem in patch 2
    KPenalty(BSplinePatches{ID2}.EFTPatches,BSplinePatches{ID1}.EFTPatches) = ...
        KPenalty(BSplinePatches{ID2}.EFTPatches,BSplinePatches{ID1}.EFTPatches) + ...
        C1PenaltyDisplacements' + C1PenaltyRotations';
end

end
