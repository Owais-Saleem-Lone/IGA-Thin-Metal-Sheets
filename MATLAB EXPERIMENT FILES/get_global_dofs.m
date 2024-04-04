%% 8. Rearrange the vectors containing the global numbering of the DOFs where Dirichlet boundary conditions are applied to account for the coupled system
sizePrevious = 0;
homDOFs = [];
inhomDOFs = [];
valuesInhomDOFs = [];
masterDOFs = [];
slaveDOFs = [];
for iPatches = 1:noPatches
    if iPatches ~= 1
        % Determine the number of the DOFs of the previous patches
        sizePrevious = sizePrevious + BSplinePatches{iPatches-1}.noDOFs;
        % Add the numbering DOFs where homogeneous Dirichlet boundary
        % conditions are applied from the patch level
        homDOFsPatch = sizePrevious + BSplinePatches{iPatches}.homDOFs;
        homDOFs = mergesorted(homDOFs,homDOFsPatch);

    else
        homDOFs = BSplinePatches{iPatches}.homDOFs;
    end
end
freeDOFs = 1:noDOFsPatches;
freeDOFs(ismember(freeDOFs,homDOFs)) = [];
freeDOFs(ismember(freeDOFs,inhomDOFs)) = [];


%% 10. Create a DOF numbering for each local patch

% For the patches :
% _________________

for iPatches = 1:noPatches
    % Get the number of Control Points in xi-direction
    nxi = length(BSplinePatches{iPatches}.CP(:,1,1));
    
    % Get the number of Control Points in eta-direction
    neta = length(BSplinePatches{iPatches}.CP(1,:,1));
    
    % Initialize the DOF numbering array
    BSplinePatches{iPatches}.DOFNumbering = zeros(nxi,neta,3);
    
    % Compute the entries of the DOF numbering array
    k = 1;
    for cpj = 1:neta
        for cpi = 1:nxi
            BSplinePatches{iPatches}.DOFNumbering(cpi,cpj,1) = k;
            BSplinePatches{iPatches}.DOFNumbering(cpi,cpj,2) = k + 1;
            BSplinePatches{iPatches}.DOFNumbering(cpi,cpj,3) = k + 2;

            % Update counter
            k = k + 3;
        end
    end
end