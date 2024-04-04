%% 2. Get the number of DOFs for each patch and make an element freedom table for each Lagrange Multipliers field within the patch
for iPatches = 1:noPatches
    BSplinePatches{iPatches}.noDOFs = 3*BSplinePatches{iPatches}.noCPs;
end
%% 3. Compute the total number of DOFs for the multipatch structure
noDOFsPatches = 0;
for iPatches = 1:noPatches
    noDOFsPatches = noDOFsPatches + BSplinePatches{iPatches}.noDOFs;
end
%% 4. Create a freedom table for each patch in the multipatch geometry
for iPatches = 1:noPatches
    noDOFsPatch = BSplinePatches{iPatches}.noDOFs;
    if iPatches == 1
        BSplinePatches{iPatches}.EFTPatches = 1:noDOFsPatch;
    else
        BSplinePatches{iPatches}.EFTPatches = ...
            BSplinePatches{iPatches-1}.EFTPatches(length(BSplinePatches{iPatches-1}.EFTPatches)) + ...
            1:BSplinePatches{iPatches-1}.EFTPatches(length(BSplinePatches{iPatches-1}.EFTPatches)) + ...
            noDOFsPatch;
    end
end