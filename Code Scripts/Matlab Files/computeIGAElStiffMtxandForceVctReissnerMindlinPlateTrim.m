function [stiffMtxEl, loadVctEl] = computeIGAElStiffMtxandForceVctReissnerMindlinPlateTrim ...
                (app,propGeom, detJac_tri_vec, propStr,GP_T_par,GW_T_par)
            numBFsXi = numel(propGeom.knotVct{1}) - propGeom.polOrd(1, 1) - 1;
            numBFsEta = numel(propGeom.knotVct{2}) - propGeom.polOrd(2, 1) - 1;
            numCPsEl = prod(propGeom.polOrd + 1);
            numDOFsEl = 3*numCPsEl;
            stiffMtxEl = zeros(numDOFsEl);
            loadVctEl = zeros(numDOFsEl, 1);
            Cmtx = [propStr.alpha*propStr.G*propStr.t 0                                 0                    0                     0
                0                                 propStr.alpha*propStr.G*propStr.t 0                    0                     0
                0                                 0                                 propStr.D            propStr.nu*propStr.D  0
                0                                 0                                 propStr.nu*propStr.D propStr.D             0
                0                                 0                                 0                    0                     propStr.D/2*(1 - propStr.nu)];
            for jj = 1:size(GP_T_par,1)

                xi = GP_T_par(jj,1);
                eta = GP_T_par(jj,2);
                knotSpanXi = app.findKnotSpan(xi, propGeom.knotVct{1}, numBFsXi);
                knotSpanEta = app.findKnotSpan(eta, propGeom.knotVct{2}, numBFsEta);

                ordDervs = 1;
                dR = app.computeNURBSBasisFunctions2d(propGeom.knotVct, propGeom.polOrd, [xi; eta], propGeom.weights, ordDervs);
                rangeXi = knotSpanXi-propGeom.polOrd(1, 1):knotSpanXi;
                rangeEta = knotSpanEta-propGeom.polOrd(2, 1):knotSpanEta;
                [kk, ll] = ndgrid(rangeXi, rangeEta);
                dR = dR(:, (ll(:)-1)*numBFsXi + kk(:));
                controlPoints = propGeom.controlPoints(:, (ll(:)-1)*numBFsXi + kk(:));

                Rmtx = zeros(3, numDOFsEl);
                Rmtx(1, 3*(1:numCPsEl) - 2) = dR(1, 1:numCPsEl);
                Rmtx(2, 3*(1:numCPsEl) - 1) = dR(1, 1:numCPsEl);
                Rmtx(3, 3*(1:numCPsEl)) = dR(1, 1:numCPsEl);

                JmtxT = dR(2:3, :)*transpose(controlPoints(1:2, :));
                dRdX = JmtxT\dR(2:3, :);
                Bmtx = zeros(5, numDOFsEl);
                Bmtx(1, 3*(1:numCPsEl) - 2) = dRdX(1, 1:numCPsEl);
                Bmtx(1, 3*(1:numCPsEl)) = dR(1, 1:numCPsEl);
                Bmtx(2, 3*(1:numCPsEl)- 2) = dRdX(2, 1:numCPsEl);
                Bmtx(2, 3*(1:numCPsEl) - 1) = -dR(1, 1:numCPsEl);
                Bmtx(3, 3*(1:numCPsEl)) = dRdX(1, 1:numCPsEl);
                Bmtx(4, 3*(1:numCPsEl) - 1) = -dRdX(2, 1:numCPsEl);
                Bmtx(5, 3*(1:numCPsEl) - 1) = -dRdX(1, 1:numCPsEl);
                Bmtx(5, 3*(1:numCPsEl)) = dRdX(2, 1:numCPsEl);
                stiffMtxEl = stiffMtxEl + transpose(Bmtx)*Cmtx*Bmtx*detJac_tri_vec(jj)*det(JmtxT)*GW_T_par(jj);
                loadVctEl = loadVctEl + transpose(Rmtx)*[propStr.pBar; propStr.mxBar; propStr.myBar]* ...
                    detJac_tri_vec(jj)*det(JmtxT)*GW_T_par(jj);
            end
