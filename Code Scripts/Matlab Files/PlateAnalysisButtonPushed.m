 function PlateAnalysisButtonPushed(app, event)
            cla(app.UIAxes5);
            propStr.t=app.ThicknessmEditField.Value;
            propStr.alpha=app.ShearFactoralphaEditField.Value;
            propStr.mxBar = 0;
            propStr.myBar = 0;
            propStr.E = app.YoungsModuluskNm2EditField_2.Value;
            propStr.nu = app.PoissonsRatioEditField_2.Value;
            propStr.G = propStr.E/2/(1 + propStr.nu);
            propStr.D = propStr.E*propStr.t^3/12/(1 - propStr.nu^2);
            propStr.pBar = -3e4*(propStr.t)^3;
            app.propStructure=propStr;
            propGeomHRef=app.changeProp;
            numBFsRef=[ numel(propGeomHRef.knotVct{1})-propGeomHRef.polOrd(1)-1;
                numel(propGeomHRef.knotVct{2})-propGeomHRef.polOrd(2)-1];
            numBFsRef=double(numBFsRef);
            numBFsRefTotal = prod(numBFsRef, "all");
            numCPsEl = prod(propGeomHRef.polOrd + 1);
            numDOFsEl = 3*numCPsEl;
            numDOFs = 3*numBFsRefTotal;
            stiffMtx = zeros(numDOFs);
            loadVct = zeros(numDOFs, 1);
            active_DOF=[];
            if app.propGeomTrim.isTrim==1
                xi_vecTrim=app.propGeomTrim.poly(:,2);
                eta_vecTrim=app.propGeomTrim.poly(:,1);
                trim_poly= polyshape(xi_vecTrim,eta_vecTrim);
            else
                trim_poly=[];
            end
            inactive_elements=0;

            for jj = propGeomHRef.polOrd(2, 1)+1:numel(propGeomHRef.knotVct{2})-propGeomHRef.polOrd(2, 1)-1
                for ii = propGeomHRef.polOrd(1, 1)+1:numel(propGeomHRef.knotVct{1})-propGeomHRef.polOrd(1, 1)-1
                    if diff(propGeomHRef.knotVct{1}(ii:ii + 1)) && diff(propGeomHRef.knotVct{2}(jj:jj + 1)) % Put this for trim plate as well
                        el = [propGeomHRef.knotVct{1}(ii) propGeomHRef.knotVct{1}(ii + 1) ...
                            propGeomHRef.knotVct{2}(jj) propGeomHRef.knotVct{2}(jj + 1)];
                        eft = zeros(1, numDOFsEl);
                        count = 1;
                        for ll = jj-propGeomHRef.polOrd(2, 1):jj
                            for kk = ii-propGeomHRef.polOrd(1, 1):ii
                                eft(count) = 3*((ll - 1)*numBFsRef(1, 1) + kk) - 2;
                                eft(count + 1) = 3*((ll - 1)*numBFsRef(1, 1) + kk) - 1;
                                eft(count + 2) = 3*((ll - 1)*numBFsRef(1, 1) + kk);
                                count = count + 3;
                            end
                        end
                        x_rect = [el(1); el(1); el(2); el(2); el(1)];
                        y_rect = [el(3); el(4); el(4); el(3); el(3)];
                        el_rect = polyshape(x_rect,y_rect);
                        if app.propGeomTrim.isTrim==1
                            intersectionPoly = intersect(el_rect, trim_poly);
                            TF = overlaps(el_rect,trim_poly);
                            if ~isempty(intersectionPoly.Vertices)

                                rectVertices=[x_rect,y_rect];
                                [xIntersect, yIntersect] = polyxpoly(x_rect,y_rect, xi_vecTrim,eta_vecTrim);
                                if ~isempty(xIntersect)

                                    active_DOF=[active_DOF eft];
                                    rectOutsideVertices = rectVertices(~inpolygon(rectVertices(:, 1), rectVertices(:, 2), xi_vecTrim, eta_vecTrim), :);
                                    new_meshVert= [rectOutsideVertices;[xIntersect yIntersect]];
                                    triangulation = delaunay(new_meshVert(:,1), new_meshVert(:,2));
                                    numGPs = 4;
                                    [xiGP, GWxi] = app.getGaussPointsAndWeightsOverUnitDomain(numGPs);
                                    [etaGP, GWeta] = app.getGaussPointsAndWeightsOverUnitDomain(numGPs);
                                    GP_T_par=zeros(numel(xiGP)*numel(etaGP)*size(triangulation,1),2);% 16 for each triangle
                                    GW_T_par=zeros(numel(xiGP)*numel(etaGP)*size(triangulation,1),1);
                                    det_Jac_vec=zeros(numel(xiGP)*numel(etaGP)*size(triangulation,1),1);
                                    counter=1;

                                    for tr= 1:size(triangulation,1)
                                        verTriangle=[new_meshVert(triangulation(tr,:),1) new_meshVert(triangulation(tr,:),2)];
                                        s1=verTriangle(1,1); s2=verTriangle(2,1);s3=verTriangle(3,1);
                                        t1=verTriangle(1,2); t2=verTriangle(2,2);t3=verTriangle(3,2);
                                        detJac_Triangle= det([-s1+s2 -t1+t2; -s1+s3 -t1+t3]);
                                        for gpa=1:numel(etaGP)
                                            for gpb=1:numel(xiGP)
                                                xiGTri= (1+xiGP(gpb))/2;
                                                etaGTri = (1-xiGP(gpb))*(1+etaGP(gpa))/4;
                                                gw_tr= (1-xiGP(gpb))*GWxi(gpb)*GWeta(gpa)/8;

                                                GP_T_par(counter,:)=GP_T_par(counter,:)+...
                                                    [(1-xiGTri-etaGTri)*s1+xiGTri*s2+etaGTri*s3 (1-xiGTri-etaGTri)*t1+xiGTri*t2+etaGTri*t3];
                                                GW_T_par(counter)=GW_T_par(counter)+gw_tr;
                                                det_Jac_vec(counter)=det_Jac_vec(counter)+detJac_Triangle;
                                                counter=counter+1;
                                            end
                                        end
                                    end
                                    [stiffMtxEl, loadVctEl] = app.computeIGAElStiffMtxandForceVctReissnerMindlinPlateTrim ...
                                        (propGeomHRef, det_Jac_vec, propStr,GP_T_par,GW_T_par);

                                else
                                    stiffMtxEl = zeros(numDOFsEl);
                                    loadVctEl = zeros(numDOFsEl, 1);
                                    inactive_elements=inactive_elements+1;
                                end
                            else
                                active_DOF=[active_DOF eft];

                                detJxiTildexi = (el(2) - el(1))*(el(4) - el(3))/4;

                                [stiffMtxEl, loadVctEl] = app.computeIGAElStiffMtxandForceVctReissnerMindlinPlate ...
                                    (propGeomHRef, el, detJxiTildexi, propStr);
                            end


                        else
                            active_DOF=[active_DOF eft];

                            detJxiTildexi = (el(2) - el(1))*(el(4) - el(3))/4;

                            [stiffMtxEl, loadVctEl] = app.computeIGAElStiffMtxandForceVctReissnerMindlinPlate ...
                                (propGeomHRef, el, detJxiTildexi, propStr);
                        end

                        stiffMtx(eft, eft) = stiffMtx(eft, eft) + stiffMtxEl;
                        loadVct(eft, 1) = loadVct(eft, 1) + loadVctEl;
                    end
                end
            end

                homDOFs=[];
                if app.CheckBox_1.Value==1 % clamp edge xi =0
                    xiSpan = [propGeomHRef.knotVct{1}(1) propGeomHRef.knotVct{1}(1)];
                    etaSpan = [propGeomHRef.knotVct{2}(1) propGeomHRef.knotVct{2}(end)];
                    numDOFsPerCP = 3;
                    idField = 1:3;
                    homDOFs = app.findDOFsNURBSSurface(homDOFs, propGeomHRef, xiSpan, etaSpan, numDOFsPerCP, idField);
                end
                if app.CheckBox_2.Value==1 % clamp edge xi =1
                    xiSpan = [propGeomHRef.knotVct{1}(end) propGeomHRef.knotVct{1}(end)];
                    etaSpan = [propGeomHRef.knotVct{2}(1) propGeomHRef.knotVct{2}(end)];
                    numDOFsPerCP = 3;
                    idField = 1:3;
                    homDOFs = app.findDOFsNURBSSurface(homDOFs, propGeomHRef, xiSpan, etaSpan, numDOFsPerCP, idField);
                end
                if app.CheckBox_3.Value==1 % clamp edge eta =0
                    xiSpan = [propGeomHRef.knotVct{1}(1) propGeomHRef.knotVct{1}(end)];
                    etaSpan = [propGeomHRef.knotVct{2}(1) propGeomHRef.knotVct{2}(1)];
                    numDOFsPerCP = 3;
                    idField = 1:3;
                    homDOFs = app.findDOFsNURBSSurface(homDOFs, propGeomHRef, xiSpan, etaSpan,numDOFsPerCP, idField);
                end


                if app.CheckBox_4.Value==1 % clamp edge eta =1
                    xiSpan = [propGeomHRef.knotVct{1}(1) propGeomHRef.knotVct{1}(end)];
                    etaSpan = [propGeomHRef.knotVct{2}(end) propGeomHRef.knotVct{2}(end)];
                    numDOFsPerCP = 3;
                    idField = 1:3;
                    homDOFs = app.findDOFsNURBSSurface(homDOFs, propGeomHRef, xiSpan,etaSpan, numDOFsPerCP, idField);
                end

                dofId = find(vecnorm(stiffMtx, 2, 2)==0);
                freeDOFs = 1:numDOFs;
                active_DOF=unique(active_DOF);
                inactive_DOF=setdiff(freeDOFs,active_DOF);
                freeDOFs(ismember(freeDOFs,homDOFs)) = [];
                freeDOFs(ismember(freeDOFs,inactive_DOF))=[];
                u = zeros(numDOFs, 1);
                u(freeDOFs) = stiffMtx(freeDOFs, freeDOFs)\loadVct(freeDOFs); % solving for displacements
                app.globDisp=u;

                propGeomHRefDef = propGeomHRef;
                propGeomHRefDef.controlPoints(3, :) = propGeomHRefDef.controlPoints(3, :) + transpose(1*u(1:3:end));
                app.deformProp=propGeomHRefDef;
                propFig.isPlotControlPolygon = false;
                propFig.isPlotKnots = true;
                propFig.isPlotNumsAndIds = false;
                propFig.isPlotJacobian=false;
                app.plotNURBSSurface(app.UIAxes5,propGeomHRefDef, 100, 100, numBFsRef, propFig);
                hold (app.UIAxes5,'on');
                grid(app.UIAxes5,'on');
                X0 = app.controlPoints(1,1);
                XLx = app.controlPoints(1,end);
                Y0 = app.controlPoints(2,1);
                YLy = app.controlPoints(2,end);
                syms wEx(X) betayEx(X) qxEx(X) Lsym psym Gsym Aqsym ...
                    Esym Isym tsym alphaScf
                wEx(X) = (psym*Lsym*X-psym*X^2/2)/Gsym/(Aqsym) - ...
                    (-psym*X^4/24+psym*Lsym*X^3/6-psym*Lsym^2*X^2/4)/Esym/Isym;
                wExSubs(X) = subs(wEx,{Lsym, psym, Gsym, Esym, Isym, Aqsym}, ...
                    {XLx, propStr.pBar, propStr.G, propStr.E, ...
                    YLy*(propStr.t^3)/12, propStr.alpha*YLy*propStr.t});

                % fsurf(app.UIAxes5,wExSubs, [X0 XLx Y0 YLy], FaceColor=[203 65 84]/255, ...
                %     FaceAlpha=0.5, LineStyle="none");
                % hold (app.UIAxes5,'off');
                % legend( app.UIAxes5,["IGA", "Analytical"], Location="best")
                %title(app.UIAxes5,"Plate's deformation using isogeometric elements")
                title(app.UIAxes5,"Deformation of Plate")
        end