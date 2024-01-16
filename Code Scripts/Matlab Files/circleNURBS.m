classdef circleNURBS < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure               matlab.ui.Figure
        PolyOrderSpinner       matlab.ui.control.Spinner
        PolyOrderSpinnerLabel  matlab.ui.control.Label
        weightEditField        matlab.ui.control.NumericEditField
        weightEditFieldLabel   matlab.ui.control.Label
        DisplayGeometryButton  matlab.ui.control.Button
        UIAxes                 matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        polyOrd=2;
        weight
        knotVector=[0 0 0 1 1 2 2 3 3 4 4 4];
    end
    
    methods (Access = private)
        
        function updatePlot(app)
            ax= app.UIAxes;
            propGeomCircle.dim=1;
            propGeomCircle.polOrd = app.polyOrd;
            propGeomCircle.knotVct = app.knotVector;
            numBFsCircle = numel(propGeomCircle.knotVct) - propGeomCircle.polOrd - 1;
            w = app.weight;
            propGeomCircle.weights = [1 w 1 w 1 w 1 w 1];
            CP1 = [0,1, 0];
            CP2 = [1,1, 0];
            CP3 = [1,0, 0];
            CP4 = [1,-1, 0];
            CP5 = [0,-1, 0];
            CP6 = [-1,-1, 0];
            CP7 = [-1,0, 0];
            CP8 = [-1,1, 0];
            CP9 = [0,1, 0];
            propGeomCircle.controlPoints = [CP1(1) CP2(1) CP3(1) CP4(1) CP5(1) CP6(1) CP7(1) CP8(1) CP9(1)
                                CP1(2) CP2(2) CP3(2) CP4(2) CP5(2) CP6(2) CP7(2) CP8(2) CP9(2)
                                0      0      0      0      0      0      0      0      0];
            numPtsCircle=100;
            isPlotNumsAndIds = true;
            app.plotNURBSCurve(ax,propGeomCircle, numPtsCircle, isPlotNumsAndIds)

        end

        function plotNURBSCurve(app,ax,propGeom, numPts, isPlotNumsAndIds)
            numDerivs = 0;
            numBFs = size(propGeom.controlPoints, 2);
            minGeo = min(propGeom.controlPoints, [], "all");
            maxGeo = max(propGeom.controlPoints, [], "all");
            charLen = (minGeo - maxGeo)/20;
            xi = linspace(propGeom.knotVct(1), propGeom.knotVct(end), numPts);
            X = zeros(numel(xi), 3);
            for ii = 1:numel(xi)
                bfs = app.computeNURBSBasisFunctions1d(propGeom.knotVct, propGeom.polOrd, xi(ii), propGeom.weights, numDerivs);
                X(ii, :) = bfs*transpose(propGeom.controlPoints);
            end
            plot(ax,X(:, 1), X(:, 2), LineWidth=4, Color=([217, 218, 219] - 20)./255);
            hold(ax,'on')
            if isPlotNumsAndIds
                for ii = 1:numBFs
                    text(propGeom.controlPoints(1, ii) + charLen, propGeom.controlPoints(2, ii) + ...
                    charLen, sprintf("CP$_{%d}$", ii), Interpreter="latex");
                end
            end
            plot(ax,propGeom.controlPoints(1, :), propGeom.controlPoints(2, :), "r--", Marker='o', MarkerSize=10, ...
            MarkerFaceColor='r' , LineWidth=2)
        end

        function colmatw = computeNURBSBasisFunctions1d(app,knotVct, polOrd, xi, wmat, ordDerv)
            colmat = app.computeBSplineBasisFunctions1d(knotVct, polOrd, brk2knt(xi, ordDerv + 1));
            knotSpanId = app.findKnotSpan(xi, knotVct, numel(knotVct) - polOrd - 1);
            F = sum(colmat(:, knotSpanId-polOrd:knotSpanId) .* wmat(knotSpanId-polOrd:knotSpanId), 2);
            colmatw = zeros(size(colmat));
            for ii = 1:polOrd + 1
                for jj = 0:ordDerv
                    indexBF = knotSpanId - polOrd + ii - 1;
                    v = colmat(jj + 1, indexBF)*wmat(1, indexBF);
                    for kk = 1:jj
                        v = v - nchoosek(jj, kk)*F(kk + 1, 1)*colmatw(jj - kk + 1, indexBF);
                    end
                    colmatw(jj + 1, knotSpanId - polOrd + ii - 1) = v/F(1, 1);
                end
            end
        end

        function colmat = computeBSplineBasisFunctions1d(app,knotVct, polOrd, xi)
            colmat = spcol(knotVct, polOrd + 1, xi);

        end

        function idKnotSpan = findKnotSpan(app,xi, knotVct, numCPs)
            
            eps = 1e-7;
            numKnots = length(knotVct);
            polOrd = numKnots - numCPs - 1;
            if any(diff(knotVct(1, 1:polOrd + 1)) ~= 0) || any(diff(knotVct(1, end-polOrd:end)) ~= 0)
                mystr = sprintf("%1.1f ", knotVct);
                error(strcat("Provided knot vector knotVct = [ ", mystr, "]", "is not open as polynomial order is %d"), polOrd)
            end
            if norm(xi) < knotVct(1) - eps
                xi = knotVct(1);
            end
            if abs(xi - knotVct(numCPs + 1)) < eps
                idKnotSpan = numCPs;
                return
            end
            for idKnotSpan = 1:numKnots-1
                if xi < knotVct(idKnotSpan + 1)
                    return
                end
             end  

            error('The provided parametric location xi = %1.1f lies outside of the knot vector knotVct = [%1.1f]', xi, knotVct);


        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        function PolyOrderSpinnerValueChanged(app, event)
            app.polyOrd = app.PolyOrderSpinner.Value;
            
        end

        % Value changed function: weightEditField
        function weightEditFieldValueChanged(app, event)
            app.weight = app.weightEditField.Value;
            
        end
        function DisplayGeometryButtonPushed(app, event)
            app.updatePlot
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'IGA GEOMETRY';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            xlabel(app.UIAxes, 'X')
            ylabel(app.UIAxes, 'Y')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [304 268 300 185];

            % Create DisplayGeometryButton
            app.DisplayGeometryButton = uibutton(app.UIFigure, 'push');
            app.DisplayGeometryButton.ButtonPushedFcn = createCallbackFcn(app, @DisplayGeometryButtonPushed, true);
            app.DisplayGeometryButton.Position = [266 63 110 23];
            app.DisplayGeometryButton.Text = 'Display Geometry';

            % Create weightEditFieldLabel
            app.weightEditFieldLabel = uilabel(app.UIFigure);
            app.weightEditFieldLabel.HorizontalAlignment = 'right';
            app.weightEditFieldLabel.Position = [48 187 40 22];
            app.weightEditFieldLabel.Text = 'weight';

            % Create weightEditField
            app.weightEditField = uieditfield(app.UIFigure, 'numeric');
            app.weightEditField.Limits = [0 10];
            app.weightEditField.Position = [103 187 100 22];
            app.weightEditField.ValueDisplayFormat = '%.10f';
            app.weightEditField.ValueChangedFcn = createCallbackFcn(app, @weightEditFieldValueChanged, true);
            app.weightEditField.Value = 1;

            % Create PolyOrderSpinnerLabel
            app.PolyOrderSpinnerLabel = uilabel(app.UIFigure);
            app.PolyOrderSpinnerLabel.HorizontalAlignment = 'right';
            app.PolyOrderSpinnerLabel.Position = [36 239 59 22];
            app.PolyOrderSpinnerLabel.Text = 'PolyOrder';

            % Create PolyOrderSpinner
            app.PolyOrderSpinner = uispinner(app.UIFigure);
            app.PolyOrderSpinner.Limits = [1 5];
            app.PolyOrderSpinner.Position = [110 239 100 22];
            app.PolyOrderSpinner.Value = 2;
            app.PolyOrderSpinner.ValueChangedFcn = createCallbackFcn(app, @PolyOrderSpinnerValueChanged, true);

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = circleNURBS

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end