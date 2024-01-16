classdef aSinxApp < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        TypeCoeffEditField       matlab.ui.control.NumericEditField
        TypeCoeffEditFieldLabel  matlab.ui.control.Label
        CoeffSpinner             matlab.ui.control.Spinner
        CoeffSpinnerLabel        matlab.ui.control.Label
        Button                   matlab.ui.control.Button
        UIAxes                   matlab.ui.control.UIAxes
    end

    
    properties (Access = private)
        Coeff=3; % Description
    end
    
    methods (Access = private)
        
        function updatePlot(app)
            ax= app.UIAxes;
            xV= linspace(0,2*pi,100);
            yV=app.Coeff*sin(xV);
            axis(ax,"tight")
            plot(ax,xV,yV);
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Button pushed function: Button
        function ButtonPushed(app, event)
            app.updatePlot;
        end

        % Value changed function: CoeffSpinner
        function CoeffSpinnerValueChanged(app, event)
            app.Coeff = app.CoeffSpinner.Value;
            %ButtonPushed(app,event)
        end

        % Value changed function: TypeCoeffEditField
        function TypeCoeffEditFieldValueChanged(app, event)
            app.Coeff = app.TypeCoeffEditField.Value;
            app.updatePlot;
            %ButtonPushed(app,event)
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';

            % Create UIAxes
            app.UIAxes = uiaxes(app.UIFigure);
            title(app.UIAxes, 'Plot of a*Sinx')
            xlabel(app.UIAxes, 'Angle (radians)')
            ylabel(app.UIAxes, 'Y = Sin(Angle)')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.Position = [304 268 300 185];

            % Create Button
            app.Button = uibutton(app.UIFigure, 'push');
            app.Button.ButtonPushedFcn = createCallbackFcn(app, @ButtonPushed, true);
            app.Button.Position = [271 63 100 23];

            % Create CoeffSpinnerLabel
            app.CoeffSpinnerLabel = uilabel(app.UIFigure);
            app.CoeffSpinnerLabel.HorizontalAlignment = 'right';
            app.CoeffSpinnerLabel.Position = [40 349 37 22];
            app.CoeffSpinnerLabel.Text = 'Coeff ';

            % Create CoeffSpinner
            app.CoeffSpinner = uispinner(app.UIFigure);
            app.CoeffSpinner.ValueChangedFcn = createCallbackFcn(app, @CoeffSpinnerValueChanged, true);
            app.CoeffSpinner.Position = [92 349 100 22];

            % Create TypeCoeffEditFieldLabel
            app.TypeCoeffEditFieldLabel = uilabel(app.UIFigure);
            app.TypeCoeffEditFieldLabel.HorizontalAlignment = 'right';
            app.TypeCoeffEditFieldLabel.Position = [23 290 63 22];
            app.TypeCoeffEditFieldLabel.Text = 'Type Coeff';

            % Create TypeCoeffEditField
            app.TypeCoeffEditField = uieditfield(app.UIFigure, 'numeric');
            app.TypeCoeffEditField.ValueChangedFcn = createCallbackFcn(app, @TypeCoeffEditFieldValueChanged, true);
            app.TypeCoeffEditField.Position = [101 290 100 22];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = aSinxApp

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