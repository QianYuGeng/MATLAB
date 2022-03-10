function Smax = TurbineBlade_LSF(x)
% see: 'Thermal Stress Analysis of Jet Engine Turbine Blade' in matlab help
% 2021-07-06

    [N,~] = size(x);
    for j = 1:N
        E   = x(j,1)*1E9; % 227E9; % in Pa
        CTE = x(j,2)*1E-6; % 12.7E-6; % coefficient of thermal expansion, in 1/K
        nu  = x(j,3); % 0.27;   % 

        p1 = x(j,4)*1E3; % 5e5; %in Pa
        p2 = x(j,5)*1E3; % 4.5e5; %in Pa

        kapp = x(j,6); % 11.5; % in W/m/K thermal conductivity

        Tin  = x(j,7); % 150; The temperature of the interior cooling air
        Tout = x(j,8); % 1000; The temperature on the pressure and suction sides 

        tic
        %% Compute the stress caused only by this pressure.
        smodel = createpde('structural','static-solid');  % create a static structural model.
        importGeometry(smodel,'Blade.stl');
        % figure
        % pdegplot(smodel,'FaceLabels','on','FaceAlpha',0.8);
        msh = generateMesh(smodel,'Hmax',0.02);   % Note: it affects the results.
        structuralProperties(smodel,'YoungsModulus',E, 'PoissonsRatio',nu, 'CTE',CTE);
        structuralBC(smodel,'Face',3,'Constraint','fixed');

        structuralBoundaryLoad(smodel,'Face',11,'Pressure',p1); % Pressure side
        structuralBoundaryLoad(smodel,'Face',10,'Pressure',p2);  % Suction side 

        % Solve the structural problem.
        Rs = solve(smodel);
        % figure
        % pdeplot3D(smodel,'ColorMapData',Rs.VonMisesStress, 'Deformation',Rs.Displacement, 'DeformationScaleFactor',100);
        % view([116,25]);

        % toc

        %% Thermal Stress
        tmodel = createpde('thermal','steadystate');  % create a thermal model for steady-state thermal analysis.
        importGeometry(tmodel,'Blade.stl');
        tmodel.Mesh = msh;
        % kapp = 11.5; % in W/m/K
        thermalProperties(tmodel,'ThermalConductivity',kapp);

        % Interior cooling
        thermalBC(tmodel,'Face',[15 12 14], ...
                         'ConvectionCoefficient',30, ...
                         'AmbientTemperature',Tin);
        % Pressure side
        thermalBC(tmodel,'Face',11, ...
                         'ConvectionCoefficient',50, ...
                         'AmbientTemperature',Tout);
        % Suction side             
        thermalBC(tmodel,'Face',10, ...
                         'ConvectionCoefficient',40, ...
                         'AmbientTemperature',Tout);
        % Tip
        thermalBC(tmodel,'Face',13, ...
                         'ConvectionCoefficient',20, ... 
                         'AmbientTemperature',Tout);
        % Base (exposed to hot gases)
        thermalBC(tmodel,'Face',1, ...
                         'ConvectionCoefficient',40, ...
                         'AmbientTemperature',800);
        % Root in contact with hot gases
        thermalBC(tmodel,'Face',[6 9 8 2 7], ...
                         'ConvectionCoefficient',15, ...
                         'AmbientTemperature',400);

        % Root in contact with metal
        thermalBC(tmodel,'Face',[3 4 5], ...
                         'ConvectionCoefficient',1000, ...
                         'AmbientTemperature',300);

        Rt = solve(tmodel);

        % figure
        % pdeplot3D(tmodel,'ColorMapData',Rt.Temperature)
        % view([130,-20]);

        tsmodel = createpde('structural','static-solid'); % create a static structural model to compute the stress and deformation due to thermal expansion.
        importGeometry(tsmodel,'Blade.stl');
        tsmodel.Mesh = msh;
        structuralProperties(tsmodel,'YoungsModulus',E, ...
                                     'PoissonsRatio',nu, ...
                                     'CTE',CTE);
        tsmodel.ReferenceTemperature = 300; %in degrees C
        structuralBodyLoad(tsmodel,'Temperature',Rt);             
        structuralBC(tsmodel,'Face',3,'Constraint','fixed');

        Rts = solve(tsmodel);

        % figure('units','normalized','outerposition',[0 0 1 1]);
        % pdeplot3D(tsmodel,'ColorMapData',Rts.VonMisesStress, ...
        %                   'Deformation',Rts.Displacement, ...
        %                   'DeformationScaleFactor',100)
        % caxis([0, 200e6])
        % view([116,25]);
        % max(Rts.Displacement.Magnitude)

        %% Combined Pressure Loading and Thermal Stress
        structuralBoundaryLoad(tsmodel,'Face',11,'Pressure',p1); % Pressure side
        structuralBoundaryLoad(tsmodel,'Face',10,'Pressure',p2);  % Suction side 

        Rc = solve(tsmodel);

        % figure('units','normalized','outerposition',[0 0 1 1]);
        % pdeplot3D(tsmodel,'ColorMapData',Rc.VonMisesStress, ...
        %                  'Deformation',Rc.Displacement, ...
        %                  'DeformationScaleFactor',100)
        % caxis([0, 250e6])
        % view([110,20]);  % view([116,25]);

        % Evaluate the maximum stress and maximum displacement. 
        Smax(j,:) = max(Rc.VonMisesStress);
        toc
        % max(Rc.Displacement.Magnitude)
    end
    
    