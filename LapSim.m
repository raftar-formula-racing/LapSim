function Log = trial(Track, Car, Param)
    % Run Setup.m before this
    % Have Track as a double array with X, Y and Z as columns in that order

    % TODO
    % Implement Voltage vs SoC
    % Implement going over the 102Nm limit for some instants
    % Implement Load transfer and tire data
    
    NoOfDataPoints = length(Track);
    iter = 1;
    
    Log = struct();
    Log.X = Track(:,1);
    Log.Y = Track(:,2);
    Log.RoC = Track(:,3);
    Log.Vel = 12.*ones(NoOfDataPoints,1); % V_init = 18;
    Log.LongAcc = zeros(NoOfDataPoints,1);
    Log.LatAcc = zeros(NoOfDataPoints,1);
    Log.Torque = zeros(NoOfDataPoints,1);
    Log.dist = zeros(NoOfDataPoints,1);
    Log.E_useful = zeros(NoOfDataPoints, 1);
    Log.E_drag = zeros(NoOfDataPoints, 1);
    Log.E_rollres = zeros(NoOfDataPoints, 1);
    
    E = 0;

    while(iter<=Param.NoOfLaps)
        dist = 0;
        for i = 1:NoOfDataPoints-2
            curr_latacc = Log.Vel(i)^2/Log.RoC(i);
            if curr_latacc >= LatGlimit(Log.Vel(i),Car)
                Log.Vel(i) = 0.99*sqrt(Log.RoC(i)*LatGlimit(Log.Vel(i),Car));
                curr_latacc = Log.Vel(i)^2/Log.RoC(i);
            end
            avlbl_decel = -sqrt((1-(curr_latacc/LatGlimit(Log.Vel(i),Car))^2))*Car.decel_limit;
            Lookahead = Log.Vel(i)^2/(2*abs(avlbl_decel));
            Lookahead_idx = floor(Lookahead/(Param.LookaheadFactor)) + i;
            if Lookahead_idx > NoOfDataPoints -1
                Lookahead_idx = max(1,rem(Lookahead_idx,NoOfDataPoints-1));
            end
        
            if Log.Vel(Lookahead_idx) < Log.Vel(i)
                dv_dt = avlbl_decel;
                Torque = 0;
            else
                [dv_dt,vmax, Torque] = accel(Log.Vel(i), Car, curr_latacc);
            end
            Log.LongAcc(i) = dv_dt;
            s = sqrt((Track(i,1) - Track(i+1,1))^2 + (Track(i,2) - Track(i+1,2))^2);
            dist = dist + s;
            Log.dist(i) = dist;
            if (Log.Vel(i) < -2*dv_dt*s)
                Log.Vel(i+1) = 0;
            else
                Log.Vel(i+1) = max(min(sqrt(Log.Vel(i)^2 + 2*dv_dt*s), sqrt(Log.RoC(i+1)*LatGlimit(Log.Vel(i+1),Car))),0);
            end
            if Log.Vel(i+1) > vmax
                Log.Vel(i+1) = vmax;
            end
            Log.Torque(i) = Torque;
            Log.LatAcc(i) = curr_latacc;
        end
        iter = iter + 1;
    end
    
    Log.time =  dist/mean(Log.Vel(:));
    
    % Filter Velocities
    Log.Vel = lowpass(Log.Vel,0.01);
    % Energy Calculations
    for i = 1:NoOfDataPoints-1
        Log.E_useful(i) = max(0.5*Car.Mass*(Log.Vel(i)^2 - Log.Vel(i+1)^2),0); % Useful Energy
        Log.E_drag(i) = 0.5*1.23*Car.FrntArea*Car.Cd*(Log.Vel(i))^2*s; % Energy lost to drag
        Log.E_rollres(i) = 9.81*0.05*Car.Mass*s; % Energy lost to Rolling Resistance
        E = E + Log.E_useful(i) + Log.E_drag(i) + Log.E_rollres(i);
    end
    Log.E = 18*(1/0.9)*(4/3)*E/3.6e6; % Energy in KWh for 28 lap endurance at 75% DoD
    scatter3(Log.X(:),Log.Y(:),Log.Vel(:)*(18/5),10,Log.Vel(:)*(18/5),'filled')
    
    function latg = LatGlimit(v,Car)
        LatForce = Car.Mu_lat*(Car.Mass*9.81 + 0.5*1.23*Car.FrntArea*Car.Cl*v^2);
        latg = LatForce/Car.Mass;
    end
    
    
    function [a,vmax, Torque] = accel(v,Car,curr_latacc)
        Torque = 102;
        Fdrag = 0.5*1.23*Car.Cd*Car.FrntArea*v^2;
        Frollres = 0.05*(9.81*Car.Mass + 0.5*1.23*Car.Cl*Car.FrntArea*v^2);
        a = Torque*Car.FDR/(Car.Rw*Car.Mass) - (Fdrag + Frollres)/Car.Mass;
        a = min(a, Car.accel_limit);
        vmax = (Car.AccuVoltage*Car.SpecificLoadSpeedConstant*Car.Rw*2*pi)/(60*Car.FDR);
        if a > (1-(curr_latacc/Car.LatG_limit)^2)*Car.accel_limit
            a = (1-(curr_latacc/Car.LatG_limit)^2)*Car.accel_limit;
            Torque = (Car.Mass*a + (Fdrag + Frollres))*Car.Rw/Car.FDR;
        end
    end
end