%J Grajek, 2022

function mySystem = initializeFromBall(mySystem, r)
%Initialize cell locations (sphere in center of domain with radius r)
if mySystem.is3D
    ballIdx = FindBallIdx3D(mySystem, r);
else
    ballIdx = FindBallIdx(mySystem, r);
end

Nseeded = length(ballIdx);

%Initialize cells with default values
mySystem.TU.TUcells = int32(ballIdx(:)');
mySystem.TU.TUprop.isStem = true(1,Nseeded);
mySystem.TU.TUprop.Pcap = uint8(mySystem.params.TUpmax*ones(1,Nseeded));
mySystem.TU.TUprop.Antigen = single(0.3*ones(1,Nseeded)); %we want lymphocytes to be able to attack the initial cells. After inocculation they would have acuires some mutations
mySystem.TU.TUprop.damage = uint8(0*ones(1, Nseeded));
mySystem.TU.TUprop.glyc_rate = ones(1, Nseeded);
mySystem.TU.TUprop.PDL = logical(zeros(1, Nseeded));
mySystem.TU.TUprop.isAcidResistant = logical(zeros(1, Nseeded));

pos_CA9 = randsample(length(ballIdx), round(length(ballIdx)*mySystem.params.CA9freq));
mySystem.TU.TUprop.isAcidResistant(pos_CA9) = 1;

pos_PDL1 = randsample(length(ballIdx), round(length(ballIdx)*mySystem.params.PDL1freq));
mySystem.TU.TUprop.PDL(pos_PDL1) = 1;

%% HELPER FUNCTIONS
%FUNCTION FINDING INDICES OF BALL WITH CENTER AT PLANE CENTER AND
%RADIUS R in 2D

    function ballIdx = FindBallIdx(mySystem, r)
        ballIdx = [];
        N1 = mySystem.grid.N;
        N2 = mySystem.grid.M;
        
        planeCenter = floor(N1/ 2) + floor((N2 / 2))*N1;
        [planeCenterX, planeCenterY] = ind2sub([N1, N2],planeCenter);
        
        for k=planeCenter-r*N1-r:planeCenter+r*N1+r
            [x,y] = ind2sub([N1, N2], k);
            if (x-planeCenterX)^2+(y-planeCenterY)^2<=r^2
                ballIdx(end+1)=k;
            end
        end
    end

%FUNCTION FINDING INDICES OF BALL WITH CENTER AT PLANE CENTER AND
%RADIUS R in 3D

    function ballIdx = FindBallIdx3D(mySystem, r)
        ballIdx = [];
        N1 = mySystem.grid.N;
        N2 = mySystem.grid.M;
        N3 = mySystem.grid.P;
        
        planeCenter = floor(N1/ 2) + floor((N2 / 2))*N1+floor(N3/2)*N2*N1+1;
        [planeCenterX, planeCenterY, planeCenterZ] = ind2sub([N1, N2, N3],planeCenter);
        
        for k=planeCenter-r*N1-r-r*N1*N2:planeCenter+r*N1+r+r*N1*N2
            [x,y, z] = ind2sub([N1, N2, N3], k);
            if (x-planeCenterX)^2+(y-planeCenterY)^2 +(z-planeCenterZ)^2<=r^2
                ballIdx(end+1)=k;
            end
        end
    end

end