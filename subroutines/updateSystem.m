
function [mySystem,currentSummary] = updateSystem(mySystem,state,i,cnst)

    % --- copy all variables back to mySystem ---
    % tumor cell coordinates and properties ------------------------------
    
    mySystem.TU.killed=state.TUcells.killed;
    mySystem.TU.TUcells = state.TUcells.position;         % list of tumor cells
    mySystem.TU.TUprop.isStem = state.TUcells.isStem;      % tumor cell stemness
    mySystem.TU.TUprop.Pcap = state.TUcells.Pcap;          % remaining cell divs
    mySystem.TU.TUprop.Antigen = state.TUcells.Antigen;    % cell antigenicity
    mySystem.TU.TUprop.damage = state.TUcells.damage;    % cell antigenicity
    mySystem.TU.TUprop.glyc_rate = state.TUcells.glyc_rate; %rate of glucose uptake
    mySystem.TU.TUprop.isAcidResistant = state.TUcells.isAcidResistant; %indicator of acid resistant cells
    mySystem.TU.TUprop.PDL = state.TUcells.PDL; %indicator of PDL1
    % lymphocyte coordinates and properties ------------------------------
    mySystem.IM.IMcells = state.Lymphocytes.position;                  % list of lymphocytes
    mySystem.IM.IMprop.Kcap = state.Lymphocytes.Kcap;          % remaining kills
    mySystem.IM.IMprop.Pcap = state.Lymphocytes.Pcap;          % remaining cell divs
    mySystem.IM.IMprop.engaged = state.Lymphocytes.engaged;    % lymphocyte engagement
    mySystem.IM.IMprop.quiescent = state.Lymphocytes.quiescent;

    % global system maps -------------------------------------------------
    mySystem.grid.ChtaxMap = state.env.ChtaxMap;              % chemotaxis map
    mySystem.grid.HypoxMap = state.env.NecroMap;              % hypoxia map
    mySystem.grid.IFNgMap = state.env.IFNgMap;                % IFNgamma map
    mySystem.grid.GlucMap = state.env.GlucMap;                % glucose map
    mySystem.grid.ProtMap = state.env.ProtMap;                % Proton Map
    mySystem.grid.ATPMap = max(0,state.env.GlucMap*mySystem.params.glucosePointConsumption-1/5*state.env.NecroMap*mySystem.params.oxygenPointConsumption)*2+state.env.NecroMap*29/5*mySystem.params.oxygenPointConsumption;
    % necrosis and fibrosis ----------------------------------------------
    mySystem.grid.Ln = state.env.Ln;                          % necrosis map
    mySystem.grid.Lf = state.env.Lf;                          % fibrosis (stroma) map
    mySystem.grid.Lh = state.env.Lh;                          % hypoxia map
    % internal grids for faster computation of occupancy -----------------
    mySystem.grid.Lt = false(size(state.env.Ln));
    mySystem.grid.Lt(mySystem.TU.TUcells) = true;                          % tumor cell map
    mySystem.grid.StepsDone = i;
    
    % create immune grid
    mySystem.grid.Li = false(size(state.env.Ln));
    mySystem.grid.Li(mySystem.IM.IMcells) = true;
    
    currentSummary = summarizeSystem_3D(mySystem,cnst);
    
    
end
