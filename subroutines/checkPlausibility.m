% JN Kather, NCT Heidelberg, 2017
% updated J Grajek, 2022

function plausible = checkPlausibility(mySystem)
plausible = true;

if sum([mySystem.params.TUpprol,mySystem.params.TUpmig,mySystem.params.TUpdeath])>1
    warning('Tumor parameter error');
    plausible = false;
end

if sum([mySystem.params.IMpprol,mySystem.params.IMpmig,mySystem.params.IMpdeath])>1
   warning('Immune cell (lymphocyte) parameter error');
       plausible = false;
end

end