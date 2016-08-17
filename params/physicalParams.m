%some physical params (heat source etc.)

%Assign heat source field
physical.heatSourceField = zeros(dFine.nEl, 1);
%Force contributions due to heat flux and source
physical.fs = get_heat_source(physical.heatSourceField, dFine);
physical.fh = get_flux_force(dFine, physical);