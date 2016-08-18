function physical = physicalParams(domain, physical)
%add some physical params (heat source etc.) to physical struct

%Assign heat source field
physical.heatSourceField = zeros(domain.nEl, 1);
%Force contributions due to heat flux and source
physical.fs = get_heat_source(physical.heatSourceField, domain);
physical.fh = get_flux_force(domain, physical);
end