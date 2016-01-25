function [dcell, icell, pcell, pcell_split, reconstructions, reconstructions_raw, gegfmat, gegsmat, gegtmat, gegrmat, filters, params, exemplar_hists] = makefblockcolorize(colorizetrain, params)

[dcell, icell, exemplar_hists] = readdatacolorize(colorizetrain, params);

[pcell, pcell_split, gegfmat, gegsmat, gegtmat, gegrmat, reconstructions, reconstructions_raw, filters, params] = makeretfeatsgrawinscolorize(icell, dcell, params);

params.initial_model = @initialmodel_colorize;

end
