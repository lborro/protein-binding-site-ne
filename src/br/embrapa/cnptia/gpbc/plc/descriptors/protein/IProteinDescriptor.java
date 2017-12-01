package br.embrapa.cnptia.gpbc.plc.descriptors.protein;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;

public interface IProteinDescriptor {
	public String[] getDescriptorNames();
	public Double[] getDescriptorValues(AminoAcid aminoAcid);
}