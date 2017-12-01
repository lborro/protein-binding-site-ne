package br.embrapa.cnptia.gpbc.plc.descriptors.protein;

import java.io.File;

import br.embrapa.cnptia.gpbc.plc.structure.Protein;

public class ElectrostaticPotential extends AbstractProteinDescriptor {
	protected ElectrostaticPotential(Protein protein, String epDir) throws Exception {
		super(protein, new String [] {"EP_Surface"});
	}

	@SuppressWarnings("unused")
	private void readEP(File epFile) {
	}
}
