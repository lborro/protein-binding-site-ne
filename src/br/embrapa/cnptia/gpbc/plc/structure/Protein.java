package br.embrapa.cnptia.gpbc.plc.structure;

import java.io.File;

import br.embrapa.cnptia.cbi.sdl.core.PDBIO;
import br.embrapa.cnptia.cbi.sdl.core.Structure;

public class Protein {
	private Structure structure;
	private String proteinID;
	private ProteinLigandComplex complex;

	public Protein(File pdbFile) throws Exception {
		if(pdbFile == null) {
			throw new NullPointerException("PDB file can not be null!");
		}

		if(!pdbFile.exists()) {
			throw new java.io.FileNotFoundException("PDB File not found: " + pdbFile.getName());
		}
		structure = PDBIO.read(pdbFile, 0);
		proteinID = pdbFile.getName().substring(0, pdbFile.getName().length()-4);
		complex = null;
	}

	public Structure getStructure() {
		return structure;
	}

	public String getProteinID() {
		return proteinID;
	}

	public ProteinLigandComplex getProteinLigandComplex () {
		return complex;
	}

	public void setProteinLigandComplex(ProteinLigandComplex complex) {
		this.complex = complex;
	}
}
