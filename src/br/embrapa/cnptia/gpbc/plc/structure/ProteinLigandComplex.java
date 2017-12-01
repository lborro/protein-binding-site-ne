package br.embrapa.cnptia.gpbc.plc.structure;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.Chain;
import br.embrapa.cnptia.cbi.sdl.core.IAtom;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;

public class ProteinLigandComplex {
	private Protein protein;
	private Ligand ligand;
	private String complexID;
	private LigandNeighborhood ligandNeighborhood;


	public ProteinLigandComplex(File pdbFile, File mol2File, double maxDist) throws Exception {
		this.protein = new Protein(pdbFile);
		int offset = protein.getStructure().getAtoms(true, false).length + 1;
		this.ligand = new Ligand(mol2File, offset);
		this.complexID = this.protein.getStructure().getPdbId() + ":" + this.ligand.getLigandID();
		this.ligandNeighborhood = new LigandNeighborhood(this, maxDist);
		this.protein.setProteinLigandComplex(this);
	}

	public Protein getProtein() {
		return protein;
	}

	protected Ligand getLigand() {
		return ligand;
	}


	public String getID() {
		return complexID;
	}

	public List<IAtom> getLigandNeighborAtoms() {
		return new ArrayList<>(ligandNeighborhood.getNeighborAtoms());
	}

	public List<IResidue> getLigandNeighborResidues() {
		return new ArrayList<>(ligandNeighborhood.getNeighborResidues());
	}

	public List<Chain> getLigandNeighborChains() {
		return new ArrayList<>(ligandNeighborhood.getNeighborChains());
	}

	public Double getDistanceFromAAToLigand(AminoAcid aa) {
		return ligandNeighborhood.getDistanceToLigand(aa);
	}

	public List<ProteinLigandAtomPair> getProteinLigandAtomPairs() {
		return new ArrayList<>(ligandNeighborhood.getProteinLigandAtomPairs());
	}
}