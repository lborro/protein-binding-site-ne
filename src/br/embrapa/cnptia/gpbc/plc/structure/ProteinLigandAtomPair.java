package br.embrapa.cnptia.gpbc.plc.structure;

import br.embrapa.cnptia.cbi.sdl.core.IAtom;

public class ProteinLigandAtomPair implements Comparable<ProteinLigandAtomPair> {
	private IAtom ligandAtom, proteinAtom;
	private double distance;

	public ProteinLigandAtomPair(IAtom ligandAtom, IAtom proteinAtom) {
		this.ligandAtom = ligandAtom;
		this.proteinAtom = proteinAtom;

		distance = 0;
		for(int i = 0; i < 3; i++) distance += Math.pow(ligandAtom.getCoords()[i] - proteinAtom.getCoords()[i], 2);
		distance = Math.sqrt(distance);
	}

	@Override
	public int compareTo(ProteinLigandAtomPair other) {
		int r = ligandAtom.compareTo(other.ligandAtom);
		if(r == 0) return proteinAtom.compareTo(other.proteinAtom);
		return r;
	}

	public double getDistance() {
		return distance;
	}

	public IAtom getLigandAtom() {
		return this.ligandAtom;
	}

	public IAtom getProteinAtom() {
		return this.proteinAtom;
	}

}
