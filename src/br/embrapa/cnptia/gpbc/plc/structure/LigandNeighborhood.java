package br.embrapa.cnptia.gpbc.plc.structure;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.Chain;
import br.embrapa.cnptia.cbi.sdl.core.IAtom;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.cbi.sdl.utils.tree.KdTree;

class LigandNeighborhood {
	private ProteinLigandComplex complex;
	private Set<IAtom> neighborAtoms;
	private Set<IResidue> neighborResidues;
	private Set<Chain> neighborChains;

	private Set<ProteinLigandAtomPair> pairs;

	private Map<AminoAcid, Double> distanceToLigand = new HashMap<>();

	public LigandNeighborhood(ProteinLigandComplex complex, double maxDistance) throws NullPointerException {
		if(complex == null) {
			throw new NullPointerException("Protein Ligand Complex can not be null!");
		}
		this.complex = complex;

		this.neighborChains = new HashSet<>();
		this.neighborResidues = new HashSet<>();
		this.neighborAtoms = new HashSet<>();
		this.pairs = new HashSet<>();

		defineBindingSite(maxDistance);
	}

	private void defineBindingSite(double maxDistance) {
		KdTree<IAtom> kdTree =  new KdTree<>(3);
		for(Chain chain : complex.getProtein().getStructure())
			for(IResidue residue : chain)
				for(IAtom atom : residue){ 
					//only heavy atoms
					if(atom.getElement().isHydrogen()) continue;
					kdTree.addPoint(atom.getCoords(), atom);
				}

		for(IAtom atom : complex.getLigand().getHetGroup()){
			if(atom.getElement().isHydrogen()) continue;
			List<IAtom> atomList = kdTree.findNearestNeighbors(atom.getCoords(), maxDistance).getAll();
			for(IAtom neighbor : atomList) {
				double dist = Utils.euclideanDistance(atom.getCoords(), neighbor.getCoords());
				if(dist <= maxDistance) {
					neighborAtoms.add(neighbor);
					neighborResidues.add(neighbor.getResidue());
					neighborChains.add(neighbor.getResidue().getChain());
					pairs.add(new ProteinLigandAtomPair(atom, neighbor));

					if(Utils.isAminoAcid(neighbor.getResidue().getName())) {
						AminoAcid aa = (AminoAcid) neighbor.getResidue();
						Double minDist = distanceToLigand.get(aa);
						if(minDist == null || dist < minDist) {
							distanceToLigand.put(aa, new Double(dist));
						}
					}
				}
			}
		}
	}

	public ProteinLigandComplex getProteinLigandComplex() {
		return this.complex;
	}

	public Set<IAtom> getNeighborAtoms() {
		return this.neighborAtoms;
	}

	public Set<IResidue> getNeighborResidues() {
		return this.neighborResidues;
	}

	public Set<Chain> getNeighborChains() {
		return neighborChains;
	}

	public Double getDistanceToLigand(AminoAcid aa) {
		Double value = distanceToLigand.get(aa);
		if(value != null)
			return this.distanceToLigand.get(aa);
		return 0d;
	}

	public Set<ProteinLigandAtomPair> getProteinLigandAtomPairs() {
		return this.pairs;
	}
}
