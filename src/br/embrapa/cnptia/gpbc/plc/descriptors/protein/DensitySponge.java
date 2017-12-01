package br.embrapa.cnptia.gpbc.plc.descriptors.protein;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.Chain;
import br.embrapa.cnptia.cbi.sdl.core.IAtom;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.cbi.sdl.utils.tree.KdTree;
import br.embrapa.cnptia.gpbc.plc.structure.Protein;
import br.embrapa.cnptia.gpbc.plc.utils.Utilities;

public class DensitySponge extends AbstractProteinDescriptor {

	private static final double[] radii = {3.0, 4.0, 5.0, 6.0, 7.0};

	public DensitySponge(Protein protein) throws Exception {
		super(protein, new String[]{
				"Density3",
				"Density4", 
				"Density5",
				"Density6",
				"Density7", 
				"Sponge3", 
				"Sponge4", 
				"Sponge5", 
				"Sponge6", 
				"Sponge7"
		});

		calculate();
	}

	private void calculate() throws Exception {
		KdTree<IAtom> kdTree = new KdTree<>(3);
		for(Chain chain : this.getProtein().getStructure().getChains())
			for(IResidue residue : chain.getResidues()){
				if(!Utils.isAminoAcid(residue.getName())) continue;
				for(IAtom atom : residue.getAtoms()) 
					if(atom.getElement().isHeavyAtom()) kdTree.addPoint(atom.getCoords(), atom);
			}

		for(IResidue residue : this.getProtein().getProteinLigandComplex().getLigandNeighborResidues()) {
			if(!Utils.isAminoAcid(residue.getName())) continue;
			AminoAcid aa = (AminoAcid) residue;

			double center[];
			try{
				center = aa.getLHA().getCoords();
			} catch(Exception e) {
				continue;
			}

			List<IAtom> nearestAtoms = kdTree.findNearestNeighbors(center, radii[radii.length-1] + 2.0).getAll();
			Double values[] = new Double[2*radii.length];
			for(int r = 0; r < radii.length; r++) {
				double radius = radii[r];
				double totalMass = 0.0, totalVolume = 0.0;
				Map<IAtom, Double> intersectAtomsVolumes = new HashMap<>();
				for(IAtom atom : nearestAtoms){
					double dist = Utils.euclideanSquaredDistance(center, atom.getCoords());
					if(dist < (radius + atom.getElement().getVDWRadius()) * (radius + atom.getElement().getVDWRadius())){
						final double vdwRadii = atom.getElement().getVDWRadius();
						double atomVolume = 4.0/3.0 * Math.PI * vdwRadii * vdwRadii * vdwRadii;
						double volume = computeIntersectionVolume(atomVolume, radius, vdwRadii, dist);
						intersectAtomsVolumes.put(atom, volume);
						totalVolume += calculateAtomsIntersection(intersectAtomsVolumes, atom, volume, radius, center);

						double volumeRatio = volume / atomVolume;
						totalMass += atom.getElement().getAtomicMass() * volumeRatio;
					}
				}
				double sphereVolume = 4.0/3.0 * Math.PI * radius * radius * radius;
				//Density
				values[r] = Utils.round(totalMass / sphereVolume, 4);
				//Sponge
				values[radii.length + r] = Utils.round((sphereVolume - totalVolume)/sphereVolume, 4);
			}
			this.addDescriptorValues(aa, values);
		}
	}

	private double calculateAtomsIntersection(
			Map<IAtom, Double> atoms, IAtom lastAtom,
			double volume, double sphereRadius, double[] centerXYZ) 
	{

		double finalVolume = volume;
		for(IAtom atom : atoms.keySet()){
			if (!atom.equals(lastAtom) && (
					Utilities.isWithinSphere(lastAtom, centerXYZ, sphereRadius) || 
					Utilities.isWithinSphere(atom, centerXYZ, sphereRadius)
					))
			{
				double dist = Utils.euclideanSquaredDistance(lastAtom.getCoords(), atom.getCoords());									
				double atomVolume = atoms.get(atom);
				double intersectVolume;
				if(lastAtom.getElement().getVDWRadius() > atom.getElement().getVDWRadius())
					intersectVolume = computeIntersectionVolume(
							atomVolume, 
							lastAtom.getElement().getVDWRadius(), 
							atom.getElement().getVDWRadius(), 
							dist);
				else
					intersectVolume = computeIntersectionVolume(
							volume, 
							atom.getElement().getVDWRadius(),
							lastAtom.getElement().getVDWRadius(), 
							dist);

				if(intersectVolume > 0) finalVolume -= intersectVolume;
			}
		}
		return finalVolume;
	}

	private double computeIntersectionVolume(double totalVolume, double R, double r, double dd) {
		if( dd > (R+r)*(R+r) )
			return 0.0;
		// Atom is totally inside the Sphere
		if (dd < (R-r)*(R-r))
			return totalVolume;

		final double d = Math.sqrt(dd);
		double volume = Math.PI * (R + r - d)*(R + r -d)*(d*d + 2*d*r - 3*r*r + 2*d*R + 6*r*R - 3*R*R)/(12.0*d);
		return volume;				
	}		
}
