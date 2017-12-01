package br.embrapa.cnptia.gpbc.plc.descriptors.protein;

import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.IAtom;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.descriptors.contacts.Contact;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.cbi.sdl.utils.tree.KdTree;
import br.embrapa.cnptia.gpbc.plc.structure.Protein;

public class EnergyDensity extends AbstractProteinDescriptor {

	private static final double[] radii = {3.0, 4.0, 5.0, 6.0, 7.0};

	private ContactsPerResidue conRes;

	protected EnergyDensity(Protein protein, ContactsPerResidue conRes) throws Exception {	
		super(protein, new String[]{
				"EnergyDensity3",
				"EnergyDensity4", 
				"EnergyDensity5",
				"EnergyDensity6",
				"EnergyDensity7"
		}); 

		if(conRes == null) throw new NullPointerException("ContactsPerResidue cannot be null!");
		this.conRes = conRes;

		calculate();
	}

	private void calculate() throws Exception {
		KdTree<Set<Contact>> tree =  new KdTree<>(3);
		for(IResidue residue : this.getProtein().getProteinLigandComplex().getLigandNeighborResidues()) {
			if(!Utils.isAminoAcid(residue.getName())) continue;
			AminoAcid aa = (AminoAcid) residue;
			for(IAtom atom : aa.getAtoms()) {
				Set<Contact> set = conRes.getContacts(atom);
				if(set != null) tree.addPoint(atom.getCoords(), set);
			}
		}

		for(IResidue residue : this.getProtein().getProteinLigandComplex().getLigandNeighborResidues()) {
			if(!Utils.isAminoAcid(residue.getName())) continue;
			AminoAcid aa = (AminoAcid) residue;
			IAtom center = aa.getLHA();
			if(center == null) continue;

			List<Set<Contact>> nearestContacts = tree.findNearestNeighbors(center.getCoords(), radii[radii.length-1]).getAll();
			Set<Contact> allContacts = new HashSet<Contact>();
			for(Set<Contact> set : nearestContacts)
				allContacts.addAll(set);

			Double[] energy = calculateEnergyForAtom(center, allContacts);
			Double[] ced = new Double[radii.length];			
			Arrays.fill(ced, 0.0);	
			// Volume = 4/3*pi*R^3;
			final double volume = 4.0/3.0 * Math.PI; //constant
			for(int r = 0; r < radii.length; r++)
				ced[r] = Utils.round(energy[r] / (volume * radii[r] * radii[r] * radii[r]), 4);

			this.addDescriptorValues(aa, ced);
		}

	}

	private Double[] calculateEnergyForAtom(IAtom centralAtom, Set<Contact> nearestContacts) {
		// Contact energy 
		Double[] energy = new Double[radii.length];
		Arrays.fill(energy, 0d);

		// Iterate over all protein contacts and only considers those
		// within the sphere
		for(Contact contact : nearestContacts){
			// Distance form central atom and the first atom of this contact
			double dist1 = Utils.euclideanDistance(centralAtom.getCoords(), contact.getAtom1().getCoords());
			// Distance form central atom and the second atom of this contact
			double dist2 = Utils.euclideanDistance(centralAtom.getCoords(), contact.getAtom2().getCoords());

			for(int r = 0; r < radii.length; r++){
				// If any of the atoms is inside the sphere of radii[i] sum the energy of this contact
				if (dist1 <= radii[r] || dist2 <= radii[r]){
					if(contact.getDistance() <= 4.0)
						energy[r] += contact.getContactType().getEnergy();
					else if(contact.getDistance() <= 6.0)
						energy[r] += contact.getContactType().getEnergy() * 0.5;
					else if(contact.getDistance() <= 9.0)
						energy[r] += contact.getContactType().getEnergy() * 0.2;
					else 
						energy[r] += contact.getContactType().getEnergy() * 0.01;
				}
			}
		}

		return energy;
	}
}
