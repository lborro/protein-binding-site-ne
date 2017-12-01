package br.embrapa.cnptia.gpbc.plc.structure;

import java.io.File;
import java.io.FileInputStream;
import java.util.HashMap;
import java.util.Map;

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.io.Mol2Reader;

import org.openscience.cdk.tools.manipulator.ChemFileManipulator;

import br.embrapa.cnptia.cbi.sdl.core.Atom;
import br.embrapa.cnptia.cbi.sdl.core.Bond;
import br.embrapa.cnptia.cbi.sdl.core.Chain;
import br.embrapa.cnptia.cbi.sdl.core.Hetatm;
import br.embrapa.cnptia.cbi.sdl.core.IAtom;
import br.embrapa.cnptia.cbi.sdl.core.Structure;
import br.embrapa.cnptia.gpbc.plc.utils.Utilities;

public class Ligand {
	private Hetatm hetGroup;
	private String ligandID;

	private double massCenter[] = {0, 0, 0};
	private double geoCenter[]  = {0, 0, 0};

	public Ligand(File mol2File) throws Exception {
		this(mol2File, 0);
	}

	public Ligand(File mol2File, int offset) throws Exception {
		if(mol2File == null) {
			throw new NullPointerException("MOL2 file can not be null!");
		}

		if(!mol2File.exists()) {
			throw new java.io.FileNotFoundException("MOL2 File not found: " + mol2File.getName());
		}

		Mol2Reader reader = new Mol2Reader(
				new FileInputStream(mol2File)
				);

		IAtomContainer molecule = null;
		ChemFile chemFile = (ChemFile) reader.read((ChemObject)new ChemFile());
		molecule  = ChemFileManipulator.getAllAtomContainers(chemFile).get(0);

		hetGroup = null;
		if (molecule != null) {
			hetGroup = convert(molecule, offset);
			ligandID = mol2File.getName().substring(0,mol2File.getName().length()-5);
		}
		reader.close();
		geoCenter = Utilities.calculateGeoCenter(hetGroup);
		massCenter = Utilities.calculateMassCenter(hetGroup);
	}

	public Hetatm getHetGroup() {
		return hetGroup;
	}

	public String getLigandID() {
		return ligandID;
	}

	public double[] getMassCenter() {
		return this.massCenter;
	}

	public double[] getGeoCenter() {
		return this.geoCenter;
	}

	private Hetatm convert(IAtomContainer molecule, int offset) {
		Hetatm hetAtm = new Hetatm();
		hetAtm.setName("MOL");
		hetAtm.setPDBName("MOL");
		hetAtm.setInsCode(' ');
		hetAtm.setNumber(1);

		Chain chain = new Chain();
		chain.setChainId(" ");
		chain.addResidue(hetAtm);

		Structure structure = new Structure();
		structure.setPdbId("");
		structure.addChain(chain);

		Map<org.openscience.cdk.interfaces.IAtom, Integer> atom2serial = new HashMap<>();

		for(int j = 0; j < molecule.getAtomCount(); j++) {
			org.openscience.cdk.interfaces.IAtom atom = molecule.getAtom(j);
			Atom newAtom = new Atom();
			newAtom.setAtomSerial(j + 1 + offset);
			atom2serial.put(atom, j+1 + offset);
			newAtom.setName(atom.getID().trim());

			double coord[] = new double[3]; 
			atom.getPoint3d().get(coord);
			newAtom.setCoords(coord);
			newAtom.setElement(Utilities.symbolToElement(atom.getSymbol().trim()));
			hetAtm.addAtom(newAtom);
		}

		for(IBond bond : molecule.bonds()) {
			IAtom a1 = hetAtm.getAtomBySerial(atom2serial.get(bond.getAtom(0)));
			IAtom a2 = hetAtm.getAtomBySerial(atom2serial.get(bond.getAtom(1)));
			new Bond(a1, a2, bond.getOrder().numeric(), true);
		}
		return hetAtm;
	}
}