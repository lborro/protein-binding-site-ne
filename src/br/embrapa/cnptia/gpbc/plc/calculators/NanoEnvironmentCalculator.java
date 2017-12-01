package br.embrapa.cnptia.gpbc.plc.calculators;
import java.io.File;
import java.util.concurrent.Callable;

import br.embrapa.cnptia.gpbc.plc.data.MaxContactsTable;
import br.embrapa.cnptia.gpbc.plc.descriptors.protein.ProteinNanoEnvironment;
import br.embrapa.cnptia.gpbc.plc.structure.ProteinLigandComplex;

public class NanoEnvironmentCalculator implements Callable<ProteinNanoEnvironment> {
	private final File protein, ligand;
	private final MaxContactsTable maxCon;
	private final double dist;

	public NanoEnvironmentCalculator (File pdbFile, File mol2File, MaxContactsTable maxCon, double dist) throws Exception {
		if (pdbFile == null) 
			throw new NullPointerException("PDB file describing the protein cannot be null !");
		if (mol2File == null) 
			throw new NullPointerException("Mol2 file describing the ligand cannot be null !");
		if (maxCon == null) 
			throw new NullPointerException("Table of maximum contacts per AA cannot be null!");
		if(dist <= 0 || dist > 12)
			throw new Exception("Distance must be greater than zero and less or equal than 12.0");

		this.protein = pdbFile;
		this.ligand = mol2File;
		this.maxCon = maxCon;
		this.dist = dist;
	}

	@Override
	public ProteinNanoEnvironment call() throws Exception {
		try{
			ProteinLigandComplex complex = new ProteinLigandComplex(protein, ligand, 12);
			ProteinNanoEnvironment nano = new ProteinNanoEnvironment(complex, maxCon);
			nano.calculate(dist);
			return nano;
		} catch (Exception e) {
			return null;
		}
	}
}
