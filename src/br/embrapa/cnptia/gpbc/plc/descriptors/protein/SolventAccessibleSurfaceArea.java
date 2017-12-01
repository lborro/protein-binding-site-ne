package br.embrapa.cnptia.gpbc.plc.descriptors.protein;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.io.FileUtils;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.Chain;
import br.embrapa.cnptia.cbi.sdl.core.IAtom;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.descriptors.hydrophobicity.HydrophobicityScale;
import br.embrapa.cnptia.cbi.sdl.descriptors.hydrophobicity.KyteDoolitteHydrophobicityScale;
import br.embrapa.cnptia.cbi.sdl.utils.Constants;
import br.embrapa.cnptia.cbi.sdl.utils.Grep;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.gpbc.plc.structure.Protein;

public class SolventAccessibleSurfaceArea extends AbstractProteinDescriptor {
	private File naccessPath;
	private File radiiFile;
	private double probeRadius;

	private File temporaryDir;

	public SolventAccessibleSurfaceArea(
			Protein protein,
			String naccessPath, 
			String radiiFile, 
			double probeRadius, 
			String tmpDir) throws Exception
	{
		super(protein, new String[]{"Accessibility", "Hydrophobibity"});

		if (naccessPath == null)
			throw new NullPointerException("naccess path cannot be null!");

		if (radiiFile == null)
			throw new NullPointerException("Atoms radii file cannot be null.");

		this.naccessPath = new File(naccessPath); 
		if (!this.naccessPath.exists()) 
			throw new FileNotFoundException("naccess path does not exist: " + naccessPath);				

		this.radiiFile = new File(radiiFile);
		if (!this.radiiFile.exists()) 
			throw new FileNotFoundException("Atoms radiues file does not exist: " + radiiFile);

		this.probeRadius = probeRadius;

		temporaryDir = new File(tmpDir + File.separator + System.nanoTime());
		if (!temporaryDir.mkdirs()){
			throw new FileNotFoundException("Could not create temporary working directory: " + temporaryDir.getAbsolutePath());
		}

		calculate();
	}

	private void calculate() throws Exception{
		final String pdbCode = this.getProtein().getProteinID();
		File pdbFile = new File(this.getProtein().getStructure().getPdbFilename());
		File proteinFile = null;
		try{
			List<String> chains = Grep.split(pdbFile, Constants.ATOM_PATTERN, Constants.MODEL_PATTERN, Constants.TER_PATTERN);	
			StringBuilder strBuilder = new StringBuilder();
			for(String chStr : chains)
				if (chStr.length() > 0)	strBuilder.append(chStr);

			proteinFile = new File(temporaryDir + File.separator + pdbCode + "_protein.pdb");
			BufferedWriter bw = new BufferedWriter(new FileWriter(proteinFile));
			bw.write(strBuilder.toString()); 
			bw.close();
		} catch(Exception e) {
			throw new Exception("Surface Accessible Surface Area: error creating temporary files");
		}

		Map<IAtom, Double> atomSasaMap = null;	
		try{
			executeNaccess(temporaryDir, proteinFile);
			String fileName = proteinFile.getName();
			File sasaFile =  new File(temporaryDir + File.separator + fileName.substring(0, fileName.length()-4) + ".asa");

			atomSasaMap = readProteinSASAFile(sasaFile);
		} catch(Exception e) {		
			throw new Exception("Surface Accessible Surface Area: error reading SASA files");
		}  finally{
			FileUtils.deleteDirectory(temporaryDir);
		}


		HydrophobicityScale hydroScale = new KyteDoolitteHydrophobicityScale();

		for(Chain chain: this.getProtein().getStructure()) {
			for(IResidue residue : chain){				
				if(!Utils.isAminoAcid(residue.getName())) continue;
				AminoAcid aa = (AminoAcid) residue;
				double sasa = 0;
				for(IAtom atom : aa) {
					Double value = atomSasaMap.get(atom);
					if(value != null)  sasa += value;
				}
				Double hydro = sasa*hydroScale.getHydrophobicity(aa.getAminoAcidType())/aa.getAminoAcidType().getMaxAccessibility();	
				this.addDescriptorValues(aa, new Double[]{sasa, hydro});
			}
		}
	}

	private void executeNaccess(File temporaryDir, File pdbFile) throws Exception {
		File tmpNaccessPath = new File(temporaryDir + File.separator +naccessPath.getName());
		File tmpRadiiFile = new File(temporaryDir + File.separator +radiiFile.getName());		

		Files.copy(naccessPath.toPath(), tmpNaccessPath.toPath(), 
				StandardCopyOption.COPY_ATTRIBUTES, StandardCopyOption.REPLACE_EXISTING);

		Files.copy(radiiFile.toPath(), tmpRadiiFile.toPath(),
				StandardCopyOption.COPY_ATTRIBUTES, StandardCopyOption.REPLACE_EXISTING);	

		String[] cmdArray = new String[]{"./"+tmpNaccessPath.getName(), pdbFile.getName(), 
				"-p", ""+probeRadius,"-r",tmpRadiiFile.getName(),"-h" ,"-a"}; 				
		int exitCode;

		try {
			exitCode = Utils.callWithTimeout(cmdArray, temporaryDir, null, 1);
		} catch (Exception e) {					
			throw e;
		}

		if (exitCode != 0)			
			throw new Exception(
					String.format("Error to calculate solvent accessible surface for %s. Program return with code: %d", 
							pdbFile.getAbsolutePath(), 
							exitCode
							));

	}

	private Map<IAtom, Double> readProteinSASAFile(File SASAFile) throws Exception {
		if (!SASAFile.exists()) throw new FileNotFoundException("Solvent accessible surface file not found: " + SASAFile);
		Map<Integer,IAtom> atomMap = this.getProtein().getStructure().getAllAtomMap();
		Map<IAtom, Double> sasaMap = new HashMap<>();

		try (BufferedReader br = new BufferedReader(new FileReader(SASAFile)); ) {
			String line = null;					
			line = br.readLine();

			//if line is null already in the first time, the buffered reader
			// has a problem
			if (line == null) {
				br.close();
				throw new IOException("Could not read file, BufferedReader returns null!");
			}
			while(line != null){
				// ignore empty lines
				if (line.equals("") || (line.equals(System.getProperty("line.separator")))){
					line = br.readLine();
					continue;
				}

				int atomSerial = Integer.parseInt(line.substring(6, 11).trim());
				double acc = Double.parseDouble(line.substring(54,62).trim());
				IAtom atom = atomMap.get(atomSerial);
				if(atom != null) sasaMap.put(atom, new Double(acc));
				line = br.readLine();
			}
		}
		return sasaMap;
	}
}
