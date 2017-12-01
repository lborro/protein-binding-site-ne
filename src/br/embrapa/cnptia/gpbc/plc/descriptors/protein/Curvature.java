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
import java.util.Map;
import java.util.concurrent.ExecutionException;

import org.apache.commons.io.FileUtils;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.Chain;
import br.embrapa.cnptia.cbi.sdl.core.FileConvert;
import br.embrapa.cnptia.cbi.sdl.core.IAtom;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.core.Structure;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.gpbc.plc.structure.Protein;

public class Curvature extends AbstractProteinDescriptor {

	public static final int RICHARDS_VDW_RADII_SET = 1;

	public static final int CHOTIA_VDW_RADII_SET = 2;

	private static final String RADII_FILENAME = "radii.txt";

	private File surfracePath;

	private File radiiFile;

	private int radiiSet;

	private double probeRadius;

	private String tmpDir;

	public Curvature(Protein protein, 
			String surfracePath, 
			String radiiFile, 
			int radiiSet, 
			double probeRadius,
			String tmpDir) throws Exception 
	{
		super(protein, new String [] {"curvature"});

		if (surfracePath == null)
			throw new NullPointerException("SURFACE RACER path cannot be null!");

		if (radiiFile == null)
			throw new NullPointerException("Atoms radii file cannot be null.");

		this.surfracePath = new File(surfracePath);
		if (!this.surfracePath.exists())
			throw new FileNotFoundException("SURFACE RACER path does not exist: " + surfracePath);				

		this.radiiFile = new File(radiiFile);
		if (!this.radiiFile.exists())
			throw new FileNotFoundException("Atoms radiues file does not exist: "+ radiiFile);

		if (radiiSet != RICHARDS_VDW_RADII_SET && radiiSet != CHOTIA_VDW_RADII_SET)
			throw new IllegalArgumentException("Unknown Radii set. Choose 1 for RICHARDS or 2 for CHOTIA radius set.");

		this.radiiSet = radiiSet;
		this.probeRadius = probeRadius;
		this.tmpDir = tmpDir;

		calculate();

	}

	private void calculate() throws Exception {
		File temporaryDir = new File(tmpDir + File.separator + System.nanoTime());
		if (!temporaryDir.mkdirs()){
			throw new FileNotFoundException("Could not create temporary working directory: " + temporaryDir.getAbsolutePath());
		}

		final String pdbCode = this.getProtein().getProteinID();

		try {
			File pdbFile = new File(temporaryDir + File.separator + pdbCode + ".pdb");
			BufferedWriter bw = new BufferedWriter(new FileWriter(pdbFile));
			for(Chain chain : this.getProtein().getProteinLigandComplex().getLigandNeighborChains()) {
				for(IResidue residue : chain) {
					if(Utils.isAminoAcid(residue.getName())){
						for(IAtom atom : residue) {
							if(atom.getElement().isHeavyAtom()) bw.write(FileConvert.toPDB(atom));
						}
					}
				}
			}
			bw.close();

			executeSurfaceRacer(temporaryDir, pdbFile);
			File curvFile = new File(temporaryDir + File.separator + pdbCode + ".txt");
			Map<IAtom, double[]> curvMap = readCurvatureFile(this.getProtein().getStructure(), curvFile);

			for(Chain chain: this.getProtein().getStructure()) {
				for(IResidue residue : chain){				
					if(!Utils.isAminoAcid(residue.getName())) continue;
					AminoAcid aa = (AminoAcid) residue;

					Double curvature = 0d, surfArea = 0d;
					for(IAtom atom : aa) {
						double values[]  = curvMap.get(atom);
						if(values != null){
							surfArea += values[0];
							curvature += values[0]*values[1];
						}
					}
					if(surfArea < 0.001) 
						this.addDescriptorValues(aa, new Double[] {0d});
					else
						this.addDescriptorValues(aa, new Double[] {curvature/surfArea});
				}
			}

		} catch(Exception e) {
			throw new Exception("Error calculating amino acid residues curvature for the following structure: " + pdbCode);
		} finally{
			FileUtils.deleteDirectory(temporaryDir);
		}

	}

	private void executeSurfaceRacer(File temporaryDir, File pdbFile) throws Exception {

		File tmpSurfacePath = new File(temporaryDir + File.separator + surfracePath.getName());
		File tmpRadiiFile   = new File(temporaryDir + File.separator + RADII_FILENAME);		

		Files.copy(surfracePath.toPath(), tmpSurfacePath.toPath(), 
				StandardCopyOption.COPY_ATTRIBUTES, StandardCopyOption.REPLACE_EXISTING);

		Files.copy(radiiFile.toPath(), tmpRadiiFile.toPath(),
				StandardCopyOption.COPY_ATTRIBUTES, StandardCopyOption.REPLACE_EXISTING);		


		String[] cmdArray = new String[]{tmpSurfacePath.getAbsolutePath()};
		int exitCode;
		String[] stdInput = new String[4];
		stdInput[0] = ""+radiiSet;
		stdInput[1] = pdbFile.getName();
		stdInput[2] = ""+probeRadius;			
		stdInput[3] = "3\n\n";

		try{
			exitCode = Utils.callWithTimeout(cmdArray, stdInput, 1);
		} catch (Exception e) {					
			throw e;
		}	

		if (exitCode != 0)			
			throw new ExecutionException(
					String.format("Error to calculate Curvature for %s. Program returned with code: %d", 
							pdbFile.getAbsolutePath(), exitCode)
					, null);
	}

	private Map<IAtom, double[]> readCurvatureFile(Structure structure, File curvFile) throws Exception {
		if (structure == null)	throw new NullPointerException("Structure cannot be null.");
		if (!curvFile.exists()) throw new FileNotFoundException("Curvature file not found: " + curvFile);
		Map<Integer,IAtom> atomMap = structure.getAllAtomMap();

		Map<IAtom, double[]> curvMap = new HashMap<>();

		try(BufferedReader br = new BufferedReader(new FileReader(curvFile));) {
			String line = null;						
			line = br.readLine();

			if (line == null) {
				br.close();
				throw new IOException("Could not read curvature file:" + curvFile);
			}			

			while(line != null){
				// ignore empty lines
				if (line.equals("") || (line.equals(System.getProperty("line.separator")))){
					line = br.readLine();
					continue;
				}
				if (line.startsWith("CAVITY")) break;	
				int atomSerial = Integer.parseInt(line.substring(6,11).trim());
				double surfValue = Double.parseDouble(line.substring(70,76).trim());
				double curvValue = Double.parseDouble(line.substring(76).trim());
				IAtom atom = atomMap.get(atomSerial);
				if(atom != null) {
					double values[] = {surfValue, curvValue};
					curvMap.put(atom,values);
				}
				line = br.readLine();
			}
			br.close();
		}
		return curvMap;
	}

}
