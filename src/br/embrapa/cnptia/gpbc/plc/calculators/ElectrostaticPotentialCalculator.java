package br.embrapa.cnptia.gpbc.plc.calculators;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;

import br.embrapa.cnptia.cbi.sdl.core.Chain;
import br.embrapa.cnptia.cbi.sdl.core.FileConvert;
import br.embrapa.cnptia.cbi.sdl.core.IAtom;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.core.PDBIO;
import br.embrapa.cnptia.cbi.sdl.core.SSBond;
import br.embrapa.cnptia.cbi.sdl.core.Structure;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;

public class ElectrostaticPotentialCalculator implements Callable<Object>{

	private static final String CHARGE_POTENTIAL_PARAMS = 
			"!fort.10\n"
					+ "perfil=80\n"
					+ "scale=2.5\n"
					+ "in(pdb,file=\"%pdbCode\")\n"
					+ "in(siz,file=\"%radiiFile\")\n"
					+ "in(crg,file=\"%chargeFile\")\n"
					+ "in(frc,frm=SELF)\n"
					+ "out(frc,file=\"%pdbCode.frc\")\n"
					+ "out(scrg,frm=PDB,file=\"%pdbCode.chrg\")\n"
					+ "atompotdist=0.6\n"
					+ "site(a,q,p,atpo)\n"
					+ "indi=2.0\n"
					+ "exdi=80.0\n"
					+ "prbrad=1.4\n"
					+ "ionrad=2.0\n"
					+ "salt=0.0000\n"
					+ "bndcon=4\n"
					+ "linit=500\n"
					+ "maxc=.0002";

	private static final String POTENTIAL_PARAMS = 
			"!fort.10.1\n"
					+ "perfil=80\n"
					+ "scale=2.5\n"
					+ "in(pdb,file=\"%pdbCode\")\n"
					+ "in(siz,file=\"%radiiFile\")\n"
					+ "in(crg,file=\"%chargeFile\")\n"
					+ "in(frc,file=\"%pdbCode.chrg\")\n"
					+ "out(frc,file=\"%pdbCode.pot\")\n"
					+ "site(a,x,p)\n"
					+ "indi=2.0\n"
					+ "exdi=80.0\n"
					+ "prbrad=1.4\n"
					+ "ionrad=2.0\n"
					+ "salt=0.0000\n"
					+ "bndcon=4\n"
					+ "linit=500\n"
					+ "maxc=.0002";

	private final String delphiPath;

	private final String reducePath;

	private final String parseRadiiFilePath;

	private final String parseChargesFilePath;

	private final String reduceDictionary;

	private final String tmpDir;

	private final String outDir;

	private final Structure structure;


	public ElectrostaticPotentialCalculator(File pdbFile, String delphiPath, String reducePath, String parseRadiiFilePath, String parseChargesFilePath, 
			String dictionaryPath, String tmpDir, String outDir) throws Exception {

		this.structure = PDBIO.read(pdbFile, 0);
		this.delphiPath = delphiPath;
		this.reducePath = reducePath;
		this.reduceDictionary = dictionaryPath;
		this.parseChargesFilePath = parseChargesFilePath;
		this.parseRadiiFilePath = parseRadiiFilePath;
		this.tmpDir = tmpDir;
		this.outDir = outDir;

		//delphi
		if(delphiPath == null || !(new File(delphiPath).exists()))
			throw new FileNotFoundException("Delphi path is null or not exist: " + delphiPath);
		//reduce
		if(reducePath == null || !(new File(reducePath).exists()))
			throw new FileNotFoundException("Reduce path is null or not exist: " + reducePath);			
		//parse radii file
		if(parseRadiiFilePath == null || !(new File(parseRadiiFilePath).exists()))
			throw new FileNotFoundException("Delphi radii file path is null or not exist: " + parseRadiiFilePath);
		//parse charges file
		if(parseChargesFilePath == null || !(new File(parseChargesFilePath).exists()))
			throw new FileNotFoundException("Delphi charge file path is null or not exist: " + parseChargesFilePath);
		//reduce dictionary
		if (dictionaryPath == null || !(new File(dictionaryPath).exists()))
			throw new FileNotFoundException("Reduce dictionary file path is null or not exist: " + dictionaryPath);	
	}

	@Override
	public Object call() throws Exception {
		// TODO Auto-generated method stub
		File temporaryDir = new File(tmpDir + File.separator + "ep_" + System.nanoTime());		
		if (!temporaryDir.mkdirs()){
			throw new FileNotFoundException("Could not create temporary working directory: " + temporaryDir.getAbsolutePath());
		}

		File pdb = new File(structure.getPdbFilename());
		String pdbCode = pdb.getName().substring(0,pdb.getName().lastIndexOf('.'));

		try{
			File cleanPdb = new File(temporaryDir + File.separator + pdbCode + ".pdb");
			BufferedWriter bw = new BufferedWriter(new FileWriter(cleanPdb));

			for(Chain chain : structure) {
				for(IResidue residue : chain) {
					if(Utils.isAminoAcid(residue.getName())){
						for(IAtom atom : residue) {
							if(atom.getElement().isHeavyAtom()) bw.write(FileConvert.toPDB(atom));
						}
					}
				}
			}

			bw.close();

			File protonatedPdb = addHydrogens(cleanPdb, temporaryDir);
			preparePDB(structure, protonatedPdb);
			generateAEP(protonatedPdb);
			Map<IResidue, double[]> potMap = generateSEP(protonatedPdb);

			List<IResidue> residues = new ArrayList<>(potMap.keySet());
			Collections.sort(residues);

			File outFile = new File(outDir + File.separator + pdbCode + ".ep");
			bw = new BufferedWriter(new FileWriter(outFile));

			for(IResidue residue : residues) {
				double[] potential = potMap.get(residue);
				double value;
				if (potential != null && potential[1] > 0)			
					value = potential[0]/potential[1];					
				else
					value = 0d;
				bw.write(residue.getChainId()  + "\t");
				bw.write(residue.getName()     + "\t");
				bw.write(residue.getNumber()   + "\t"); 
				bw.write(residue.getInsCode()  + "\t");
				bw.write(Utils.round(value, 4) + "\n");
			}
			bw.close();
		}catch(Exception e){
			throw e;
		}finally{
			//	FileUtils.deleteQuietly(temporaryDir);	
		}
		return null;
	}

	private void generateAEP(File protonatedPdb) throws Exception {
		String pdbPath = protonatedPdb.getAbsolutePath();
		String paramFileContent = CHARGE_POTENTIAL_PARAMS.replaceAll("%pdbCode", pdbPath)
				.replaceAll("%chargeFile", new File(parseChargesFilePath).getAbsolutePath())
				.replaceAll("%radiiFile", new File(parseRadiiFilePath).getAbsolutePath());

		File paramFile = new File(pdbPath + ".fort10");   

		try(BufferedWriter bw = new BufferedWriter(new FileWriter(paramFile));){
			bw.write(paramFileContent);
		} catch (IOException e) {
			throw e;
		}

		String[] cmdArray = new String[]{new File(delphiPath).getAbsolutePath(),paramFile.getAbsolutePath()};
		int exitCode = -1;

		try{
			exitCode = Utils.callWithTimeout(cmdArray, 5L);
		}catch(Exception e){
			throw e;
		}

		if (exitCode != 0)                      
			throw new ExecutionException(
					String.format("Error calculating electrostatic potential: %s. Program returned code %d", 
							protonatedPdb.getAbsolutePath(), exitCode)
					, null);

	}

	private Map<IResidue, double[]> generateSEP(File protonatedPdb) throws Exception {
		String pdbPath = protonatedPdb.getAbsolutePath();
		String paramFileContent = POTENTIAL_PARAMS.replaceAll("%pdbCode", pdbPath)
				.replaceAll("%chargeFile", new File(parseChargesFilePath).getAbsolutePath())
				.replaceAll("%radiiFile", new File(parseRadiiFilePath).getAbsolutePath());
		File paramFile = new File(protonatedPdb.getAbsolutePath()+".fort101");

		try(BufferedWriter bw = new BufferedWriter(new FileWriter(paramFile));){
			bw.write(paramFileContent);
		}catch (IOException e) {
			return null;
		}

		String[] cmdArray = new String[]{new File(delphiPath).getAbsolutePath(), paramFile.getAbsolutePath()};
		int exitCode = -1;

		try{
			exitCode = Utils.callWithTimeout(cmdArray, 5L);
		}catch(Exception e){
			throw e;
		}

		if(exitCode != 0)
			throw new ExecutionException(
					String.format("Error calculating surface electrostatic potential for %s. Program returned code %d",
							protonatedPdb.getAbsolutePath(), exitCode)
					, null);

		Structure protStructure = PDBIO.read(protonatedPdb, 0);
		Map<IResidue, double[]> potMap = new HashMap<>();

		Map<Integer,IAtom> atomMap = protStructure.getAllAtomMap();

		File potFile = new File(protonatedPdb.getAbsolutePath() + ".pot");
		try(BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(potFile)))) {
			String line = null;
			line = br.readLine();
			if (line == null) {
				IOException e = new IOException("Could not read file, BufferedReader returns null!");			
				throw e; 
			}			

			while(!line.matches("^ATOM.*")) 
				line = br.readLine();
			line = br.readLine();					

			while(line != null){
				try{
					String[] f = line.split("\\s+");					
					int serial;
					try{
						serial = Integer.parseInt(f[0]);	
					}catch(NumberFormatException e){
						line = br.readLine();
						continue;
					}				

					double potential = Double.parseDouble(f[f.length-1]);
					IResidue residue = atomMap.get(serial).getResidue();

					double[] potentials = potMap.get(residue);
					if (potentials == null){
						potentials = new double[]{0.0, 0};
						potMap.put(residue, potentials);
					}

					potentials[0] += potential;
					++potentials[1];										

				}catch(Exception e){	
				}

				line = br.readLine();
			}

		}catch (Exception e) {
			throw e;
		}
		return potMap;
	}

	private void preparePDB(Structure structure, File protonatedPdb) throws Exception {
		Structure protStructure = PDBIO.read(protonatedPdb, 0);

		for(SSBond ssbond : structure.getSSBonds()){
			try{
				IResidue cys1 = protStructure.getChain(ssbond.getChainID1()).getResidue(ssbond.getResnum1(), ssbond.getInsCode1());
				IResidue cys2 = protStructure.getChain(ssbond.getChainID2()).getResidue(ssbond.getResnum2(), ssbond.getInsCode2());
				cys1.setPDBName("CSS");
				cys2.setPDBName("CSS");
			}catch(Exception e){
				continue;
			}
		}

		BufferedWriter bw = new BufferedWriter(new FileWriter(protonatedPdb));
		int serial = 1;

		for(Chain chain : protStructure) {

			Set<IResidue> NGaps = new HashSet<IResidue>();
			Set<IResidue> chargedNGaps = new HashSet<IResidue>();
			Set<IResidue> chargedCGaps = new HashSet<IResidue>();

			Collections.sort(chain.getResidues());

			IResidue prevResidue = null;
			for(IResidue residue : chain.getResidues()){
				if (residue.getPDBName().equals("HIS"))
					residue.setPDBName("HIA");
				//Looking for gaps: C-gap is the last residue BEFORE gap; N-gap is the first residue AFTER gap
				if(prevResidue != null && (residue.getNumber()-prevResidue.getNumber() > 1)) {
					//Once you have located a C-gap you should look and see whether it is charged or not
					//you do this by searching within the residue for atoms OXT (or O2, which is a synonym).
					if(prevResidue.getAtom("OXT") != null || prevResidue.getAtom("O2") != null) {
						chargedCGaps.add(prevResidue);
						chargedNGaps.add(residue);
					} else {
						NGaps.add(residue);
					}
				}
				prevResidue = residue;          
			}		

			IResidue nterminal = chain.getResidues().remove(0);				
			IResidue cterminal = chain.getResidues().remove(chain.getResidues().size()-1);

			String originalResName = nterminal.getName();
			String nterminalStr = "", cterminalStr = "", pdbStr = "";
			boolean isProline = originalResName.equals("PRO");

			for(IAtom atom : nterminal.getAtoms()){
				final String name = atom.getName();
				atom.setAtomSerial(serial);

				if(name.equals("N") || name.equals("CA") || name.equals("H1") || name.equals("H2") || name.equals("H3") 
						|| name.equals("1H") || name.equals("2H") || name.equals("3H")){
					if(isProline)
						nterminal.setPDBName("PR+");
					else
						nterminal.setPDBName("BK+");				
					nterminalStr += FileConvert.toPDB(atom);								
				}
				else{
					nterminal.setPDBName(originalResName);
					nterminalStr += FileConvert.toPDB(atom);
				}				
				++serial;
			}

			IAtom[] atoms = chain.getAtoms(false, true);
			Arrays.sort(atoms);		
			for(int i = 0; i < atoms.length; i++){
				IResidue residue = atoms[i].getResidue();
				atoms[i].setAtomSerial(serial);
				++serial;

				String name = atoms[i].getName();
				originalResName = residue.getPDBName();

				//if you have a charged C-gap:
				//you change residue names to C, O, (or O1), O2 (or OXT) moving it to BK- as if they were C-termini.
				//Then you go to the corresponding N-gap and rename the N atom to N+1.
				if(name.equals("C") || name.equals("O") || name.equals("O1") || name.equals("O2") || name.equals("OXT")) {
					if(chargedCGaps.contains(residue))
						residue.setPDBName("BK-");
				} else if (name.equals("N")) {
					if(chargedNGaps.contains(residue))
						atoms[i].setName("N+1");
					//When you have a NOT charged C-gap, don't do anything to it
					//The corresponding N-gap and rename the N atom to N0
					else if(NGaps.contains(residue))
						atoms[i].setName("N0");
				}
				pdbStr += FileConvert.toPDB(atoms[i]);
				residue.setPDBName(originalResName);       
			}

			originalResName = cterminal.getName();		
			isProline = originalResName.equals("PRO");
			for(IAtom atom : cterminal.getAtoms()){			
				final String name = atom.getName();
				atom.setAtomSerial(serial);

				if(name.equals("OXT") || name.equals("O") || name.equals("C") || name.equals("O1")
						|| name.equals("O2")){
					cterminal.setPDBName("BK-");				
					cterminalStr += FileConvert.toPDB(atom);								
				}
				else{
					cterminal.setPDBName(originalResName);
					cterminalStr += FileConvert.toPDB(atom);
				}
				++serial;
			}

			bw.write(nterminalStr);
			bw.write(pdbStr);
			bw.write(cterminalStr);									
		}
		bw.close();
	}	

	private File addHydrogens(File pdbFile, File temporaryDir) throws Exception {		
		File protonatedPdb = new File(temporaryDir + File.separator+pdbFile.getName() + ".prot");	
		File reduceFile = new File(reducePath);		
		File dictionaryFile = new File(reduceDictionary);
		String[] cmdArray = new String[]{
				reduceFile.getAbsolutePath(),
				"-HIS",
				"-NOADJust",
				"-Quiet",
				"-DB",
				dictionaryFile.getAbsolutePath(), 
				pdbFile.getAbsolutePath()};

		int exitCode = -1;
		try {		
			exitCode = Utils.callWithTimeout(cmdArray, null, 5L, new FileWriter(protonatedPdb));			
		} catch (Exception e) {			
			throw e;
		}

		if (exitCode != 0)			
			throw new ExecutionException(
					String.format("Error to Adding Hydrogen to %s. Program return with code: %d", 
							pdbFile.getAbsolutePath(), exitCode)
					, null);

		return protonatedPdb;
	}
}

