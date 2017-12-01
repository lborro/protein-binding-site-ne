package br.embrapa.cnptia.gpbc.plc.main;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.LinkedHashSet;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.gpbc.plc.calculators.NanoEnvironmentCalculator;
import br.embrapa.cnptia.gpbc.plc.config.Config;
import br.embrapa.cnptia.gpbc.plc.data.MaxContactsTable;
import br.embrapa.cnptia.gpbc.plc.descriptors.protein.ProteinNanoEnvironment;


public class ProteinLigandNanoEnvironmentCalc {

	public static void main(String[] args) throws Exception {
		System.setProperty("http.proxyHost", "proxy.cnptia.embrapa.br");
		System.setProperty("http.proxyPort", "3128");

		ExecutorService pool = Executors.newFixedThreadPool(64);
		Set<Future<ProteinNanoEnvironment>> results = new LinkedHashSet<Future<ProteinNanoEnvironment>>();

		MaxContactsTable maxCon = new MaxContactsTable(Config.MAX_CONTACTS_PATH);

		//List of protein structures (PDB files) and their respective ligands (mol2 files)
		//Each line define a protein ligand complex (tab separated)
		//for instance: 10gs_protein.pdb	10gs_ligand.mol2
		BufferedReader br = new BufferedReader(new FileReader(args[0]));
		//directory where the PDB files are stored
		String proteinDir = args[1];
		//directory where the mol2 files are stored
		String ligandDir =  args[2];
		//cutoff for the nano-environment (angstrom)
		double cutoff = new Double(args[3].trim());
		
		String line;
		while ((line = br.readLine()) != null) {
			String[] tk = line.split("\t");
			File protein = null, ligand = null;
			try{
				protein = new File(proteinDir +  File.separator + tk[0].trim());
				ligand = new File(ligandDir + File.separator + tk[1].trim());
				Callable<ProteinNanoEnvironment> callable = new NanoEnvironmentCalculator(protein, ligand, maxCon, cutoff) ;
				Future<ProteinNanoEnvironment> future = pool.submit(callable);
				results.add(future);
			} catch (Exception e) {
				//System.out.println("Error calculating the protein nano-environment for the complex: " + line);
				//e.printStackTrace();
			}
		}
		br.close();

		int idx = 0;
		for (Future<ProteinNanoEnvironment> future : results) {
			try {
				ProteinNanoEnvironment nanoEnvironment = future.get();
				if(nanoEnvironment == null) continue;
				if(idx++ == 0) {
					String names[] = nanoEnvironment.getDescriptorsNames();
					System.out.print("ID");
					for(int i = 0; i < names.length; i++) System.out.print("\t" + names[i]);
				}
				System.out.println();
				Double values[] = nanoEnvironment.getDescriptorsValues();
				System.out.print(nanoEnvironment.getProteinLigandComplex().getID());
				for(int i = 0; i < values.length; i++) System.out.print("\t" + Utils.round(values[i],4));
			} catch(Exception e) {
			}

		}
		pool.shutdown();
		System.exit(0);
	}
}
