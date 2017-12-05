package br.embrapa.cnptia.gpbc.plc.main;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import br.embrapa.cnptia.gpbc.plc.calculators.ElectrostaticPotentialCalculator;
import br.embrapa.cnptia.gpbc.plc.config.Config;

public class ElectrostaticPotentialCalc {
	public static void main(String[] args) throws Exception {
		System.setProperty("http.proxyHost", "proxy.cnptia.embrapa.br");
		System.setProperty("http.proxyPort", "3128");

		ExecutorService pool = Executors.newFixedThreadPool(64);

		//List of protein structures (PDB files)
		BufferedReader br = new BufferedReader(new FileReader(args[0]));
		//directory where PDB files are stored
		String proteinDir = args[1];
		//directory where mol2  file will be stored
		String ligandDir  = args[2];
		//directory where EP file will be stored
		String outDir     = args[3]; 

		String line;
		while ((line = br.readLine()) != null) {
			String[] tk = line.split("\t");
			File pdbFile = null, mol2File = null;
			try{
				pdbFile  = new File(proteinDir +  File.separator + tk[0].trim());
				mol2File = new File(ligandDir  +  File.separator + tk[1].trim());

				Callable<Object> callable = new ElectrostaticPotentialCalculator(
						pdbFile,
						mol2File,
						Config.DELPHI_PATH, 
						Config.REDUCE_PATH, 
						Config.PARSE_RADII_PATH, 
						Config.PARSE_CHARGES_PATH, 
						Config.REDUCE_DIC_PATH, 
						Config.TMP_DIR, 
						outDir);

				pool.submit(callable);
			} catch(Exception e) {
				System.out.println("Error calculating the surface electrostatic potential for the protein-ligand complex: " + line);
			}
		}

		br.close();
		pool.shutdown();
		while (!pool.awaitTermination(60, TimeUnit.SECONDS)) {}
		System.exit(0);
	}
}
