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
		//directory where the PDB files are stored
		String proteinDir = args[1];
		//directory where the EP file will be stored
		String outDir = args[2]; 

		String line;
		while ((line = br.readLine()) != null) {
			File pdbFile = new File(proteinDir +  File.separator + line.trim());
			Callable<Object> callable = new ElectrostaticPotentialCalculator(
					pdbFile, 
					Config.DELPHI_PATH, 
					Config.REDUCE_PATH, 
					Config.PARSE_RADII_PATH, 
					Config.PARSE_CHARGES_PATH, 
					Config.REDUCE_DIC_PATH, 
					Config.TMP_DIR, 
					outDir);
			
			pool.submit(callable);

		}
		br.close();
		pool.shutdown();
		while (!pool.awaitTermination(60, TimeUnit.SECONDS)) {}
		System.exit(0);
	}
}
