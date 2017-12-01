package br.embrapa.cnptia.gpbc.plc.data;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;

import br.embrapa.cnptia.cbi.sdl.core.types.AminoAcidType;

public class MaxContactsTable {

	private File maxConFile;
	private int maxCon[][];

	public MaxContactsTable(String maxConPath) throws NullPointerException, FileNotFoundException{
		if (maxConPath == null) throw new NullPointerException("Path for maxcon file can not be null!");
		maxConFile = new File(maxConPath);
		if (!maxConFile.exists()) throw new FileNotFoundException("maxcon file does not exist: " + maxConPath);
		maxCon = new int[20][14];
		readMaxconFile();
	}

	private void readMaxconFile() {
		try {
			BufferedReader br = new BufferedReader(new FileReader(maxConFile));
			String line;
			int i = 0, j = 0;
			while((line = br.readLine()) != null) {
				if(line.startsWith("residue")){
					AminoAcidType aminoAcid = AminoAcidType.getByName(line.substring(15).trim());
					i = aminoAcid.ordinal(); 
					j = 0;
				} else {
					maxCon[i][j++] = new Integer(line.substring(15).split(" ")[0]);
				}
			}
			br.close();
		} catch(Exception e) {
			e.printStackTrace();
			System.out.println("Error reading maxcon file!");
		}
	}

	public int[] getMaxContacts(AminoAcidType aaType) {
		return maxCon[aaType.ordinal()];
	}
}

