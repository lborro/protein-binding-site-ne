package br.embrapa.cnptia.gpbc.plc.descriptors.protein;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.Chain;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.gpbc.plc.structure.Protein;

public class ElectrostaticPotential extends AbstractProteinDescriptor {
	protected ElectrostaticPotential(Protein protein, String epDir) throws Exception {
		super(protein, new String [] {"EP_Surface"});
		File epFile = new File(epDir + File.separator + protein.getProteinID() + ".ep");
		readEP(epFile);
	}

	private void readEP(File epFile) throws IOException {
		BufferedReader br = new BufferedReader(new FileReader(epFile));
		String line = null;
		line = br.readLine();
		while(line != null) {
			try{
				String tk[] = line.split("\t");
				Chain chain = getProtein().getStructure().getChain(tk[0].trim());		
				IResidue residue = (AminoAcid) chain.getResidue(new Integer(tk[1]), tk[2].trim().charAt(0));
				Double value = new Double(tk[3]);
				if(!Utils.isAminoAcid(residue.getName())) continue;
				this.addDescriptorValues((AminoAcid)residue, new Double[] {value});
				line = br.readLine();
			} catch (Exception e) {
				line = br.readLine();
			}
		}
		br.close();
	}
}

