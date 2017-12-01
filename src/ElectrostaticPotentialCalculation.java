import java.io.File;

import br.embrapa.cnptia.cbi.sdl.core.PDBIO;
import br.embrapa.cnptia.cbi.sdl.core.Structure;

public class ElectrostaticPotentialCalculation {
        public static void main(String[] args) throws Exception {
                System.setProperty("http.proxyHost", "proxy.cnptia.embrapa.br");
                System.setProperty("http.proxyPort", "3128");

                Structure structure = PDBIO.read(new File("db/protein/1a4w_protein.pdb"), 0);
                ElectrostaticPotentialCalculator epCalculator = new ElectrostaticPotentialCalculator("programs/delphi/delphi", "programs/REDUCE/reduce_old", "data/parse3_red.siz", "data/parse3_newn_reduce.crg",
                        "data/reduce_wwPDB_het_dict.txt", "tmp");
                
                epCalculator.calculate(structure);
                System.exit(0);
        }
}

