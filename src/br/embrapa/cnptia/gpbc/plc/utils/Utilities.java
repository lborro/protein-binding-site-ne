package br.embrapa.cnptia.gpbc.plc.utils;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

import org.biojava.nbio.structure.Element;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.Chain;
import br.embrapa.cnptia.cbi.sdl.core.IAtom;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.core.types.AminoAcidType;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;

public class Utilities {

	public static Element symbolToElement(String symbol){
		if(symbol.equals("C"))
			return 	Element.C;
		else if(symbol.equals("N")) 
			return Element.N;
		else if(symbol.equals("O")) 
			return Element.O;
		else if(symbol.equals("H")) 
			return Element.H;
		else if(symbol.equals("S")) 
			return Element.S;
		else if(symbol.equals("Cl")) 
			return Element.Cl;
		else if(symbol.equals("Br")) 
			return Element.Br;
		else if(symbol.equals("I")) 
			return Element.I;
		else if(symbol.equals("F")) 
			return Element.F;
		else if(symbol.equals("P")) 
			return Element.P;
		else if(symbol.equals("B")) 
			return Element.B;
		else if(symbol.equals("Fe"))
			return Element.Fe;
		else if(symbol.equals("As"))
			return Element.As;

		return null;
	}

	public static double[] calculateGeoCenter(IResidue res){
		double geoCenter[] = {0.0d, 0.0d, 0.0d};
		int n = res.getAtoms().size();
		for(IAtom atom: res.getAtoms()) {
			geoCenter[0] += atom.getX()/n;
			geoCenter[1] += atom.getY()/n;
			geoCenter[2] += atom.getZ()/n;
		}
		return geoCenter;
	}

	public static double[] calculateMassCenter(IResidue res){
		double coord[] = {0.0d, 0.0d, 0.0d};
		double totalMass = 0;
		for(IAtom atom: res.getAtoms()) {
			coord[0] += atom.getElement().getAtomicMass()*atom.getX();
			coord[1] += atom.getElement().getAtomicMass()*atom.getY();
			coord[2] += atom.getElement().getAtomicMass()*atom.getZ();
			totalMass += atom.getElement().getAtomicMass();
		}

		return new double[]{coord[0]/totalMass, coord[1]/totalMass, coord[2]/totalMass};
	}


	public static boolean isWithinSphere(IAtom atom, double[] centerXYZ, double sphereRadius) {
		double sqDist = (sphereRadius - atom.getElement().getVDWRadius()) * (sphereRadius - atom.getElement().getVDWRadius());
		return Utils.euclideanSquaredDistance(centerXYZ, atom.getCoords()) < sqDist;
	}

	public static double sphereSphereIntersectionVolume(double totalVolume, double R, double r, double dd) {
		if( dd > (R+r)*(R+r) )
			return 0.0;
		// Atom is totally inside the Sphere
		if (dd < (R-r)*(R-r))
			return totalVolume;

		final double d = Math.sqrt(dd);
		double volume = Math.PI * (R + r - d)*(R + r -d)*(d*d + 2*d*r - 3*r*r + 2*d*R + 6*r*R - 3*R*R)/(12.0*d);
		return volume;				
	}		

	public static boolean isHydrophobic(String residueName) {
		return (residueName.equals("ALA") || residueName.equals("CYS") || 
				residueName.equals("ILE") || residueName.equals("LEU") || 
				residueName.equals("MET") || residueName.equals("PHE") ||
				residueName.equals("VAL"));
	}

	public static boolean isPolar(String residueName) {
		return (residueName.equals("SER") || residueName.equals("THR") || 
				residueName.equals("TYR") || residueName.equals("ASN") || 
				residueName.equals("GLN") || residueName.equals("HIS") ||
				residueName.equals("TRP") || residueName.equals("PRO"));
	}

	public static boolean isCharged(String residueName) {
		return (residueName.equals("ASP") || residueName.equals("GLU") || 
				residueName.equals("ARG") || residueName.equals("LYS"));
	}


	public static String chainToFast(Chain chain) {
		String fasta = "";
		for (IResidue res : chain) {
			if(Utils.isAminoAcid(res.getName())) {
				AminoAcid aa = (AminoAcid) res;
				fasta += aa.getAminoAcidType().getOneLetterName();
			} else {
				String modRes =  modifiedResidues.get(res.getName());
				if(modRes != null) fasta += getOneLetterCode(modRes);
			}
		}
		return fasta;
	}

	public static char getOneLetterCode(String threeLetterCode) {
		char oneLetterCode = ' ';
		for(AminoAcidType aaType : AminoAcidType.values()) {
			if(aaType.getName().equals(threeLetterCode)) {
				oneLetterCode = aaType.getOneLetterName();
				break;
			}
		}
		return oneLetterCode;
	}

	public static final Map<String, String> modifiedResidues;
	static {
		Map<String, String> aMap = new HashMap<>();
		aMap.put("0CS","ALA");                    //  0CS ALA  3-[(S)-HYDROPEROXYSULFINYL]-L-ALANINE
		aMap.put("1AB","PRO");                    //  1AB PRO  1,4-DIDEOXY-1,4-IMINO-D-ARABINITOL
		aMap.put("1LU","LEU");                    //  1LU LEU  4-METHYL-PENTANOIC ACID-2-OXYL GROUP
		aMap.put("1PA","PHE");                    //  1PA PHE  PHENYLMETHYLACETIC ACID ALANINE
		aMap.put("1TQ","TRP");                    //  1TQ TRP  6-(FORMYLAMINO)-7-HYDROXY-L-TRYPTOPHAN
		aMap.put("1TY","TYR");                    //  1TY TYR
		aMap.put("23F","PHE");                    //  23F PHE  (2Z)-2-AMINO-3-PHENYLACRYLIC ACID
		aMap.put("23S","TRP");                    //  23S TRP  MODIFIED TRYPTOPHAN
		aMap.put("2BU","ALA");                    //  2BU ADE
		aMap.put("2ML","LEU");                    //  2ML LEU  2-METHYLLEUCINE
		aMap.put("2MR","ARG");                    //  2MR ARG  N3, N4-DIMETHYLARGININE
		aMap.put("2MT","PRO");                    //  2MT PRO
		aMap.put("2OP","ALA");                    //  2OP (2S  2-HYDROXYPROPANAL
		aMap.put("2TY","TYR");                    //  2TY TYR
		aMap.put("32S","TRP");                    //  32S TRP  MODIFIED TRYPTOPHAN
		aMap.put("32T","TRP");                    //  32T TRP  MODIFIED TRYPTOPHAN
		aMap.put("3AH","HIS");                    //  3AH HIS
		aMap.put("3MD","ASP");                    //  3MD ASP  2S,3S-3-METHYLASPARTIC ACID
		aMap.put("3TY","TYR");                    //  3TY TYR  MODIFIED TYROSINE
		aMap.put("4DP","TRP");                    //  4DP TRP
		aMap.put("4F3","ALA");                    //  4F3 ALA  CYCLIZED
		aMap.put("4FB","PRO");                    //  4FB PRO  (4S)-4-FLUORO-L-PROLINE
		aMap.put("4FW","TRP");                    //  4FW TRP  4-FLUOROTRYPTOPHANE
		aMap.put("4HT","TRP");                    //  4HT TRP  4-HYDROXYTRYPTOPHAN
		aMap.put("4IN","TRP");                    //  4IN TRP  4-AMINO-L-TRYPTOPHAN
		aMap.put("4PH","PHE");                    //  4PH PHE  4-METHYL-L-PHENYLALANINE
		aMap.put("5CS","CYS");                    //  5CS CYS
		aMap.put("6CL","LYS");                    //  6CL LYS  6-CARBOXYLYSINE
		aMap.put("6CW","TRP");                    //  6CW TRP  6-CHLORO-L-TRYPTOPHAN
		aMap.put("A0A","ASP");                    //  A0A ASP  ASPARTYL-FORMYL MIXED ANHYDRIDE
		aMap.put("AA4","ALA");                    //  AA4 ALA  2-AMINO-5-HYDROXYPENTANOIC ACID
		aMap.put("AAR","ARG");                    //  AAR ARG  ARGININEAMIDE
		aMap.put("AB7","GLU");                    //  AB7 GLU  ALPHA-AMINOBUTYRIC ACID
		aMap.put("ABA","ALA");                    //  ABA ALA  ALPHA-AMINOBUTYRIC ACID
		aMap.put("ACB","ASP");                    //  ACB ASP  3-METHYL-ASPARTIC ACID
		aMap.put("ACL","ARG");                    //  ACL ARG  DEOXY-CHLOROMETHYL-ARGININE
		aMap.put("ACY","GLY");                    //  ACY GLY  POST-TRANSLATIONAL MODIFICATION
		aMap.put("AEI","THR");                    //  AEI THR  ACYLATED THR
		aMap.put("AFA","ASN");                    //  AFA ASN  N-[7-METHYL-OCT-2,4-DIENOYL]ASPARAGINE
		aMap.put("AGM","ARG");                    //  AGM ARG  4-METHYL-ARGININE
		aMap.put("AGT","CYS");                    //  AGT CYS  AGMATINE-CYSTEINE ADDUCT
		aMap.put("AHB","ASN");                    //  AHB ASN  BETA-HYDROXYASPARAGINE
		aMap.put("AHO","ALA");                    //  AHO ALA  N-ACETYL-N-HYDROXY-L-ORNITHINE
		aMap.put("AHP","ALA");                    //  AHP ALA  2-AMINO-HEPTANOIC ACID
		aMap.put("AIB","ALA");                    //  AIB ALA  ALPHA-AMINOISOBUTYRIC ACID
		aMap.put("AKL","ASP");                    //  AKL ASP  3-AMINO-5-CHLORO-4-OXOPENTANOIC ACID
		aMap.put("ALA","ALA");                    //  ALA ALA
		aMap.put("ALC","ALA");                    //  ALC ALA  2-AMINO-3-CYCLOHEXYL-PROPIONIC ACID
		aMap.put("ALG","ARG");                    //  ALG ARG  GUANIDINOBUTYRYL GROUP
		aMap.put("ALM","ALA");                    //  ALM ALA  1-METHYL-ALANINAL
		aMap.put("ALN","ALA");                    //  ALN ALA  NAPHTHALEN-2-YL-3-ALANINE
		aMap.put("ALO","THR");                    //  ALO THR  ALLO-THREONINE
		aMap.put("ALS","ALA");                    //  ALS ALA  2-AMINO-3-OXO-4-SULFO-BUTYRIC ACID
		aMap.put("ALT","ALA");                    //  ALT ALA  THIOALANINE
		aMap.put("ALY","LYS");                    //  ALY LYS  N(6)-ACETYLLYSINE
		aMap.put("AME","MET");                    //  AME MET  ACETYLATED METHIONINE
		aMap.put("AP7","ALA");                    //  AP7 ADE
		aMap.put("APH","ALA");                    //  APH ALA  P-AMIDINOPHENYL-3-ALANINE
		aMap.put("API","LYS");                    //  API LYS  2,6-DIAMINOPIMELIC ACID
		aMap.put("APK","LYS");                    //  APK LYS
		aMap.put("AR2","ARG");                    //  AR2 ARG  ARGINYL-BENZOTHIAZOLE-6-CARBOXYLIC ACID
		aMap.put("AR4","GLU");                    //  AR4 GLU
		aMap.put("ARG","ARG");                    //  ARG ARG
		aMap.put("ARM","ARG");                    //  ARM ARG  DEOXY-METHYL-ARGININE
		aMap.put("ARO","ARG");                    //  ARO ARG  C-GAMMA-HYDROXY ARGININE
		aMap.put("ASA","ASP");                    //  ASA ASP  ASPARTIC ALDEHYDE
		aMap.put("ASB","ASP");                    //  ASB ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
		aMap.put("ASI","ASP");                    //  ASI ASP  L-ISO-ASPARTATE
		aMap.put("ASK","ASP");                    //  ASK ASP  DEHYDROXYMETHYLASPARTIC ACID
		aMap.put("ASL","ASP");                    //  ASL ASP  ASPARTIC ACID-4-CARBOXYETHYL ESTER
		aMap.put("ASN","ASN");                    //  ASN ASN
		aMap.put("ASP","ASP");                    //  ASP ASP
		aMap.put("AYA","ALA");                    //  AYA ALA  N-ACETYLALANINE
		aMap.put("AYG","ALA");                    //  AYG ALA
		aMap.put("AZK","LYS");                    //  AZK LYS  (2S)-2-AMINO-6-TRIAZANYLHEXAN-1-OL
		aMap.put("B2A","ALA");                    //  B2A ALA  ALANINE BORONIC ACID
		aMap.put("B2F","PHE");                    //  B2F PHE  PHENYLALANINE BORONIC ACID
		aMap.put("B2I","ILE");                    //  B2I ILE  ISOLEUCINE BORONIC ACID
		aMap.put("B2V","VAL");                    //  B2V VAL  VALINE BORONIC ACID
		aMap.put("B3A","ALA");                    //  B3A ALA  (3S)-3-AMINOBUTANOIC ACID
		aMap.put("B3D","ASP");                    //  B3D ASP  3-AMINOPENTANEDIOIC ACID
		aMap.put("B3E","GLU");                    //  B3E GLU  (3S)-3-AMINOHEXANEDIOIC ACID
		aMap.put("B3K","LYS");                    //  B3K LYS  (3S)-3,7-DIAMINOHEPTANOIC ACID
		aMap.put("B3S","SER");                    //  B3S SER  (3R)-3-AMINO-4-HYDROXYBUTANOIC ACID
		aMap.put("B3X","ASN");                    //  B3X ASN  (3S)-3,5-DIAMINO-5-OXOPENTANOIC ACID
		aMap.put("B3Y","TYR");                    //  B3Y TYR
		aMap.put("BAL","ALA");                    //  BAL ALA  BETA-ALANINE
		aMap.put("BBC","CYS");                    //  BBC CYS
		aMap.put("BCS","CYS");                    //  BCS CYS  BENZYLCYSTEINE
		aMap.put("BCX","CYS");                    //  BCX CYS  BETA-3-CYSTEINE
		aMap.put("BFD","ASP");                    //  BFD ASP  ASPARTATE BERYLLIUM FLUORIDE
		aMap.put("BG1","SER");                    //  BG1 SER
		aMap.put("BHD","ASP");                    //  BHD ASP  BETA-HYDROXYASPARTIC ACID
		aMap.put("BIF","PHE");                    //  BIF PHE
		aMap.put("BLE","LEU");                    //  BLE LEU  LEUCINE BORONIC ACID
		aMap.put("BLY","LYS");                    //  BLY LYS  LYSINE BORONIC ACID
		aMap.put("BMT","THR");                    //  BMT THR
		aMap.put("BNN","ALA");                    //  BNN ALA  ACETYL-P-AMIDINOPHENYLALANINE
		aMap.put("BOR","ARG");                    //  BOR ARG
		aMap.put("BOR","ARG");                    //  BOR ARG
		aMap.put("BPE","CYS");                    //  BPE CYS
		aMap.put("BTR","TRP");                    //  BTR TRP  6-BROMO-TRYPTOPHAN
		aMap.put("BUC","CYS");                    //  BUC CYS  S,S-BUTYLTHIOCYSTEINE
		aMap.put("BUG","LEU");                    //  BUG LEU  TERT-LEUCYL AMINE
		aMap.put("C12","ALA");                    //  C12 ALA
		aMap.put("C1X","LYS");                    //  C1X LYS  MODIFIED LYSINE
		aMap.put("C3Y","CYS");                    //  C3Y CYS  MODIFIED CYSTEINE
		aMap.put("C5C","CYS");                    //  C5C CYS  S-CYCLOPENTYL THIOCYSTEINE
		aMap.put("C6C","CYS");                    //  C6C CYS  S-CYCLOHEXYL THIOCYSTEINE
		aMap.put("C99","ALA");                    //  C99 ALA
		aMap.put("CAB","ALA");                    //  CAB ALA  4-CARBOXY-4-AMINOBUTANAL
		aMap.put("CAF","CYS");                    //  CAF CYS  S-DIMETHYLARSINOYL-CYSTEINE
		aMap.put("CAS","CYS");                    //  CAS CYS  S-(DIMETHYLARSENIC)CYSTEINE
		aMap.put("CCS","CYS");                    //  CCS CYS  CARBOXYMETHYLATED CYSTEINE
		aMap.put("CGU","GLU");                    //  CGU GLU  CARBOXYLATION OF THE CG ATOM
		aMap.put("CH6","ALA");                    //  CH6 ALA
		aMap.put("CH7","ALA");                    //  CH7 ALA
		aMap.put("CHG","GLY");                    //  CHG GLY  CYCLOHEXYL GLYCINE
		aMap.put("CHP","GLY");                    //  CHP GLY  3-CHLORO-4-HYDROXYPHENYLGLYCINE
		aMap.put("CHS","PHE");                    //  CHS PHE  4-AMINO-5-CYCLOHEXYL-3-HYDROXY-PENTANOIC AC
		aMap.put("CIR","ARG");                    //  CIR ARG  CITRULLINE
		aMap.put("CLB","ALA");                    //  CLB ALA
		aMap.put("CLD","ALA");                    //  CLD ALA
		aMap.put("CLE","LEU");                    //  CLE LEU  LEUCINE AMIDE
		aMap.put("CLG","LYS");                    //  CLG LYS
		aMap.put("CLH","LYS");                    //  CLH LYS
		aMap.put("CLV","ALA");                    //  CLV ALA
		aMap.put("CME","CYS");                    //  CME CYS  MODIFIED CYSTEINE
		aMap.put("CML","CYS");                    //  CML CYS
		aMap.put("CMT","CYS");                    //  CMT CYS  O-METHYLCYSTEINE
		aMap.put("CQR","ALA");                    //  CQR ALA
		aMap.put("CR2","ALA");                    //  CR2 ALA  POST-TRANSLATIONAL MODIFICATION
		aMap.put("CR5","ALA");                    //  CR5 ALA
		aMap.put("CR7","ALA");                    //  CR7 ALA
		aMap.put("CR8","ALA");                    //  CR8 ALA
		aMap.put("CRK","ALA");                    //  CRK ALA
		aMap.put("CRO","THR");                    //  CRO THR  CYCLIZED
		aMap.put("CRQ","TYR");                    //  CRQ TYR
		aMap.put("CRW","ALA");                    //  CRW ALA
		aMap.put("CRX","ALA");                    //  CRX ALA
		aMap.put("CS1","CYS");                    //  CS1 CYS  S-(2-ANILINYL-SULFANYL)-CYSTEINE
		aMap.put("CS3","CYS");                    //  CS3 CYS
		aMap.put("CS4","CYS");                    //  CS4 CYS
		aMap.put("CSA","CYS");                    //  CSA CYS  S-ACETONYLCYSTEIN
		aMap.put("CSB","CYS");                    //  CSB CYS  CYS BOUND TO LEAD ION
		aMap.put("CSD","CYS");                    //  CSD CYS  3-SULFINOALANINE
		aMap.put("CSE","CYS");                    //  CSE CYS  SELENOCYSTEINE
		aMap.put("CSI","ALA");                    //  CSI ALA
		aMap.put("CSO","CYS");                    //  CSO CYS  INE S-HYDROXYCYSTEINE
		aMap.put("CSR","CYS");                    //  CSR CYS  S-ARSONOCYSTEINE
		aMap.put("CSS","CYS");                    //  CSS CYS  1,3-THIAZOLE-4-CARBOXYLIC ACID
		aMap.put("CSU","CYS");                    //  CSU CYS  CYSTEINE-S-SULFONIC ACID
		aMap.put("CSW","CYS");                    //  CSW CYS  CYSTEINE-S-DIOXIDE
		aMap.put("CSX","CYS");                    //  CSX CYS  OXOCYSTEINE
		aMap.put("CSY","ALA");                    //  CSY ALA  MODIFIED TYROSINE COMPLEX
		aMap.put("CSZ","CYS");                    //  CSZ CYS  S-SELANYL CYSTEINE
		aMap.put("CTH","THR");                    //  CTH THR  4-CHLOROTHREONINE
		aMap.put("CWR","ALA");                    //  CWR ALA
		aMap.put("CXM","MET");                    //  CXM MET  N-CARBOXYMETHIONINE
		aMap.put("CY0","CYS");                    //  CY0 CYS  MODIFIED CYSTEINE
		aMap.put("CY1","CYS");                    //  CY1 CYS  ACETAMIDOMETHYLCYSTEINE
		aMap.put("CY3","CYS");                    //  CY3 CYS  2-AMINO-3-MERCAPTO-PROPIONAMIDE
		aMap.put("CY4","CYS");                    //  CY4 CYS  S-BUTYRYL-CYSTEIN
		aMap.put("CY7","CYS");                    //  CY7 CYS  MODIFIED CYSTEINE
		aMap.put("CYD","CYS");                    //  CYD CYS
		aMap.put("CYF","CYS");                    //  CYF CYS  FLUORESCEIN LABELLED CYS380 (P14)
		aMap.put("CYG","CYS");                    //  CYG CYS
		aMap.put("CYJ","LYS");                    //  CYJ LYS  MODIFIED LYSINE
		aMap.put("CYQ","CYS");                    //  CYQ CYS
		aMap.put("CYR","CYS");                    //  CYR CYS
		aMap.put("CYS","CYS");                    //  CYS CYS
		aMap.put("CYX","CYS");                    //  CYX CYS  Not an official code, but used by some programs as disulfide residues
		aMap.put("CZ2","CYS");                    //  CZ2 CYS  S-(DIHYDROXYARSINO)CYSTEINE
		aMap.put("CZZ","CYS");                    //  CZZ CYS  THIARSAHYDROXY-CYSTEINE
		aMap.put("DA2","ARG");                    //  DA2 ARG  MODIFIED ARGININE
		aMap.put("DAB","ALA");                    //  DAB ALA  2,4-DIAMINOBUTYRIC ACID
		aMap.put("DAH","PHE");                    //  DAH PHE  3,4-DIHYDROXYDAHNYLALANINE
		aMap.put("DAL","ALA");                    //  DAL ALA  D-ALANINE
		aMap.put("DAM","ALA");                    //  DAM ALA  N-METHYL-ALPHA-BETA-DEHYDROALANINE
		aMap.put("DAR","ARG");                    //  DAR ARG  D-ARGININE
		aMap.put("DAS","ASP");                    //  DAS ASP  D-ASPARTIC ACID
		aMap.put("DBU","ALA");                    //  DBU ALA  (2E)-2-AMINOBUT-2-ENOIC ACID
		aMap.put("DBY","TYR");                    //  DBY TYR  3,5 DIBROMOTYROSINE
		aMap.put("DBZ","ALA");                    //  DBZ ALA  3-(BENZOYLAMINO)-L-ALANINE
		aMap.put("DCL","LEU");                    //  DCL LEU  2-AMINO-4-METHYL-PENTANYL GROUP
		aMap.put("DCY","CYS");                    //  DCY CYS  D-CYSTEINE
		aMap.put("DDE","HIS");                    //  DDE HIS
		aMap.put("DGL","GLU");                    //  DGL GLU  D-GLU
		aMap.put("DGN","GLN");                    //  DGN GLN  D-GLUTAMINE
		aMap.put("DHA","ALA");                    //  DHA ALA  2-AMINO-ACRYLIC ACID
		aMap.put("DHI","HIS");                    //  DHI HIS  D-HISTIDINE
		aMap.put("DHL","SER");                    //  DHL SER  POST-TRANSLATIONAL MODIFICATION
		aMap.put("DIL","ILE");                    //  DIL ILE  D-ISOLEUCINE
		aMap.put("DIV","VAL");                    //  DIV VAL  D-ISOVALINE
		aMap.put("DLE","LEU");                    //  DLE LEU  D-LEUCINE
		aMap.put("DLS","LYS");                    //  DLS LYS  DI-ACETYL-LYSINE
		aMap.put("DLY","LYS");                    //  DLY LYS  D-LYSINE
		aMap.put("DMH","ASN");                    //  DMH ASN  N4,N4-DIMETHYL-ASPARAGINE
		aMap.put("DMK","ASP");                    //  DMK ASP  DIMETHYL ASPARTIC ACID
		aMap.put("DNE","LEU");                    //  DNE LEU  D-NORLEUCINE
		aMap.put("DNG","LEU");                    //  DNG LEU  N-FORMYL-D-NORLEUCINE
		aMap.put("DNL","LYS");                    //  DNL LYS  6-AMINO-HEXANAL
		aMap.put("DNM","LEU");                    //  DNM LEU  D-N-METHYL NORLEUCINE
		aMap.put("DPH","PHE");                    //  DPH PHE  DEAMINO-METHYL-PHENYLALANINE
		aMap.put("DPL","PRO");                    //  DPL PRO  4-OXOPROLINE
		aMap.put("DPN","PHE");                    //  DPN PHE  D-CONFIGURATION
		aMap.put("DPP","ALA");                    //  DPP ALA  DIAMMINOPROPANOIC ACID
		aMap.put("DPQ","TYR");                    //  DPQ TYR  TYROSINE DERIVATIVE
		aMap.put("DPR","PRO");                    //  DPR PRO  D-PROLINE
		aMap.put("DSE","SER");                    //  DSE SER  D-SERINE N-METHYLATED
		aMap.put("DSG","ASN");                    //  DSG ASN  D-ASPARAGINE
		aMap.put("DSN","SER");                    //  DSN SER  D-SERINE
		aMap.put("DTH","THR");                    //  DTH THR  D-THREONINE
		aMap.put("DTR","TRP");                    //  DTR TRP  D-TRYPTOPHAN
		aMap.put("DTY","TYR");                    //  DTY TYR  D-TYROSINE
		aMap.put("DVA","VAL");                    //  DVA VAL  D-VALINE
		aMap.put("DYG","ALA");                    //  DYG ALA
		aMap.put("DYS","CYS");                    //  DYS CYS
		aMap.put("EFC","CYS");                    //  EFC CYS  S,S-(2-FLUOROETHYL)THIOCYSTEINE
		aMap.put("ESB","TYR");                    //  ESB TYR
		aMap.put("ESC","MET");                    //  ESC MET  2-AMINO-4-ETHYL SULFANYL BUTYRIC ACID
		aMap.put("FCL","PHE");                    //  FCL PHE  3-CHLORO-L-PHENYLALANINE
		aMap.put("FGL","ALA");                    //  FGL ALA  2-AMINOPROPANEDIOIC ACID
		aMap.put("FGP","SER");                    //  FGP SER
		aMap.put("FHL","LYS");                    //  FHL LYS  MODIFIED LYSINE
		aMap.put("FLE","LEU");                    //  FLE LEU  FUROYL-LEUCINE
		aMap.put("FLT","TYR");                    //  FLT TYR  FLUOROMALONYL TYROSINE
		aMap.put("FME","MET");                    //  FME MET  FORMYL-METHIONINE
		aMap.put("FOE","CYS");                    //  FOE CYS
		aMap.put("FOG","PHE");                    //  FOG PHE  PHENYLALANINOYL-[1-HYDROXY]-2-PROPYLENE
		aMap.put("FOR","MET");                    //  FOR MET
		aMap.put("FRF","PHE");                    //  FRF PHE  PHE FOLLOWED BY REDUCED PHE
		aMap.put("FTR","TRP");                    //  FTR TRP  FLUOROTRYPTOPHANE
		aMap.put("FTY","TYR");                    //  FTY TYR  DEOXY-DIFLUOROMETHELENE-PHOSPHOTYROSINE
		aMap.put("GHG","GLN");                    //  GHG GLN  GAMMA-HYDROXY-GLUTAMINE
		aMap.put("GHP","GLY");                    //  GHP GLY  4-HYDROXYPHENYLGLYCINE
		aMap.put("GL3","GLY");                    //  GL3 GLY  POST-TRANSLATIONAL MODIFICATION
		aMap.put("GLH","GLN");                    //  GLH GLN
		aMap.put("GLN","GLN");                    //  GLN GLN
		aMap.put("GLU","GLU");                    //  GLU GLU
		aMap.put("GLY","GLY");                    //  GLY GLY
		aMap.put("GLZ","GLY");                    //  GLZ GLY  AMINO-ACETALDEHYDE
		aMap.put("GMA","GLU");                    //  GMA GLU  1-AMIDO-GLUTAMIC ACID
		aMap.put("GMU","ALA");                    //  GMU 5MU
		aMap.put("GPL","LYS");                    //  GPL LYS
		aMap.put("GT9","CYS");                    //  GT9 CYS  SG ALKYLATED
		aMap.put("GVL","SER");                    //  GVL SER  SERINE MODIFED WITH PHOSPHOPANTETHEINE
		aMap.put("GYC","CYS");                    //  GYC CYS
		aMap.put("GYS","GLY");                    //  GYS GLY
		aMap.put("H5M","PRO");                    //  H5M PRO  TRANS-3-HYDROXY-5-METHYLPROLINE
		aMap.put("HHK","ALA");                    //  HHK ALA  (2S)-2,8-DIAMINOOCTANOIC ACID
		aMap.put("HIA","HIS");                    //  HIA HIS  L-HISTIDINE AMIDE
		aMap.put("HIC","HIS");                    //  HIC HIS  4-METHYL-HISTIDINE
		aMap.put("HID","HIS");                    //  HID HIS  Not an official code, but used by some programs as delta-protonated histidine
		aMap.put("HIE","HIS");                    //  HIE HIS  Not an official code, but used by some programs as epsilon-protonated histidine
		aMap.put("HIP","HIS");                    //  HIP HIS  ND1-PHOSPHONOHISTIDINE, also non-officially as protonated histidine
		aMap.put("HIQ","HIS");                    //  HIQ HIS  MODIFIED HISTIDINE
		aMap.put("HIS","HIS");                    //  HIS HIS
		aMap.put("HLU","LEU");                    //  HLU LEU  BETA-HYDROXYLEUCINE
		aMap.put("HMF","ALA");                    //  HMF ALA  2-AMINO-4-PHENYL-BUTYRIC ACID
		aMap.put("HMR","ARG");                    //  HMR ARG  BETA-HOMOARGININE
		aMap.put("HPE","PHE");                    //  HPE PHE  HOMOPHENYLALANINE
		aMap.put("HPH","PHE");                    //  HPH PHE  PHENYLALANINOL GROUP
		aMap.put("HPQ","PHE");                    //  HPQ PHE  HOMOPHENYLALANINYLMETHANE
		aMap.put("HRG","ARG");                    //  HRG ARG  L-HOMOARGININE
		aMap.put("HSE","SER");                    //  HSE SER  L-HOMOSERINE
		aMap.put("HSL","SER");                    //  HSL SER  HOMOSERINE LACTONE
		aMap.put("HSO","HIS");                    //  HSO HIS  HISTIDINOL
		aMap.put("HTI","CYS");                    //  HTI CYS
		aMap.put("HTR","TRP");                    //  HTR TRP  BETA-HYDROXYTRYPTOPHANE
		aMap.put("HY3","PRO");                    //  HY3 PRO  3-HYDROXYPROLINE
		aMap.put("HYP","PRO");                    //  HYP PRO  4-HYDROXYPROLINE
		aMap.put("IAM","ALA");                    //  IAM ALA  4-[(ISOPROPYLAMINO)METHYL]PHENYLALANINE
		aMap.put("IAS","ASP");                    //  IAS ASP  ASPARTYL GROUP
		aMap.put("IGL","ALA");                    //  IGL ALA  ALPHA-AMINO-2-INDANACETIC ACID
		aMap.put("IIL","ILE");                    //  IIL ILE  ISO-ISOLEUCINE
		aMap.put("ILE","ILE");                    //  ILE ILE
		aMap.put("ILG","GLU");                    //  ILG GLU  GLU LINKED TO NEXT RESIDUE VIA CG
		aMap.put("ILX","ILE");                    //  ILX ILE  4,5-DIHYDROXYISOLEUCINE
		aMap.put("IML","ILE");                    //  IML ILE  N-METHYLATED
		aMap.put("IPG","GLY");                    //  IPG GLY  N-ISOPROPYL GLYCINE
		aMap.put("IT1","LYS");                    //  IT1 LYS
		aMap.put("IYR","TYR");                    //  IYR TYR  3-IODO-TYROSINE
		aMap.put("KCX","LYS");                    //  KCX LYS  CARBAMOYLATED LYSINE
		aMap.put("KGC","LYS");                    //  KGC LYS
		aMap.put("KOR","CYS");                    //  KOR CYS  MODIFIED CYSTEINE
		aMap.put("KST","LYS");                    //  KST LYS  N~6~-(5-CARBOXY-3-THIENYL)-L-LYSINE
		aMap.put("KYN","ALA");                    //  KYN ALA  KYNURENINE
		aMap.put("LA2","LYS");                    //  LA2 LYS
		aMap.put("LAL","ALA");                    //  LAL ALA  N,N-DIMETHYL-L-ALANINE
		aMap.put("LCK","LYS");                    //  LCK LYS
		aMap.put("LCX","LYS");                    //  LCX LYS  CARBAMYLATED LYSINE
		aMap.put("LDH","LYS");                    //  LDH LYS  N~6~-ETHYL-L-LYSINE
		aMap.put("LED","LEU");                    //  LED LEU  POST-TRANSLATIONAL MODIFICATION
		aMap.put("LEF","LEU");                    //  LEF LEU  2-5-FLUOROLEUCINE
		aMap.put("LET","LYS");                    //  LET LYS  ODIFIED LYSINE
		aMap.put("LEU","LEU");                    //  LEU LEU
		aMap.put("LLP","LYS");                    //  LLP LYS
		aMap.put("LLY","LYS");                    //  LLY LYS  NZ-(DICARBOXYMETHYL)LYSINE
		aMap.put("LME","GLU");                    //  LME GLU  (3R)-3-METHYL-L-GLUTAMIC ACID
		aMap.put("LNT","LEU");                    //  LNT LEU
		aMap.put("LPD","PRO");                    //  LPD PRO  L-PROLINAMIDE
		aMap.put("LSO","LYS");                    //  LSO LYS  MODIFIED LYSINE
		aMap.put("LYM","LYS");                    //  LYM LYS  DEOXY-METHYL-LYSINE
		aMap.put("LYN","LYS");                    //  LYN LYS  2,6-DIAMINO-HEXANOIC ACID AMIDE
		aMap.put("LYP","LYS");                    //  LYP LYS  N~6~-METHYL-N~6~-PROPYL-L-LYSINE
		aMap.put("LYR","LYS");                    //  LYR LYS  MODIFIED LYSINE
		aMap.put("LYS","LYS");                    //  LYS LYS
		aMap.put("LYX","LYS");                    //  LYX LYS  
		aMap.put("LYZ","LYS");                    //  LYZ LYS  5-HYDROXYLYSINE
		aMap.put("M0H","CYS");                    //  M0H CYS  S-(HYDROXYMETHYL)-L-CYSTEINE
		aMap.put("M2L","LYS");                    //  M2L LYS
		aMap.put("M3L","LYS");                    //  M3L LYS  N-TRIMETHYLLYSINE
		aMap.put("MAA","ALA");                    //  MAA ALA  N-METHYLALANINE
		aMap.put("MAI","ARG");                    //  MAI ARG  DEOXO-METHYLARGININE
		aMap.put("MBQ","TYR");                    //  MBQ TYR
		aMap.put("MC1","SER");                    //  MC1 SER  METHICILLIN ACYL-SERINE
		aMap.put("MCL","LYS");                    //  MCL LYS  NZ-(1-CARBOXYETHYL)-LYSINE
		aMap.put("MCS","CYS");                    //  MCS CYS  MALONYLCYSTEINE
		aMap.put("MDO","ALA");                    //  MDO ALA
		aMap.put("MEA","PHE");                    //  MEA PHE  N-METHYLPHENYLALANINE
		aMap.put("MEG","GLU");                    //  MEG GLU  (2S,3R)-3-METHYL-GLUTAMIC ACID
		aMap.put("MEN","ASN");                    //  MEN ASN  GAMMA METHYL ASPARAGINE
		aMap.put("MET","MET");                    //  MET MET
		aMap.put("MEU","GLY");                    //  MEU GLY  O-METHYL-GLYCINE
		aMap.put("MFC","ALA");                    //  MFC ALA  CYCLIZED
		aMap.put("MGG","ARG");                    //  MGG ARG  MODIFIED D-ARGININE
		aMap.put("MGN","GLN");                    //  MGN GLN  2-METHYL-GLUTAMINE
		aMap.put("MHL","LEU");                    //  MHL LEU  N-METHYLATED, HYDROXY
		aMap.put("MHO","MET");                    //  MHO MET  POST-TRANSLATIONAL MODIFICATION
		aMap.put("MHS","HIS");                    //  MHS HIS  1-N-METHYLHISTIDINE
		aMap.put("MIS","SER");                    //  MIS SER  MODIFIED SERINE
		aMap.put("MLE","LEU");                    //  MLE LEU  N-METHYLATED
		aMap.put("MLL","LEU");                    //  MLL LEU  METHYL L-LEUCINATE
		aMap.put("MLY","LYS");                    //  MLY LYS  METHYLATED LYSINE
		aMap.put("MLZ","LYS");                    //  MLZ LYS  N-METHYL-LYSINE
		aMap.put("MME","MET");                    //  MME MET  N-METHYL METHIONINE
		aMap.put("MNL","LEU");                    //  MNL LEU  4,N-DIMETHYLNORLEUCINE
		aMap.put("MNV","VAL");                    //  MNV VAL  N-METHYL-C-AMINO VALINE
		aMap.put("MPQ","GLY");                    //  MPQ GLY  N-METHYL-ALPHA-PHENYL-GLYCINE
		aMap.put("MSA","GLY");                    //  MSA GLY  (2-S-METHYL) SARCOSINE
		aMap.put("MSE","MET");                    //  MSE MET  ELENOMETHIONINE
		aMap.put("MSO","MET");                    //  MSO MET  METHIONINE SULFOXIDE
		aMap.put("MTY","PHE");                    //  MTY PHE  3-HYDROXYPHENYLALANINE
		aMap.put("MVA","VAL");                    //  MVA VAL  N-METHYLATED
		aMap.put("N10","SER");                    //  N10 SER  O-[(HEXYLAMINO)CARBONYL]-L-SERINE
		aMap.put("NAL","ALA");                    //  NAL ALA  BETA-(2-NAPHTHYL)-ALANINE
		aMap.put("NAM","ALA");                    //  NAM ALA  NAM NAPTHYLAMINOALANINE
		aMap.put("NBQ","TYR");                    //  NBQ TYR
		aMap.put("NC1","SER");                    //  NC1 SER  NITROCEFIN ACYL-SERINE
		aMap.put("NCB","ALA");                    //  NCB ALA  CHEMICAL MODIFICATION
		aMap.put("NEP","HIS");                    //  NEP HIS  N1-PHOSPHONOHISTIDINE
		aMap.put("NFA","PHE");                    //  NFA PHE  MODIFIED PHENYLALANINE
		aMap.put("NIY","TYR");                    //  NIY TYR  META-NITRO-TYROSINE
		aMap.put("NLE","LEU");                    //  NLE LEU  NORLEUCINE
		aMap.put("NLN","LEU");                    //  NLN LEU  NORLEUCINE AMIDE
		aMap.put("NLO","LEU");                    //  NLO LEU  O-METHYL-L-NORLEUCINE
		aMap.put("NMC","GLY");                    //  NMC GLY  N-CYCLOPROPYLMETHYL GLYCINE
		aMap.put("NMM","ARG");                    //  NMM ARG  MODIFIED ARGININE
		aMap.put("NPH","CYS");                    //  NPH CYS
		aMap.put("NRQ","ALA");                    //  NRQ ALA
		aMap.put("NVA","VAL");                    //  NVA VAL  NORVALINE
		aMap.put("NYC","ALA");                    //  NYC ALA
		aMap.put("NYS","CYS");                    //  NYS CYS
		aMap.put("NZH","HIS");                    //  NZH HIS
		aMap.put("OAS","SER");                    //  OAS SER  O-ACETYLSERINE
		aMap.put("OBS","LYS");                    //  OBS LYS  MODIFIED LYSINE
		aMap.put("OCS","CYS");                    //  OCS CYS  CYSTEINE SULFONIC ACID
		aMap.put("OCY","CYS");                    //  OCY CYS  HYDROXYETHYLCYSTEINE
		aMap.put("OHI","HIS");                    //  OHI HIS  3-(2-OXO-2H-IMIDAZOL-4-YL)-L-ALANINE
		aMap.put("OHS","ASP");                    //  OHS ASP  O-(CARBOXYSULFANYL)-4-OXO-L-HOMOSERINE
		aMap.put("OLT","THR");                    //  OLT THR  O-METHYL-L-THREONINE
		aMap.put("OMT","MET");                    //  OMT MET  METHIONINE SULFONE
		aMap.put("OPR","ARG");                    //  OPR ARG  C-(3-OXOPROPYL)ARGININE
		aMap.put("ORN","ALA");                    //  ORN ALA  ORNITHINE
		aMap.put("ORQ","ARG");                    //  ORQ ARG  N~5~-ACETYL-L-ORNITHINE
		aMap.put("OSE","SER");                    //  OSE SER  O-SULFO-L-SERINE
		aMap.put("OTY","TYR");                    //  OTY TYR
		aMap.put("OXX","ASP");                    //  OXX ASP  OXALYL-ASPARTYL ANHYDRIDE
		aMap.put("P1L","CYS");                    //  P1L CYS  S-PALMITOYL CYSTEINE
		aMap.put("P2Y","PRO");                    //  P2Y PRO  (2S)-PYRROLIDIN-2-YLMETHYLAMINE
		aMap.put("PAQ","TYR");                    //  PAQ TYR  SEE REMARK 999
		aMap.put("PAT","TRP");                    //  PAT TRP  ALPHA-PHOSPHONO-TRYPTOPHAN
		aMap.put("PBB","CYS");                    //  PBB CYS  S-(4-BROMOBENZYL)CYSTEINE
		aMap.put("PBF","PHE");                    //  PBF PHE  PARA-(BENZOYL)-PHENYLALANINE
		aMap.put("PCA","PRO");                    //  PCA PRO  5-OXOPROLINE
		aMap.put("PCS","PHE");                    //  PCS PHE  PHENYLALANYLMETHYLCHLORIDE
		aMap.put("PEC","CYS");                    //  PEC CYS  S,S-PENTYLTHIOCYSTEINE
		aMap.put("PF5","PHE");                    //  PF5 PHE  2,3,4,5,6-PENTAFLUORO-L-PHENYLALANINE
		aMap.put("PFF","PHE");                    //  PFF PHE  4-FLUORO-L-PHENYLALANINE
		aMap.put("PG1","SER");                    //  PG1 SER  BENZYLPENICILLOYL-ACYLATED SERINE
		aMap.put("PG9","GLY");                    //  PG9 GLY  D-PHENYLGLYCINE
		aMap.put("PHA","PHE");                    //  PHA PHE  PHENYLALANINAL
		aMap.put("PHD","ASP");                    //  PHD ASP  2-AMINO-4-OXO-4-PHOSPHONOOXY-BUTYRIC ACID
		aMap.put("PHE","PHE");                    //  PHE PHE
		aMap.put("PHI","PHE");                    //  PHI PHE  IODO-PHENYLALANINE
		aMap.put("PHL","PHE");                    //  PHL PHE  L-PHENYLALANINOL
		aMap.put("PHM","PHE");                    //  PHM PHE  PHENYLALANYLMETHANE
		aMap.put("PIA","ALA");                    //  PIA ALA  FUSION OF ALA 65, TYR 66, GLY 67
		aMap.put("PLE","LEU");                    //  PLE LEU  LEUCINE PHOSPHINIC ACID
		aMap.put("PM3","PHE");                    //  PM3 PHE
		aMap.put("POM","PRO");                    //  POM PRO  CIS-5-METHYL-4-OXOPROLINE
		aMap.put("PPH","LEU");                    //  PPH LEU  PHENYLALANINE PHOSPHINIC ACID
		aMap.put("PPN","PHE");                    //  PPN PHE  THE LIGAND IS A PARA-NITRO-PHENYLALANINE
		aMap.put("PR3","CYS");                    //  PR3 CYS  INE DTT-CYSTEINE
		aMap.put("PRO","PRO");                    //  PRO PRO
		aMap.put("PRQ","PHE");                    //  PRQ PHE  PHENYLALANINE
		aMap.put("PRR","ALA");                    //  PRR ALA  3-(METHYL-PYRIDINIUM)ALANINE
		aMap.put("PRS","PRO");                    //  PRS PRO  THIOPROLINE
		aMap.put("PSA","PHE");                    //  PSA PHE
		aMap.put("PSH","HIS");                    //  PSH HIS  1-THIOPHOSPHONO-L-HISTIDINE
		aMap.put("PTH","TYR");                    //  PTH TYR  METHYLENE-HYDROXY-PHOSPHOTYROSINE
		aMap.put("PTM","TYR");                    //  PTM TYR  ALPHA-METHYL-O-PHOSPHOTYROSINE
		aMap.put("PTR","TYR");                    //  PTR TYR  O-PHOSPHOTYROSINE
		aMap.put("PYA","ALA");                    //  PYA ALA  3-(1,10-PHENANTHROL-2-YL)-L-ALANINE
		aMap.put("PYC","ALA");                    //  PYC ALA  PYRROLE-2-CARBOXYLATE
		aMap.put("PYR","SER");                    //  PYR SER  CHEMICALLY MODIFIED
		aMap.put("PYT","ALA");                    //  PYT ALA  MODIFIED ALANINE
		aMap.put("PYX","CYS");                    //  PYX CYS  S-[S-THIOPYRIDOXAMINYL]CYSTEINE
		aMap.put("R1A","CYS");                    //  R1A CYS
		aMap.put("R1B","CYS");                    //  R1B CYS
		aMap.put("R1F","CYS");                    //  R1F CYS
		aMap.put("R7A","CYS");                    //  R7A CYS
		aMap.put("RC7","ALA");                    //  RC7 ALA
		aMap.put("RCY","CYS");                    //  RCY CYS
		aMap.put("S1H","SER");                    //  S1H SER  1-HEXADECANOSULFONYL-O-L-SERINE
		aMap.put("SAC","SER");                    //  SAC SER  N-ACETYL-SERINE
		aMap.put("SAH","CYS");                    //  SAH CYS  S-ADENOSYL-L-HOMOCYSTEINE
		aMap.put("SAR","GLY");                    //  SAR GLY  SARCOSINE
		aMap.put("SBD","SER");                    //  SBD SER
		aMap.put("SBG","SER");                    //  SBG SER  MODIFIED SERINE
		aMap.put("SBL","SER");                    //  SBL SER
		aMap.put("SC2","CYS");                    //  SC2 CYS  N-ACETYL-L-CYSTEINE
		aMap.put("SCH","CYS");                    //  SCH CYS  S-METHYL THIOCYSTEINE GROUP
		aMap.put("SCS","CYS");                    //  SCS CYS  MODIFIED CYSTEINE
		aMap.put("SCY","CYS");                    //  SCY CYS  CETYLATED CYSTEINE
		aMap.put("SDP","SER");                    //  SDP SER
		aMap.put("SEB","SER");                    //  SEB SER  O-BENZYLSULFONYL-SERINE
		aMap.put("SEC","ALA");                    //  SEC ALA  2-AMINO-3-SELENINO-PROPIONIC ACID
		aMap.put("SEL","SER");                    //  SEL SER  2-AMINO-1,3-PROPANEDIOL
		aMap.put("SEP","SER");                    //  SEP SER  E PHOSPHOSERINE
		aMap.put("SER","SER");                    //  SER SER
		aMap.put("SET","SER");                    //  SET SER  AMINOSERINE
		aMap.put("SGB","SER");                    //  SGB SER  MODIFIED SERINE
		aMap.put("SGR","SER");                    //  SGR SER  MODIFIED SERINE
		aMap.put("SHC","CYS");                    //  SHC CYS  S-HEXYLCYSTEINE
		aMap.put("SHP","GLY");                    //  SHP GLY  (4-HYDROXYMALTOSEPHENYL)GLYCINE
		aMap.put("SIC","ALA");                    //  SIC ALA
		aMap.put("SLZ","LYS");                    //  SLZ LYS  L-THIALYSINE
		aMap.put("SMC","CYS");                    //  SMC CYS  POST-TRANSLATIONAL MODIFICATION
		aMap.put("SME","MET");                    //  SME MET  METHIONINE SULFOXIDE
		aMap.put("SMF","PHE");                    //  SMF PHE  4-SULFOMETHYL-L-PHENYLALANINE
		aMap.put("SNC","CYS");                    //  SNC CYS  S-NITROSO CYSTEINE
		aMap.put("SNN","ASP");                    //  SNN ASP  POST-TRANSLATIONAL MODIFICATION
		aMap.put("SOC","CYS");                    //  SOC CYS  DIOXYSELENOCYSTEINE
		aMap.put("SOY","SER");                    //  SOY SER  OXACILLOYL-ACYLATED SERINE
		aMap.put("SUI","ALA");                    //  SUI ALA
		aMap.put("SUN","SER");                    //  SUN SER  TABUN CONJUGATED SERINE
		aMap.put("SVA","SER");                    //  SVA SER  SERINE VANADATE
		aMap.put("SVV","SER");                    //  SVV SER  MODIFIED SERINE
		aMap.put("SVX","SER");                    //  SVX SER  MODIFIED SERINE
		aMap.put("SVY","SER");                    //  SVY SER  MODIFIED SERINE
		aMap.put("SVZ","SER");                    //  SVZ SER  MODIFIED SERINE
		aMap.put("SXE","SER");                    //  SXE SER  MODIFIED SERINE
		aMap.put("TBG","GLY");                    //  TBG GLY  T-BUTYL GLYCINE
		aMap.put("TBM","THR");                    //  TBM THR
		aMap.put("TCQ","TYR");                    //  TCQ TYR  MODIFIED TYROSINE
		aMap.put("TEE","CYS");                    //  TEE CYS  POST-TRANSLATIONAL MODIFICATION
		aMap.put("TH5","THR");                    //  TH5 THR  O-ACETYL-L-THREONINE
		aMap.put("THC","THR");                    //  THC THR  N-METHYLCARBONYLTHREONINE
		aMap.put("THR","THR");                    //  THR THR
		aMap.put("TIH","ALA");                    //  TIH ALA  BETA(2-THIENYL)ALANINE
		aMap.put("TMD","THR");                    //  TMD THR  N-METHYLATED, EPSILON C ALKYLATED
		aMap.put("TNB","CYS");                    //  TNB CYS  S-(2,3,6-TRINITROPHENYL)CYSTEINE
		aMap.put("TOX","TRP");                    //  TOX TRP
		aMap.put("TPL","TRP");                    //  TPL TRP  TRYTOPHANOL
		aMap.put("TPO","THR");                    //  TPO THR  HOSPHOTHREONINE
		aMap.put("TPQ","ALA");                    //  TPQ ALA  2,4,5-TRIHYDROXYPHENYLALANINE
		aMap.put("TQQ","TRP");                    //  TQQ TRP
		aMap.put("TRF","TRP");                    //  TRF TRP  N1-FORMYL-TRYPTOPHAN
		aMap.put("TRN","TRP");                    //  TRN TRP  AZA-TRYPTOPHAN
		aMap.put("TRO","TRP");                    //  TRO TRP  2-HYDROXY-TRYPTOPHAN
		aMap.put("TRP","TRP");                    //  TRP TRP
		aMap.put("TRQ","TRP");                    //  TRQ TRP
		aMap.put("TRW","TRP");                    //  TRW TRP
		aMap.put("TRX","TRP");                    //  TRX TRP  6-HYDROXYTRYPTOPHAN
		aMap.put("TTQ","TRP");                    //  TTQ TRP  6-AMINO-7-HYDROXY-L-TRYPTOPHAN
		aMap.put("TTS","TYR");                    //  TTS TYR
		aMap.put("TY2","TYR");                    //  TY2 TYR  3-AMINO-L-TYROSINE
		aMap.put("TY3","TYR");                    //  TY3 TYR  3-HYDROXY-L-TYROSINE
		aMap.put("TYB","TYR");                    //  TYB TYR  TYROSINAL
		aMap.put("TYC","TYR");                    //  TYC TYR  L-TYROSINAMIDE
		aMap.put("TYI","TYR");                    //  TYI TYR  3,5-DIIODOTYROSINE
		aMap.put("TYN","TYR");                    //  TYN TYR  ADDUCT AT HYDROXY GROUP
		aMap.put("TYO","TYR");                    //  TYO TYR
		aMap.put("TYQ","TYR");                    //  TYQ TYR  AMINOQUINOL FORM OF TOPA QUINONONE
		aMap.put("TYR","TYR");                    //  TYR TYR
		aMap.put("TYS","TYR");                    //  TYS TYR  INE SULPHONATED TYROSINE
		aMap.put("TYT","TYR");                    //  TYT TYR
		aMap.put("TYX","CYS");                    //  TYX CYS  S-(2-ANILINO-2-OXOETHYL)-L-CYSTEINE
		aMap.put("TYY","TYR");                    //  TYY TYR  IMINOQUINONE FORM OF TOPA QUINONONE
		aMap.put("TYZ","ARG");                    //  TYZ ARG  PARA ACETAMIDO BENZOIC ACID
		aMap.put("UMA","ALA");                    //  UMA ALA
		aMap.put("VAD","VAL");                    //  VAD VAL  DEAMINOHYDROXYVALINE
		aMap.put("VAF","VAL");                    //  VAF VAL  METHYLVALINE
		aMap.put("VAL","VAL");                    //  VAL VAL
		aMap.put("VDL","VAL");                    //  VDL VAL  (2R,3R)-2,3-DIAMINOBUTANOIC ACID
		aMap.put("VLL","VAL");                    //  VLL VAL  (2S)-2,3-DIAMINOBUTANOIC ACID
		aMap.put("VME","VAL");                    //  VME VAL  O- METHYLVALINE
		aMap.put("X9Q","ALA");                    //  X9Q ALA
		aMap.put("XX1","LYS");                    //  XX1 LYS  N~6~-7H-PURIN-6-YL-L-LYSINE
		aMap.put("XXY","ALA");                    //  XXY ALA
		aMap.put("XYG","ALA");                    //  XYG ALA
		aMap.put("YCM","CYS");                    //  YCM CYS  S-(2-AMINO-2-OXOETHYL)-L-CYSTEINE
		aMap.put("YOF","TYR");                    //  YOF TYR  3-FLUOROTYROSINE
		modifiedResidues = Collections.unmodifiableMap(aMap);
	}

}
