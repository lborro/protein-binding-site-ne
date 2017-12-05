package br.embrapa.cnptia.gpbc.plc.config;

public class Config {
	//Curvature
	public static final String SURFRACE_PATH       = "programs/surfrace/surfrace"; 
	public static final String SURFRACE_RADII_PATH = "programs/surfrace/radii.txt";
	
	//Accessibility
	public static final String NACCESS_PATH        = "programs/naccess/naccess"; 
	public static final String NACCESS_RADII_PATH  = "programs/naccess/vdw.radii";
	
	//Electrostatic Potential
	public static final String REDUCE_PATH         = "programs/reduce/reduce";
	public static final String DELPHI_PATH         = "programs/delphi/delphi";
	public static final String REDUCE_DIC_PATH     = "data/reduce_wwPDB_het_dict.txt";
	public static final String PARSE_RADII_PATH    = "data/parse3_red.siz";
	public static final String PARSE_CHARGES_PATH  = "data/parse3_newn_reduce.crg";
	
	//Maximum contacts table
	public static final String MAX_CONTACTS_PATH   = "data/maxcon.txt";
	
	//Directory for temporary files
	public static final String TMP_DIR = "tmp";

	//Time-out for running external programs (Minutes)
	public static final int TIMEOUT = 15; 

	//Cutoff for defining the binding-site nano-environment
	public static final double NANO_ENVIRONMENT_MAX_CUTOFF = 12.0;
}
