package br.embrapa.cnptia.gpbc.plc.descriptors.protein;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.descriptors.contacts.Contact;
import br.embrapa.cnptia.cbi.sdl.descriptors.contacts.ContactsCalculator;
import br.embrapa.cnptia.cbi.sdl.descriptors.prototypes.Descriptors;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.gpbc.plc.config.Config;
import br.embrapa.cnptia.gpbc.plc.data.MaxContactsTable;
import br.embrapa.cnptia.gpbc.plc.structure.ProteinLigandComplex;

public class ProteinNanoEnvironment {

	private final ProteinLigandComplex proteinLigandComplex;
	private final MaxContactsTable maxCon;
	private final String epDir;

	private Descriptors<Contact> contacts;
	private ContactsPerResidue conRes; 

	private List<String> descriptorsNames;
	private List<Double> descriptorsValues;
	private List<IProteinDescriptor> descriptors;

	public ProteinNanoEnvironment(ProteinLigandComplex complex, MaxContactsTable maxCon, String epDir) throws Exception {
		this.proteinLigandComplex = complex;
		this.maxCon = maxCon;
		this.epDir = epDir;
		
		if(epDir == null || !(new File(epDir).exists()))
			throw new FileNotFoundException("EP dir does not exist: " + epDir);

		if(complex == null) 
			throw new NullPointerException("ProteinLigandComplex can not be null!");

		if(maxCon == null) 
			throw new NullPointerException("MaxContactsTable can not be null!");
		
		descriptorsNames  = new ArrayList<>();
		descriptorsValues = new ArrayList<>();
		descriptors = new ArrayList<>();

		calculateResidueDescriptors();
	}

	public String[] getDescriptorsNames() {
		String[] array = new String[descriptorsNames.size()];
		array = descriptorsNames.toArray(array);
		return array;
	}

	public Double[] getDescriptorsValues() {
		Double[] array = new Double[descriptorsValues.size()];
		array = descriptorsValues.toArray(array);
		return array;
	}

	public ProteinLigandComplex getProteinLigandComplex() {
		return proteinLigandComplex;
	}

	private void calculateResidueDescriptors() throws Exception {
		ContactsCalculator conCalc = new ContactsCalculator();
		contacts = conCalc.calculate(proteinLigandComplex.getProtein().getStructure(), null);

		conRes = new ContactsPerResidue(contacts);
		descriptors.add(new DensitySponge(proteinLigandComplex.getProtein()));
		descriptors.add(new EnergyDensity(proteinLigandComplex.getProtein(), conRes));
		descriptors.add(new UnusedContacts(proteinLigandComplex.getProtein(), conRes, maxCon));
		descriptors.add(new ElectrostaticPotential(proteinLigandComplex.getProtein(), epDir));

		Curvature curvature = new Curvature(
				proteinLigandComplex.getProtein(),
				Config.SURFRACE_PATH, 
				Config.SURFRACE_RADII_PATH, 
				1, 1.4,"tmp");

		descriptors.add(curvature);


		SolventAccessibleSurfaceArea sasa = new SolventAccessibleSurfaceArea(
				proteinLigandComplex.getProtein(), 
				Config.NACCESS_PATH, 
				Config.NACCESS_RADII_PATH,
				1.4, 
				"tmp"
				);
		descriptors.add(sasa);

		for(IProteinDescriptor descriptor : descriptors) 
			descriptorsNames.addAll(Arrays.asList(descriptor.getDescriptorNames()));

	}

	public void calculate(double maxDist) throws Exception {
		Double[] values = new Double[descriptorsNames.size()];
		Arrays.fill(values, 0d);
		descriptorsValues.addAll(Arrays.asList(values));

		for(IResidue residue : proteinLigandComplex.getLigandNeighborResidues()){
			if(!Utils.isAminoAcid(residue.getName())) continue;
			AminoAcid aa = (AminoAcid) residue;
			double dist = proteinLigandComplex.getDistanceFromAAToLigand(aa);
			if(dist <= maxDist) {
				int idx = 0;
				for(IProteinDescriptor descriptor : descriptors) {
					try {
						Double[] aaDescriptorValues = descriptor.getDescriptorValues(aa);
						for(Double aaValue : aaDescriptorValues) {
							Double aux = descriptorsValues.get(idx);
							descriptorsValues.set(idx++, aux + aaValue*Math.exp(-dist));
						}

					} catch(NullPointerException e) {
						idx += descriptor.getDescriptorNames().length; 
					} 
				} 
			}
		}
	}
}
