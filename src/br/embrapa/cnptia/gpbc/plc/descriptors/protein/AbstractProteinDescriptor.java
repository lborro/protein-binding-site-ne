package br.embrapa.cnptia.gpbc.plc.descriptors.protein;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.gpbc.plc.structure.Protein;

public abstract class AbstractProteinDescriptor implements IProteinDescriptor{
	private Protein protein;
	private List<String> names;
	private Map<AminoAcid, List<Double>> valuesPerAminoAcid; 

	protected AbstractProteinDescriptor(Protein protein, String[] descriptorsNames) throws Exception  {
		if(protein == null) throw new NullPointerException("Protein cannot be null");
		this.protein = protein;

		this.names = new ArrayList<>();
		if(descriptorsNames == null) throw new NullPointerException("List of names for this descriptor cannot be null");
		this.names.addAll(Arrays.asList(descriptorsNames));
		this.valuesPerAminoAcid =  new HashMap<>();
	}

	public String[] getDescriptorNames() {
		String[] array = new String[names.size()];
		array = names.toArray(array);
		return array;
	}

	public void addDescriptorValues(AminoAcid aminoAcid, Double[] values) {
		if(values.length != names.size()) new Exception("Number of descriptor values must be :" + names.size());
		List<Double> listOfValues = new ArrayList<>();
		listOfValues.addAll(Arrays.asList(values));
		this.valuesPerAminoAcid.put(aminoAcid, listOfValues);	
	}

	public Double[] getDescriptorValues(AminoAcid aminoAcid) {
		List<Double> values = this.valuesPerAminoAcid.get(aminoAcid);
		if(values != null) {
			Double[] array = new Double[values.size()];
			array = values.toArray(array);
			return array;
		}
		return null;
	}

	public Protein getProtein() {
		return this.protein;
	}
}
