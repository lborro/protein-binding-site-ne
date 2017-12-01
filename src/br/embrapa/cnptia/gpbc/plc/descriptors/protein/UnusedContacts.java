package br.embrapa.cnptia.gpbc.plc.descriptors.protein;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import br.embrapa.cnptia.cbi.sdl.core.AminoAcid;
import br.embrapa.cnptia.cbi.sdl.core.Chain;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.core.types.ContactType;
import br.embrapa.cnptia.cbi.sdl.descriptors.contacts.ResidueContacts;
import br.embrapa.cnptia.cbi.sdl.utils.Utils;
import br.embrapa.cnptia.gpbc.plc.data.MaxContactsTable;
import br.embrapa.cnptia.gpbc.plc.structure.Protein;

public class UnusedContacts extends AbstractProteinDescriptor {

	private MaxContactsTable maxConTable;
	private ContactsPerResidue conPerResidue;

	public UnusedContacts(Protein protein, ContactsPerResidue conPerResidue, MaxContactsTable maxCon) throws Exception{
		super(protein, new String[]{
				"Uncon_Hydrophobic", 
				"Uncon_ChargedAttr", 
				"Uncon_ChargedRepu", 
				"Uncon_HB_MM", 
				"Uncon_HB_MS", 
				"Uncon_HB_SS", 
				"Uncon_Aromatic",
				"Uncon_Total"
		});

		this.maxConTable = maxCon;
		this.conPerResidue = conPerResidue;

		calculate();

	}

	private void calculate() {
		List<Integer> allowedContacts = new ArrayList<>();
		allowedContacts.add(ContactType.AROMATIC.ordinal());    //Hydrophobic
		allowedContacts.add(ContactType.CHARGE_ATTR.ordinal()); //Charged Attractive
		allowedContacts.add(ContactType.CHARGE_REPU.ordinal()); //Charged Repulsive
		allowedContacts.add(ContactType.HB_MM.ordinal());       //HB:MM
		allowedContacts.add(ContactType.HB_MS.ordinal());       //HB:MS
		allowedContacts.add(ContactType.HB_SS.ordinal());       //HB:SS
		allowedContacts.add(ContactType.AROMATIC.ordinal());    //Aromatic

		for(Chain chain: this.getProtein().getStructure()) {
			for(IResidue residue : chain) {
				if(!Utils.isAminoAcid(residue.getName())) continue;
				AminoAcid aa = (AminoAcid) residue;
				//calculate unused contacts
				ResidueContacts contactsSummary = conPerResidue.getContactsSummary(residue);
				int intContacts[], ifrContacts[]; 
				if(contactsSummary != null){
					intContacts = contactsSummary.getNumberOfInternalContacts();
					ifrContacts = contactsSummary.getNumberOfInterfaceContacts();	
				} else {
					intContacts = new int[ContactType.NUMBER_OF_PROTEIN_CONTACTS_TYPE];
					ifrContacts = new int[ContactType.NUMBER_OF_PROTEIN_CONTACTS_TYPE];
					Arrays.fill(intContacts, 0);
					Arrays.fill(ifrContacts, 0);
				}
				int max[] = maxConTable.getMaxContacts(aa.getAminoAcidType());

				Double[] unConValues = new Double[allowedContacts.size() + 1];
				double total = 0;

				int idx = 0;
				for(Integer conIdx : allowedContacts) {
					int numberOfUnCon = max[conIdx] - intContacts[conIdx] - ifrContacts[conIdx];
					if(numberOfUnCon <= 0) 
						unConValues[idx] = 0d;
					else { 
						unConValues[idx] = numberOfUnCon*ContactType.getEnergies()[conIdx];
						total += numberOfUnCon*ContactType.getEnergies()[conIdx];
					}
					idx++;
				}
				unConValues[idx] = total;
				this.addDescriptorValues(aa, unConValues);
			}
		}
	}
}