package br.embrapa.cnptia.gpbc.plc.descriptors.protein;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import br.embrapa.cnptia.cbi.sdl.core.IAtom;
import br.embrapa.cnptia.cbi.sdl.core.IResidue;
import br.embrapa.cnptia.cbi.sdl.descriptors.contacts.Contact;
import br.embrapa.cnptia.cbi.sdl.descriptors.contacts.ResidueContacts;
import br.embrapa.cnptia.cbi.sdl.descriptors.prototypes.Descriptors;

public class ContactsPerResidue {

	private Map<IResidue, ResidueContacts> contactsPerResidue;
	private Map<IAtom, Set<Contact>> contactsMap;

	public ContactsPerResidue(Descriptors<Contact> contacts) {
		calculate(contacts);
	}

	private void calculate(Descriptors<Contact>contacts) {
		contactsPerResidue = new HashMap<>();
		contactsMap = new HashMap<>();

		for(Contact c : contacts){
			IResidue r = c.getAtom1().getResidue();
			ResidueContacts con = contactsPerResidue.get(r);
			if (con == null){
				con = new ResidueContacts(r);
				contactsPerResidue.put(r, con);
			}
			if (c.isIFR())
				++con.getNumberOfInterfaceContacts()[c.getContactType().ordinal()];
			else
				++con.getNumberOfInternalContacts()[c.getContactType().ordinal()];

			r = c.getAtom2().getResidue();
			con = contactsPerResidue.get(r);
			if (con == null){
				con = new ResidueContacts(r);
				contactsPerResidue.put(r, con);
			}
			if (c.isIFR())
				++con.getNumberOfInterfaceContacts()[c.getContactType().ordinal()];
			else
				++con.getNumberOfInternalContacts()[c.getContactType().ordinal()];

			//Set of contacts per atom
			Set<Contact> set = contactsMap.get(c.getAtom1());
			if (set == null){
				set = new HashSet<>();
				contactsMap.put(c.getAtom1(), set);
			}
			set.add(c);

			set = contactsMap.get(c.getAtom2());
			if (set == null){
				set = new HashSet<>();
				contactsMap.put(c.getAtom2(), set);
			}
			set.add(c);
		}
	}

	public ResidueContacts getContactsSummary(IResidue residue){
		return this.contactsPerResidue.get(residue);
	}

	public Set<Contact> getContacts(IAtom atom) {
		return this.contactsMap.get(atom);
	}
}