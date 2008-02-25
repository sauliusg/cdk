/* $Revision: 7691 $ $Author: egonw $ $Date: 2007-01-11 12:47:48 +0100 (Thu, 11 Jan 2007) $
 * 
 * Copyright (C) 2007  Egon Willighagen <egonw@users.sf.net>
 * 
 * Contact: cdk-devel@lists.sourceforge.net
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA. 
 */
package org.openscience.cdk.tools.manipulator;

import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAminoAcid;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.NewCDKTestCase;
import org.openscience.cdk.tools.manipulator.AminoAcidManipulator;

/**
 * @cdk.module test-standard
 */
public class AminoAcidManipulatorTest extends NewCDKTestCase {
    
	private final static IChemObjectBuilder builder = NoNotificationChemObjectBuilder.getInstance();
	
    public AminoAcidManipulatorTest() {
        super();
    }

    @Test
    public void testRemoveAcidicOxygen_IAminoAcid() throws Exception {
		IAminoAcid glycine = builder.newAminoAcid();
		glycine.add(new SmilesParser(builder).parseSmiles("C(C(=O)O)N"));
		Assert.assertEquals(5, glycine.getAtomCount());
		glycine.addCTerminus(glycine.getAtom(1));
		AminoAcidManipulator.removeAcidicOxygen(glycine);
		Assert.assertEquals(4, glycine.getAtomCount());
	}
	
	/**
	 * @cdk.bug 1646861
	 */
    @Test
    public void testAddAcidicOxygen_IAminoAcid() throws Exception {
		// FIXME: I think this is the proper test, but it currently fails
		IAminoAcid glycine = builder.newAminoAcid();
		glycine.add(new SmilesParser(builder).parseSmiles("C(C=O)N"));
		Assert.assertEquals(4, glycine.getAtomCount());
		glycine.addCTerminus(glycine.getAtom(1));
		AminoAcidManipulator.addAcidicOxygen(glycine);
		Assert.assertEquals(5, glycine.getAtomCount());
	}
	
}


