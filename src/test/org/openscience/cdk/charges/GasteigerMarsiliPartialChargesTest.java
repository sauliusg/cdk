/*
 *  $RCSfile$
 *  $Author$
 *  $Date$
 *   *
 *  Copyright (C) 1997-2007  The Chemistry Development Kit (CDK) project
 *
 *  Contact: cdk-devel@list.sourceforge.net
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation; either version 2.1
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 */
package org.openscience.cdk.charges;


import org.junit.Assert;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.charges.GasteigerMarsiliPartialCharges;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.NewCDKTestCase;

/**
 *  Description of the Class
 *
 * @cdk.module test-charges
 *
 *@author     chhoppe
 *@cdk.created    2004-11-04
 */
public class GasteigerMarsiliPartialChargesTest extends NewCDKTestCase {
	



	
	/**
	 *  A unit test for JUnit with methylenfluoride
	 */
    @Test
    public void testAssignGasteigerMarsiliPartialCharges() throws Exception {
		double [] testResult={0.07915,-0.25264,0.05783,0.05783,0.05783};
		GasteigerMarsiliPartialCharges peoe=new GasteigerMarsiliPartialCharges();
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		IMolecule mol = sp.parseSmiles("CF");
		addExplicitHydrogens(mol);
		peoe.assignGasteigerMarsiliSigmaPartialCharges(mol, true);
		for (int i=0;i<mol.getAtomCount();i++){
			//logger.debug("Charge for atom:"+i+" S:"+mol.getAtomAt(i).getSymbol()+" Charge:"+mol.getAtomAt(i).getCharge());
			Assert.assertEquals(testResult[i],mol.getAtom(i).getCharge(),0.01);
		}
	}
}
