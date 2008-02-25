/* $RCSfile$
 * $Author: miguelrojasch $
 * $Date: 2006-05-11 14:25:07 +0200 (Do, 11 Mai 2006) $
 * $Revision: 6221 $
 *
 *  Copyright (C) 2004-2007  Miguel Rojas <miguel.rojas@uni-koeln.de>
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
package org.openscience.cdk.reaction.type;


import java.util.Iterator;

import junit.framework.Assert;
import junit.framework.Test;
import junit.framework.TestSuite;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IMoleculeSet;
import org.openscience.cdk.interfaces.IReactionSet;
import org.openscience.cdk.reaction.IReactionProcess;
import org.openscience.cdk.reaction.type.ElectronImpactPDBReaction;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.CDKTestCase;

/**
 * TestSuite that runs a test for the ElectronImpactPDBReactionTest.
 *
 * @cdk.module test-reaction
 */
 
public class ElectronImpactPDBReactionTest extends CDKTestCase {
	
	public  ElectronImpactPDBReactionTest() {}
    
	public static Test suite() {
		return new TestSuite(ElectronImpactPDBReactionTest.class);
	}
	/**
	 *  A unit test for JUnit
	 *
	 *@return    Description of the Return Value
	 */
	public void test1_5_Hexen_3_one() throws ClassNotFoundException, CDKException, java.lang.Exception {
		/* ionize >C=C< , set the reactive center*/
		IMolecule reactant = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("C=CCC(=O)CC");
		
		Iterator bonds = reactant.bonds();
		while (bonds.hasNext()){
			IBond bond = (IBond)bonds.next();
			IAtom atom0 = bond.getAtom(0);
			IAtom atom1 = bond.getAtom(1);
			if(bond.getOrder() == IBond.Order.DOUBLE &&
					atom0.getSymbol().equals("C")&&
					atom1.getSymbol().equals("C")){
				bond.setFlag(CDKConstants.REACTIVE_CENTER,true);
			}
		}
		
		IMoleculeSet setOfReactants = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		setOfReactants.addMolecule(reactant);
		
		IReactionProcess type  = new ElectronImpactPDBReaction();
        Object[] params = {Boolean.TRUE};
        type.setParameters(params);
        IReactionSet setOfReactions = type.initiate(setOfReactants, null);
        
        Assert.assertEquals(2, setOfReactions.getReactionCount());

        IMolecule molecule = setOfReactions.getReaction(0).getProducts().getMolecule(0);
        Assert.assertEquals(1, molecule.getAtom(1).getFormalCharge().intValue());
        Assert.assertEquals(1, molecule.getConnectedSingleElectronsCount(molecule.getAtom(0)));

        molecule = setOfReactions.getReaction(1).getProducts().getMolecule(0);
        Assert.assertEquals(1, molecule.getAtom(0).getFormalCharge().intValue());
        Assert.assertEquals(1, molecule.getConnectedSingleElectronsCount(molecule.getAtom(1)));
        
        Assert.assertEquals(3,setOfReactions.getReaction(0).getMappingCount());
		
	}
	/**
	 *  A unit test for JUnit
	 *
	 *@return    Description of the Return Value
	 */
	public void test_Propene() throws ClassNotFoundException, CDKException, java.lang.Exception {
		/* ionize all possible double bonds */
		IMolecule reactant = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("C=CC");
		
		
		IMoleculeSet setOfReactants = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		setOfReactants.addMolecule(reactant);
		
		IReactionProcess type  = new ElectronImpactPDBReaction();
//        Object[] params = {false};
//        type.setParameters(params);
		IReactionSet setOfReactions = type.initiate(setOfReactants, null);
        
		Assert.assertEquals(2, setOfReactions.getReactionCount());

        
        IMolecule molecule = setOfReactions.getReaction(0).getProducts().getMolecule(0);
        Assert.assertEquals(1, molecule.getAtom(1).getFormalCharge().intValue());
        Assert.assertEquals(1, molecule.getConnectedSingleElectronsCount(molecule.getAtom(0)));

        molecule = setOfReactions.getReaction(1).getProducts().getMolecule(0);
        Assert.assertEquals(1, molecule.getAtom(0).getFormalCharge().intValue());
        Assert.assertEquals(1, molecule.getConnectedSingleElectronsCount(molecule.getAtom(1)));
		
	}
	/**
	 *  A unit test for JUnit
	 *
	 *@return    Description of the Return Value
	 */
	public void test2_5_Hexen_3_one() throws ClassNotFoundException, CDKException, java.lang.Exception {
		/* ionize >C=C< , set the reactive center*/
		IMolecule reactant = (new SmilesParser(org.openscience.cdk.DefaultChemObjectBuilder.getInstance())).parseSmiles("C=CCC(=O)CC");
		
			
		IMoleculeSet setOfReactants = DefaultChemObjectBuilder.getInstance().newMoleculeSet();
		setOfReactants.addMolecule(reactant);
		
		IReactionProcess type  = new ElectronImpactPDBReaction();
        Object[] params = {Boolean.FALSE};
        type.setParameters(params);
        IReactionSet setOfReactions = type.initiate(setOfReactants, null);
        
        Assert.assertEquals(2, setOfReactions.getReactionCount());

        IMolecule molecule = setOfReactions.getReaction(0).getProducts().getMolecule(0);
        Assert.assertEquals(1, molecule.getAtom(1).getFormalCharge().intValue());
        Assert.assertEquals(1, molecule.getConnectedSingleElectronsCount(molecule.getAtom(0)));

        molecule = setOfReactions.getReaction(1).getProducts().getMolecule(0);
        Assert.assertEquals(1, molecule.getAtom(0).getFormalCharge().intValue());
        Assert.assertEquals(1, molecule.getConnectedSingleElectronsCount(molecule.getAtom(1)));
        
        Assert.assertEquals(3,setOfReactions.getReaction(0).getMappingCount());
		
	}
}
