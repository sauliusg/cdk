/* $RCSfile$
 * $Author$    
 * $Date$    
 * $Revision$
 * 
 * Copyright (C) 1997-2007  The Chemistry Development Kit (CDK) project
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
 * 
 */
package org.openscience.cdk.debug;

import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.debug.DebugChemObjectBuilder;
import org.openscience.cdk.AminoAcidTest;

/**
 * Checks the functionality of the AtomContainer.
 *
 * @cdk.module test-datadebug
 */
public class DebugAminoAcidTest extends AminoAcidTest {

    @BeforeClass public static void setUp() {
    	AminoAcidTest.builder = DebugChemObjectBuilder.getInstance();
    }

    @Test public void testAminoAcid() {
        super.testAminoAcid();
    }
    
    @Test public void testAddCTerminus_IAtom() {
        super.testAddCTerminus_IAtom();
    }
    
    @Test public void testGetCTerminus() {
        super.testGetCTerminus();
    }

    @Test public void testAddNTerminus_IAtom() {
        super.testAddNTerminus_IAtom();
    }
    
    @Test public void testGetNTerminus() {
        super.testGetNTerminus();
    }
    
    /**
     * Method to test whether the class complies with RFC #9.
     */
    @Test public void testToString() {
    	super.testToString();
    }

    @Test public void testClone() throws Exception {
        super.testClone();
    }
    
}
