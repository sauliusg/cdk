/* $Revision$ $Author$ $Date$    
 * 
 * Copyright (C) 2004-2007  The Chemistry Development Kit (CDK) project
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
package org.openscience.cdk.modulesuites;

import junit.framework.Test;
import junit.framework.TestSuite;

import org.openscience.cdk.AssociationTest;
import org.openscience.cdk.atomtype.MM2AtomTypeMatcherTest;
import org.openscience.cdk.atomtype.MMFF94AtomTypeMatcherTest;
import org.openscience.cdk.io.VASPReaderTest;
import org.openscience.cdk.tools.GenerateFragmentsTest;

/**
 * TestSuite that runs all the sample tests for experimental classes.
 *
 * @cdk.module test-experimental
 */
public class MexperimentalTests {

    public static Test suite () {
        TestSuite suite= new TestSuite("The cdk.experimental Tests");

        suite.addTest(AssociationTest.suite());
        suite.addTest(VASPReaderTest.suite());
        suite.addTest(GenerateFragmentsTest.suite());

        // from cdk.test.atomtype
        suite.addTest(MMFF94AtomTypeMatcherTest.suite());
        suite.addTest(MM2AtomTypeMatcherTest.suite());
        
        return suite;
    }

}
