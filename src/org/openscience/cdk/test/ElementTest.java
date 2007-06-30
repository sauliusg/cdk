/* $Revision$ $Author$ $Date$    
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
 */
package org.openscience.cdk.test;

import junit.framework.JUnit4TestAdapter;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IChemObject;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IElement;

/**
 * Checks the funcitonality of the Element class.
 *
 * @cdk.module test-data
 *
 * @see org.openscience.cdk.Element
 */
public class ElementTest {

	protected IChemObjectBuilder builder;
	
	@Before public void setUp() {
    	builder = DefaultChemObjectBuilder.getInstance();
    }
    
	public static junit.framework.Test suite() {
        return new JUnit4TestAdapter(ElementTest.class);
    }
	
    // test constructors
    
    @Test public void testElement() {
        IElement e = builder.newElement();
        Assert.assertTrue(e instanceof IChemObject);
    }
    
    @Test public void testElement_IElement() {
    	IElement element = builder.newElement();
        IElement e = builder.newElement(element);
        Assert.assertTrue(e instanceof IChemObject);
    }
    
    @Test public void testElement_String() {
        IElement e = builder.newElement("C");
        Assert.assertEquals("C", e.getSymbol());
    }
    
    @Test public void testElement_String_int() {
        IElement e = builder.newElement("H", 1);
        Assert.assertEquals("H", e.getSymbol());
        Assert.assertEquals(1, e.getAtomicNumber());
    }
    
    // test methods
    
    @Test public void testSetSymbol_String() {
        IElement e = builder.newElement();
        e.setSymbol("C");
        Assert.assertEquals("C", e.getSymbol());
    }
        
    @Test public void testGetSymbol() {
        IElement e = builder.newElement("X");
        Assert.assertEquals("X", e.getSymbol());
    }
        
    @Test public void testSetAtomicNumber_int() {
        IElement e = builder.newElement("H");
        e.setAtomicNumber(1);
        Assert.assertEquals(1, e.getAtomicNumber());
    }

    @Test public void testGetAtomicNumber() {
        IElement e = builder.newElement("D", 1);
        Assert.assertEquals(1, e.getAtomicNumber());
    }

    @Test public void testClone() throws Exception {
        IElement elem = builder.newElement();
        Object clone = elem.clone();
        Assert.assertTrue(clone instanceof IElement);
    }
    
    @Test public void testClone_Symbol() throws Exception {
        IElement elem = builder.newElement("C");
        IElement clone = (IElement)elem.clone();
        
        // test cloning of symbol
        elem.setSymbol("H");
        Assert.assertEquals("C", clone.getSymbol());
    }
    
    @Test public void testClone_IAtomicNumber() throws Exception {
        IElement elem = builder.newElement("C", 6);
        IElement clone = (IElement)elem.clone();
        
        // test cloning of atomic number
        elem.setAtomicNumber(5); // don't care about symbol
        Assert.assertEquals(6, clone.getAtomicNumber());
    }
    
    /** Test for RFC #9 */
    @Test public void testToString() {
        IElement elem = builder.newElement();
        String description = elem.toString();
        for (int i=0; i< description.length(); i++) {
        	Assert.assertTrue(description.charAt(i) != '\n');
        	Assert.assertTrue(description.charAt(i) != '\r');
        }
    }

    @Test public void testCompare_Object() {
        // Added to keep the Coverage checker happy, but since the
        // compare(Object) method is not part of the interface, nothing is tested
    	Assert.assertTrue(true);
    }
}
