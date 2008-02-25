package org.openscience.cdk.applications.undoredo;

import java.util.ArrayList;

import junit.framework.Test;
import junit.framework.TestSuite;

import org.openscience.cdk.ElectronContainer;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.SingleElectron;
import org.openscience.cdk.applications.undoredo.ConvertToRadicalEdit;
import org.openscience.cdk.templates.MoleculeFactory;
import org.openscience.cdk.CDKTestCase;

/**
 * Junit test for the ConvertToRadicalEdit class
 * 
 * @author tohel
 * @cdk.module test-extra
 */
public class ConvertToRadicalEditTest extends CDKTestCase {

	private ArrayList electronContainerList;

	/**
	 * 
	 */
	public ConvertToRadicalEditTest() {
	}

	/**
	 * @return
	 */
	public static Test suite() {
		return new TestSuite(ConvertToRadicalEditTest.class);
	}

	/*
	 * Test method for
	 * 'org.openscience.cdk.applications.undoredo.ConvertToRadicalEdit.redo()'
	 */
	public void testRedo() {
		Molecule mol = MoleculeFactory.makeAlphaPinene();
		for (int i = 0; i < mol.getAtomCount(); i++) {
			org.openscience.cdk.interfaces.IAtom atom = mol.getAtom(i);
			ElectronContainer electronContainer = new SingleElectron(atom);
			ConvertToRadicalEdit edit = new ConvertToRadicalEdit(mol,
					electronContainer);
			edit.redo();
		}
		int singleElectronContainerCount = mol.getSingleElectronCount();
//		for (int i = 0; i < mol.getElectronContainerCount(); i++) {
//			org.openscience.cdk.interfaces.IElectronContainer container = mol.getElectronContainer(i);
//			if (container.getClass() == SingleElectron.class) {
//				singleElectronContainerCount += 1;
//			}
//		}
		assertTrue(singleElectronContainerCount == mol.getAtomCount());
	}

	/*
	 * Test method for
	 * 'org.openscience.cdk.applications.undoredo.ConvertToRadicalEdit.undo()'
	 */
	public void testUndo() throws Exception {
		Molecule mol = MoleculeFactory.makeAlphaPinene();
		Molecule allRadicalsMol = createAllRadicalsMol(mol);
		for (int i = 0; i < electronContainerList.size(); i++) {
			ElectronContainer container = (ElectronContainer) electronContainerList
					.get(i);
			ConvertToRadicalEdit edit = new ConvertToRadicalEdit(
					allRadicalsMol, container);
			edit.undo();
		}
		int singleElectronContainerCount = allRadicalsMol.getSingleElectronCount();
//		for (int i = 0; i < allRadicalsMol.getElectronContainerCount(); i++) {
//			org.openscience.cdk.interfaces.IElectronContainer container = allRadicalsMol
//					.getElectronContainer(i);
//			if (container instanceof SingleElectron) {
//				singleElectronContainerCount += 1;
//			}
//		}
		assertTrue(singleElectronContainerCount == 0);
	}

	private Molecule createAllRadicalsMol(Molecule mol) throws Exception {
		Molecule allRadicalsMol = (Molecule) mol.clone();
		electronContainerList = new ArrayList();
		for (int i = 0; i < allRadicalsMol.getAtomCount(); i++) {
			org.openscience.cdk.interfaces.IAtom atom = allRadicalsMol.getAtom(i);
			SingleElectron se = new SingleElectron(atom);
			allRadicalsMol.addSingleElectron(se);
			electronContainerList.add(se);
		}
		return allRadicalsMol;

	}

}
