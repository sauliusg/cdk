/*---------------------------------------------------------------------------*\
**$Author: saulius $
**$Date: 2017-11-19 11:12:20 +0200 (Sun, 19 Nov 2017) $
**$Revision: 75 $
**$URL: svn://www.crystallography.net/smiles-scripts/trunk/src/cdkcif2smiles.java $
\*---------------------------------------------------------------------------*/

// This CLI Java program reads a CIF (Crystallographic Interchange
// File, http://www.iucr.org/resources/cif) and prints out essential
// data (unic cell parameters, space group, cartesian coordinates).

import org.openscience.cdk.ChemFile;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IChemModel;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.interfaces.ICrystal;
import org.openscience.cdk.io.CIFReader;

import javax.vecmath.Matrix3d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;
import java.io.FileReader;

@SuppressWarnings("ALL")
public class cdkcifrd {

    static final String progName =
        System.getProperty("sun.java.command").split("[ \t]+")[0];
    
    @SuppressWarnings("SingleStatementInBlock")
    public static void main(String[] argv ) {
        String filenames[] = argv;
        if( filenames == null || filenames.length == 0 ) {
            filenames = new String[] { "-" };
        }

        for( String filename : filenames ) {
            try {
                CIFReader reader;
                if (filename.equals("-")) {
                    reader = new CIFReader(System.in);
                } else {
                    reader = new CIFReader(new FileReader(filename));
                }

                ChemFile cifFile = reader.read(new ChemFile());

                assert cifFile != null : "cifFile is null";
                assert cifFile.chemSequences() != null : "cifFile.chemSequences() returned null";
                for ( IChemSequence sequence : cifFile.chemSequences() ) {
                    assert sequence != null : "sequence is null";
                    assert sequence.chemModels() != null : "sequence.chemModels() returned null";
                    for (IChemModel model : sequence.chemModels()) {
                        assert model != null : "model is null";
                        ICrystal crystal = model.getCrystal();
                        assert crystal != null : "Crystal is null";
                        //------------------------------------------------------------
                        // Print out the chemical 'XYZ' format file:
                        System.out.println( crystal.getAtomCount() );
                        System.out.print( filename );
                        System.out.print( "; CELL = " );
                        assert crystal != null;
                        System.out.print( crystal.getA() );
                        System.out.print( " " );
                        System.out.print( crystal.getB() );
                        System.out.print( " " );
                        System.out.print( crystal.getC() );
                        System.out.println();
                        assert crystal.atoms() != null;
                        computeOrthogonalCoordinates( crystal );
                        for( IAtom atom : crystal.atoms() ) {
                            assert atom != null;
                            String symbol = atom.getSymbol();
                            assert symbol != null : "symbol is null";
                            Point3d xyz = atom.getPoint3d();
                            assert xyz != null : "xyz is null";
                            System.out.printf( "%-2s ", symbol );
                            System.out.printf( "%10.6g %10.6g %10.6g\n",
                                               xyz.x, xyz.y, xyz.z );
                        } // for
                    } // for
                } // for
            } // try
            catch(Exception e) {
                System.err.println( progName + ": " + "WARNING, " + e );
            }
        } // for
    }

    public static void computeOrthogonalCoordinates( ICrystal crystal ) {
        assert crystal != null;

        Matrix3d m = changeOfBaseMatrix( crystal.getA(),
                                         crystal.getB(),
                                         crystal.getC() );

        int n = 0;
        for( IAtom atom : crystal.atoms() ) {
            n++;
            Point3d f = atom.getFractionalPoint3d();
            assert f != null : "f is null at atom No " + n;
            Point3d xyz = mpMultiply( m, f );
            assert xyz != null : "xyz is null";
            atom.setPoint3d( xyz );
        }
    }
    
    // Produce a change-of-base matrix given the components of *old*
    // basis vectors in the *new* basis:
    
    public static Matrix3d changeOfBaseMatrix( Vector3d va, Vector3d vb,
                                               Vector3d vc ) {
        Matrix3d m = new Matrix3d();
        m.setColumn(0, va);
        m.setColumn(1, vb);
        m.setColumn(2, vc);
        return m;
    }

    // Multiply a matrix with a vector (given as a Pion3d) from the
    // right. If the 'm' matrix was prodcued by the
    // 'changeOfBaseMatrix()' from this class, and the point 'p'
    // components are in the same *old* basis, then the result wil be
    // the point coordinates in the *new* basis:
    
    public static Point3d mpMultiply( Matrix3d m, Point3d p ) {
        Point3d xyz = new Point3d();
        xyz.setX( m.m00 * p.x + m.m01 * p.y + m.m02 * p.z );
        xyz.setY( m.m10 * p.x + m.m11 * p.y + m.m12 * p.z );
        xyz.setZ( m.m20 * p.x + m.m21 * p.y + m.m22 * p.z );
        return xyz;
    }
    
}
