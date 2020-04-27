//hi

import java.awt.BorderLayout;
import java.awt.EventQueue;
import java.awt.Color;
import java.io.*;
import java.nio.FloatBuffer;
import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.text.DecimalFormat;

import com.jogamp.opengl.math.Matrix4;
import org.apache.commons.lang3.math.NumberUtils; //this for checking if some inputs are indeed numeric
import com.jogamp.opengl.math.FloatUtil;

import com.jogamp.opengl.awt.GLJPanel;
import com.jogamp.opengl.util.GLBuffers;
import com.jogamp.common.nio.Buffers;
import com.jogamp.opengl.util.gl2.GLUT;
import com.jogamp.opengl.GL;
import com.jogamp.opengl.GL2;
import com.jogamp.opengl.GLAutoDrawable;
import com.jogamp.opengl.GLEventListener;
import com.jogamp.opengl.GLProfile;
import com.jogamp.opengl.GLCapabilities;
import com.jogamp.opengl.fixedfunc.GLMatrixFunc;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.border.EmptyBorder;
import java.awt.Toolkit;
import javax.swing.JTabbedPane;
import javax.swing.GroupLayout;
import javax.swing.GroupLayout.Alignment;
import javax.swing.JLabel;
import javax.swing.JTextField;
import javax.swing.LayoutStyle.ComponentPlacement;
import javax.swing.SwingConstants;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;

import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Objects;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.DefaultTableModel;
import java.awt.ComponentOrientation;
import javax.swing.border.TitledBorder;
import javax.swing.JButton;
import javax.swing.UIManager;
import javax.swing.JSeparator;
import javax.swing.JSpinner;
import javax.swing.ListSelectionModel;
import javax.swing.SpinnerNumberModel;
import javax.swing.border.EtchedBorder;
import javax.swing.JTree;
import java.awt.Component;
import java.awt.Dimension;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;
import javax.swing.tree.DefaultMutableTreeNode;
import org.eclipse.wb.swing.FocusTraversalOnArray;
import javax.swing.event.TreeSelectionListener;
import javax.swing.event.TreeSelectionEvent;
import javax.swing.ImageIcon;
import javax.swing.border.BevelBorder;
import javax.swing.event.CaretListener;
import javax.swing.event.CaretEvent;
import javax.vecmath.AxisAngle4f;
import javax.vecmath.Vector3f;
import javax.vecmath.Matrix3f;
import java.awt.event.MouseWheelListener;
import java.awt.event.MouseWheelEvent;
import java.util.Vector;
import javax.swing.JComboBox;

/*import com.sun.javafx.geom.Matrix3f;
import com.sun.javafx.geom.Quat4f;
import com.sun.javafx.geom.Vector3f;*/

import static com.jogamp.opengl.fixedfunc.GLMatrixFunc.GL_MODELVIEW_MATRIX;
import com.jogamp.opengl.GL3;
import sun.security.provider.ConfigFile;


public class GUIFRAME extends JFrame
{
    /**
     *
     */
    private static final long serialVersionUID = 1L;
    private JPanel contentPane;
    private JTextField textFieldNUMiterations;
    private JTextField textFieldTimeStep;
    private JTextField textFieldBoxX;
    private JTextField textFieldBoxY;
    private JTextField textFieldBoxZ;
    private JTextField textFieldConfOutFreq;
    private JTextField textFieldRestartOutFreq;
    private JTextField textFieldStatOutFreq;
    private JMenu mnEx;
    private JMenuItem mntmExport;
    private JMenuItem mntmExportrun;
    private JMenuItem mntmExit;
    private JTable tableMolecule;
    private JTable tableInterface;
    private JTextField textFieldInterfaceName;
    private JTextField textFieldInterfaceX;
    private JTextField textFieldInterfaceY;
    private JTextField textFieldInterfaceZ;
    private JTextField textFieldMoleculeName;
    private JTextField textFieldMoleculeCount;
    private JTextField textFieldTransDifCoeff;
    private JTextField textFieldRotDifCoeff;
    private JButton btnAddInterface;
    private JButton btnUpdateInterface;
    private JButton btnDeleteInterface;
    private JButton btnAddMolecule;
    private JButton btnUpdateMolecule;
    private JButton btnDeleteMolecule;
    private JLabel lblMoleculeNameEmptyError;
    private JCheckBox chckbxAnchoredToMembrane;
    private JCheckBox chckbxIsLipid;
    private JCheckBox chckbxIsImplicitLipid;
    private JScrollPane scrollPaneReactions;
    private JTable tableReaction;
    private JButton btnAddReaction;
    private JButton btnUpdateReaction;
    private JButton btnDeleteReaction;
    private JPanel panelInteractionDesigner;
    private JPanel panelDrawInterface;
    // updating reaction variables
    private JLabel lblTheta;
    private JSpinner SpinnerTheta;
    private JLabel lblThetab;
    private JSpinner SpinnerThetab;
    private JLabel lblPhi;
    private JSpinner SpinnerPhi;
    private JLabel lblPhib;
    private JSpinner SpinnerPhib;
    ArrayList <Molecule> MolCla = new ArrayList<Molecule>();
    ArrayList <Reaction> Reax = new ArrayList<Reaction>();
    private JCheckBox chckbxRestart;
    private JLabel lblOmega;
    private JSpinner SpinnerOmega;
    private JTextField textFieldSigmaX;
    private JLabel lblSigmanm;
    private JPanel panelReactionDesigner;
    private JLabel lblNumberOfProducts;
    private JLabel lblNumberOfReactants;
    private JSpinner spinnerNumReactant;
    private JSpinner spinnerNumProduct;
    private JScrollPane scrollPaneReactant1;
    private JScrollPane scrollPaneProduct1;
    private JScrollPane scrollPaneReactant2;
    private JScrollPane scrollPaneProduct2;
    private JCheckBox chckbxReversibleReaction;
    private JLabel lblForwardRate;
    private JLabel lblBackwardRate;
    private JTextField textFieldForwardRate;
    private JTextField textFieldBackwardRate;
    private JSeparator separator;
    private JSeparator separator_1;
    private JTree treeReactant1;
    private JTree treeReactant2;
    private JTree treeProduct1;
    private JTree treeProduct2;
    private JLabel lblErrorReaction;
    public static int[] react1data= new int[3];//initialize react1data=[moleculeindex,interfaceindex,interfacestateindex]
    public static int[] react2data= new int[3];
    public static int[] productdata= new int[6];//if numproduct is 1, check size6, else, break it into 3 and think about indiv molecs (0-2)&(3-5)
    private JCheckBox chckbx2dReaction;
    private JLabel lblX_1;
    private JLabel lblMoleculeRadiusnm;
    private JTextField textFieldRadius;
    private JLabel lblcalculatediffusioncoefficients;
    private JButton btnCalcDiff;
    private Camera camera;
    private Camera cameraReaction;
    private GLUT glut;
    private int sphereDisplayList;  // A display list for drawing one sphere.

    private FloatBuffer sphereVertexBuffer; // Holds the vertex coords, for use with glDrawArrays
    private FloatBuffer sphereNormalBuffer; // Holds the normal vectors, for use with glDrawArrays

    private float[] sphereVertexArray;  // The same data as in sphereVertexBuffer, stored in an array.
    private float[] sphereNormalArray;  // The same data as in sphereNormalBuffer, stored in an array.

    private int vertexVboId;   // identifier for the Vertex Buffer Object to hold the vertex coords
    private int normalVboId;   // identifier for the Vertex Buffer Object to hold the normal vectors
    private JPanel panel_1;
    private JPanel panel_2;

    private JLabel lblNewLabel;
    private JLabel lblNewLabel_1;
    private JPanel panelMolecule;
    private JSeparator separator_2;
    private JTextField textFieldSigmaY;
    private JTextField textFieldSigmaZ;
    private JLabel labelSigmaX;
    private JLabel labelSigmaY;
    private JLabel labelSigmaZ;
    private JButton parminput; // save button - Nomo
    private JComboBox borderTypeX; //added drop boxes - Nomo
    private JComboBox borderTypeY;
    private JComboBox borderTypeZ;
    private JCheckBox reflectiveX;
    private JCheckBox reflectiveY;
    private JCheckBox reflectiveZ;
    private JCheckBox periodicX;
    private JCheckBox periodicY;
    private JCheckBox periodicZ;

    private float scale1 = 1;
    private float scale2 = 1;

    //these variables most stay in scope for the entire reaction.
    rxnCoordSys sys1;
    rxnCoordSys sys2;
    HashMap<String, Float> angles;
    boolean first = true;

    /**
     * Launch the application.
     */
    public static void main(String[] args)
    {
        EventQueue.invokeLater(new Runnable()
        {
            public void run()
            {
                try
                {
                    GUIFRAME frame = new GUIFRAME();
                    frame.setVisible(true);
                } catch (Exception e)
                {
                    e.printStackTrace();
                }
            }
        });
    }

    /**
     * Create the frame.
     */
    public GUIFRAME()
    {
//		GLJPanel drawable = new GLJPanel();
//		drawable.setPreferredSize(new Dimension(600,600));
////		drawable.addGLEventListener(this);
//		setLayout(new BorderLayout());
//		add(drawable, BorderLayout.CENTER);

        initComponents();
        createEvents();
    }
    /////////////////////////////////////////////////////////////////////////////
    //This Method contains all of the code for creating events and
    //initializing components
    /////////////////////////////////////////////////////////////////////////////
    private void initComponents()
    {
        GLProfile glprofile = GLProfile.getDefault();
        GLCapabilities glcapabilities = new GLCapabilities( glprofile );
        GLJPanel gljpanel = new GLJPanel( glcapabilities );
        GLJPanel gljpanelReaction = new GLJPanel( glcapabilities );

        camera = new Camera();
        camera.setLimits(-10,10,-10,10,-10,10);
        camera.installTrackball(gljpanel);

        cameraReaction = new Camera();
        cameraReaction.setLimits(-15,15,-15,15,-15,15);
        cameraReaction.installTrackball(gljpanelReaction);

        gljpanel.addGLEventListener( new GLEventListener()
        {

            @Override
            /**
             * This method will be called when the GLJPanel is first
             * created.  It can be used to initialize the GL context.
             */
            public void init(GLAutoDrawable drawable)
            {

                glut = new GLUT();
                GL2 gl = drawable.getGL().getGL2();
                gl.glEnable(GL2.GL_LIGHTING);
                gl.glEnable(GL2.GL_LIGHT0);
                gl.glEnable(GL2.GL_DEPTH_TEST);
                gl.glShadeModel(GL2.GL_SMOOTH);
                gl.glEnable(GL2.GL_COLOR_MATERIAL);

                /* Create a display list that contains all the OpenGL
                 * commands that are generated when the sphere is drawn.
                 */
                sphereDisplayList = gl.glGenLists(1);
                gl.glNewList(sphereDisplayList, GL2.GL_COMPILE);
                TexturedShapes.uvSphere(gl, 0.4, 32, 16, false);//radius
                gl.glEndList();

                /* Create the data for glDrawArrays, for render modes 2 and 3
                 */
                createSphereArraysAndVBOs(gl);

            }

            /**
             * Draw the molecule in Molecules screen.
             * Draw a color cube of spheres, using
             *  different colors.  Each of the red, green, and blue components
             *  of the color varies along one dimension of the cube.
             * @param drawable
             */
            public void display(GLAutoDrawable drawable)
            {
                DefaultTableModel model = (DefaultTableModel) tableInterface.getModel();
                float x,y,z;

                GL2 gl = drawable.getGL().getGL2();
                gl.glClear(GL2.GL_COLOR_BUFFER_BIT | GL2.GL_DEPTH_BUFFER_BIT);
                gl.glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
                camera.apply(gl);

                gl.glRotatef(-45f, 0.0f, 1.0f, 0.0f );//rotate for better view
//        		DISPLAY COORDINATE AXES
                axes(gl);

                for(int i=0;i<model.getRowCount();i++){
                    x = (float) Double.parseDouble(model.getValueAt(i, 1).toString().trim());
                    y = (float) Double.parseDouble(model.getValueAt(i, 2).toString().trim());
                    z = (float) Double.parseDouble(model.getValueAt(i, 3).toString().trim());
                    if (x > 0 || y > 0 || z > 0) {
                        gl.glPushMatrix();
                        drawMolecVertex(gl, x, y, z, 1.0f, 1.0f, 0.0f);
                        gl.glPopMatrix();
                    }
                }
                gl.glFlush();  // Make sure all commands are sent to graphics card.
                gl.glFinish(); // Wait for all commands to complete, before checking time.
            }

            public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height)
            {
            }

            public void dispose(GLAutoDrawable drawable)
            {
            }

        });

        /**
         * activated when ever the 3D viewport is changed
         */
        gljpanelReaction.addGLEventListener( new GLEventListener()
        {
            @Override
            /**
             * This method will be called when the GLJPanel is first
             * created.  It can be used to initialize the GL context.
             */
            public void init(GLAutoDrawable drawable)
            {

                glut = new GLUT();
                GL2 gl = drawable.getGL().getGL2();
                gl.glEnable(GL2.GL_LIGHTING);
                gl.glEnable(GL2.GL_LIGHT0);
                gl.glEnable(GL2.GL_DEPTH_TEST);
                gl.glShadeModel(GL2.GL_SMOOTH);
                gl.glEnable(GL2.GL_COLOR_MATERIAL);

                /* Create a display list that contains all the OpenGL
                 * commands that are generated when the sphere is drawn.
                 */
                sphereDisplayList = gl.glGenLists(1);
                gl.glNewList(sphereDisplayList, GL2.GL_COMPILE);
                TexturedShapes.uvSphere(gl, 0.5, 32, 16, false);//radius
                gl.glEndList();

                /* Create the data for glDrawArrays, for render modes 2 and 3
                 */
                createSphereArraysAndVBOs(gl);

            }

            /**
             *get the angles from the spinners.
             * @return
             */
            private HashMap getAngles() {
                HashMap<String, Float> angles = new HashMap(5);
                // get angles from spinner fields
                angles.put("theta", 0.0f);
                angles.put("theta2", 0.0f);
                angles.put("phi", 0.0f);
                angles.put("phi2", 0.0f);
                angles.put("omega", 0.0f);

                if(NumberUtils.isCreatable(SpinnerTheta.getValue().toString()))
                {
                    angles.put("theta", (float) Double.parseDouble(SpinnerTheta.getValue().toString()));
                }

                if(NumberUtils.isCreatable(SpinnerTheta.getValue().toString()))
                {
                    angles.put("theta2", (float) Double.parseDouble(SpinnerThetab.getValue().toString()));
                    //convert to RAD
                }

                if(NumberUtils.isCreatable(SpinnerPhi.getValue().toString()))
                {
                    angles.put("phi", (float) Double.parseDouble(SpinnerPhi.getValue().toString()));
                    //convert to RAD

                }

                if(NumberUtils.isCreatable(SpinnerPhi.getValue().toString()))
                {
                    angles.put("phi2", (float) Double.parseDouble(SpinnerPhib.getValue().toString()));
                    //convert to RAD

                }

                if(NumberUtils.isCreatable(SpinnerOmega.getValue().toString()))
                {
                    angles.put("omega", (float) Double.parseDouble(SpinnerOmega.getValue().toString()));
                }

                return angles;
            }


            public FloatBuffer toFloatBuffer(float[] v){
                ByteBuffer buf=ByteBuffer.allocateDirect(v.length * 4);
                buf.order(ByteOrder.nativeOrder());
                FloatBuffer buffer=buf.asFloatBuffer();
                buffer.put(v);
                buffer.position(0);
                return buffer;
            }

            private boolean almostEqual(float a, float b, float eps){
                return Math.abs(a-b)<eps;
            }

            /**
             * Draw a color cube of spheres, using
             * different colors.  Each of the red, green, and blue components
             * of the color varies along one dimension of the cube.
             */
            public void display(GLAutoDrawable drawable)
            {


                GL2 gl = drawable.getGL().getGL2();
                gl.glClear(GL2.GL_COLOR_BUFFER_BIT | GL2.GL_DEPTH_BUFFER_BIT);
                gl.glClearColor(1.0f, 1.0f, 1.0f, 0.0f);
                cameraReaction.apply(gl);
//

                if(react1data[0]>-1
                        &&react1data[1]>-1
                        &&react2data[0]>-1
                        &&react2data[1]>-1
                        &&(Integer) spinnerNumProduct.getValue()==1
                        &&(Integer) spinnerNumReactant.getValue()==2
                        &&digitizeTreeComplex(1)==1)
                {
                    //Get angles

//        			//Load interface data
                    int molecule1ID = react1data[0];
                    int interface1ID = react1data[1];
                    float int1x = (float) (float) MolCla.get(molecule1ID).getinterfaceCoords(interface1ID, 0);
                    float int1y = (float) MolCla.get(molecule1ID).getinterfaceCoords(interface1ID, 1);
                    float int1z = (float) MolCla.get(molecule1ID).getinterfaceCoords(interface1ID, 2);

                    int molecule2ID = react2data[0];
                    int interface2ID = react2data[1];
                    float int2x = (float) MolCla.get(molecule2ID).getinterfaceCoords(interface2ID, 0);
                    float int2y = (float) MolCla.get(molecule2ID).getinterfaceCoords(interface2ID, 1);
                    float int2z = (float) MolCla.get(molecule2ID).getinterfaceCoords(interface2ID, 2);


                    //Get displacments
                    float sigmax = 1.0f;
                    float sigmay = 0.0f;
                    float sigmaz = 0.0f;

                    gl.glPushMatrix();

                    //Set Length of Sigma
                    if(NumberUtils.isCreatable(textFieldSigmaX.getText()))
                    {
                        sigmax = (float) Double.parseDouble(textFieldSigmaX.getText());
                    }

                    //create sigma vector (1,0,0)
                    Vector3f sigma = new Vector3f(sigmax, sigmay, sigmaz);

                    angles = getAngles();
                    float requestedT = (float) Math.toRadians(angles.get("theta"));
                    float requestedP = (float) Math.toRadians(angles.get("phi"));
                    float requestedO = (float) Math.toRadians(angles.get("omega"));
                    float requestedT2 = (float) Math.toRadians(angles.get("theta2"));
                    float requestedP2 = (float) Math.toRadians(angles.get("phi2"));

                    if (almostEqual(requestedO, 0, 1) && almostEqual(requestedT, 0, 1) &&
                            almostEqual(requestedT2, 0, 1) && almostEqual(requestedP, 0, 1) &&
                            almostEqual(requestedP2, 0, 1)) {
                        first = true;
                    }

                    //Intitialize some stuff
                    if (first) {
                        Vector3f BindSite1 = new Vector3f(int1x, int1y, int1z); //vector to bind site
                        Vector3f BindSite2 = new Vector3f(int2x, int2y, int2z); //vector to bind site 2

                        sys1 = new rxnCoordSys(BindSite1, sigma, MolCla.get(molecule1ID).interfacecoords, interface1ID, true);

                        sys2 = new rxnCoordSys(BindSite2, sigma, MolCla.get(molecule1ID).interfacecoords, interface2ID, false);

                        scale1 = BindSite1.length();
                        scale2 = BindSite2.length();


                        //ensuring display matches internal representation
                        SpinnerPhi.setValue(Math.toDegrees(sys1.curPhiAngle()));
                        SpinnerOmega.setValue(Math.toDegrees(sys1.curOmegaAngle(sys2)));
                        SpinnerTheta.setValue(Math.toDegrees(sys1.curThetaAngle()));
                        SpinnerPhib.setValue(Math.toDegrees(sys2.curPhiAngle()));
                        SpinnerThetab.setValue(Math.toDegrees(sys2.curThetaAngle()));

                    } else {
                        int1x = (sys1.scale)*sys1.curVec.x;
                        int1y = (sys1.scale)*sys1.curVec.y;
                        int1z = (sys1.scale)*sys1.curVec.z;

                        int2x = (sys2.scale)*sys2.curVec.x;
                        int2y = (sys2.scale)*sys2.curVec.y;
                        int2z = (sys2.scale)*sys2.curVec.z;

                    }

                    angles = getAngles();

                    //Vector3f omegaAxis = sys1.omegaAxis();
                    //gl.glRotatef(angles.get("omega"), omegaAxis.x, omegaAxis.y, omegaAxis.z );//y [-pi/2,pi/2]


                    requestedT = (float) Math.toRadians(angles.get("theta"));
                    float curT = sys1.curThetaAngle();
                    requestedP = (float) Math.toRadians(angles.get("phi"));
                    float curP = (float) Math.toRadians(sys1.curPhiAngle());
                    requestedO = (float) Math.toRadians(angles.get("omega"));
                    float curO = sys1.curOmegaAngle(sys2);
                    requestedT2 = (float) Math.toRadians(angles.get("theta2"));
                    float curT2 = sys2.curThetaAngle();
                    requestedP2 = (float) Math.toRadians(angles.get("phi2"));
                    float curP2 = sys2.curPhiAngle();

                    ////THETA ROTATION////
                    if (requestedT > Math.toRadians(179.5)) {
                        SpinnerTheta.setValue(179.5);
                    } else if (requestedT < Math.toRadians(-179.5)) {
                        SpinnerTheta.setValue(-179.5);
                    }else if (!almostEqual(requestedT, curT, (float) .01)) {
                        Vector3f thetaAxis1 = sys1.thetaRotationAxis();
                        float delta = (requestedT - curT);
                        boolean neg = false;
                        if (requestedT < 0) {
                            neg = true;
                        }
                        System.out.println(delta);
                        sys1.Tneg = neg;
                        sys1.updateVector(thetaAxis1, (requestedT - curT));
                    }

                    ///OMEGA ROTATION///

                    if (requestedO < Math.toRadians(-179.5)) {
                        SpinnerOmega.setValue(-179.5);
                    } else if (requestedO > Math.toRadians(179.5)) {
                        SpinnerOmega.setValue(179.5);
                    } else if (!almostEqual(requestedO, curO, (float) .01)) {
                        Vector3f omegaAxis1 = sys1.omegaAxis();
                        boolean neg = false;
                        if (requestedO < 0) {
                            neg = true;
                        }
                        sys1.Oneg = neg;
                        sys1.updateVector(omegaAxis1, (requestedO - curO));
                    }

                    ////PHI ROTATION///
                    if (requestedP > Math.toRadians(179.5)) {
                        SpinnerPhi.setValue(179.5);
                    } else if (requestedP < Math.toRadians(-179.5)) {
                        SpinnerPhi.setValue(-179.5);
                    } else {
                        if (requestedP < 0) {
                            sys1.Pneg = true;
                        } else {
                            sys1.Pneg = false;
                        }
                        sys1.phiMod(requestedP * 2);
                    }

                    //compose and perform the rotation (from initial system state)
                    //AxisAngle4f rot = sys1.getSingleRot();
                    //gl.glRotatef((float) Math.toDegrees(rot.angle), rot.x, rot.y, rot.z);//z [-pi,pi)

                    //draw stuff in the right place
                    gl.glTranslatef(-int1x, -int1y, -int1z);
                    paintSphere(gl, int1x, int1y, int1z, 0.0f, 0.0f, 0.0f);
                    drawMolecule1(gl, molecule1ID, sys1,1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f);


                    //display second molecule
                    gl.glPopMatrix();
                    gl.glPushMatrix();



                    //THETA ROTATION 2 ///
                    if (requestedT2 > Math.toRadians(179.5)) {
                        SpinnerThetab.setValue(179.5);
                    } else if (requestedT2 < Math.toRadians(.5)) {
                        SpinnerThetab.setValue(.5);
                    }else if (!almostEqual(requestedT2, curT2, (float) .01)) {
                        Vector3f thetaAxis2 = sys2.thetaRotationAxis();
                        float delta =  (requestedT2 - curT2);
                        boolean neg = false;
                        if (requestedT2 < 0) {
                            neg = true;
                        }
                        sys2.Tneg = neg;
                        sys2.updateVector(thetaAxis2, delta);
                    }

                    ////PHI ROTATION 2////
                    if (requestedP2 > Math.toRadians(179.5)) {
                        SpinnerPhib.setValue(179.5);
                    } else if (requestedP2 < Math.toRadians(-179.5)) {
                        SpinnerPhib.setValue(-179.5);
                    } else {
                        if (requestedP2 < 0) {
                            sys2.Pneg = true;
                        } else {
                            sys2.Pneg = false;
                        }
                        sys2.phiMod(requestedP2 * 2);
                    }

                    //Put int2 onto int1
                    //gl.glTranslatef(int1x, int1y, int1z);

                    //Displace int2 away from int1
                    gl.glTranslatef(sigmax, sigmay, sigmaz);

                    //translate so that interface2 becomes center of mass
                    gl.glTranslatef(-int2x, -int2y, -int2z);
                    paintSphere(gl, int2x, int2y, int2z, 0.0f, 0.0f, 0.0f);
                    drawMolecule1(gl, molecule2ID, sys2, 1.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f);

                    gl.glPopMatrix();

                    //plot interfacial dashed line
                    drawVertex(gl, 0, 0, 0, sigmax, sigmay, sigmaz, 0.0f, 0.0f, 0.0f);

                    //Set all the spinners to their appropriate values.
//					SpinnerPhi.setValue(Math.toDegrees(sys1.curPhiAngle()));
//					SpinnerOmega.setValue(Math.toDegrees(sys1.curOmegaAngle(sys2)));
//					SpinnerTheta.setValue(Math.toDegrees(sys1.curThetaAngle()));
//					SpinnerPhib.setValue(Math.toDegrees(sys2.curPhiAngle()));
//					SpinnerThetab.setValue(Math.toDegrees(sys2.curThetaAngle()));
                    //System.out.println(sys1.curOmegaAngle(sys2));
                    first = false;

                }

                gl.glFlush();  // Make sure all commands are sent to graphics card.
                gl.glFinish(); // Wait for all commands to complete, before checking time.
            }

            public void reshape(GLAutoDrawable drawable, int x, int y, int width, int height)
            {
            }

            public void dispose(GLAutoDrawable drawable)
            {
            }

        });

        setTitle("NERDSS");

        int maxlevel = 3;
        //		initialize for no selection for molecules in reaction
        for(int i=0;i<maxlevel;i++){
            react1data[i] = -1;
            react2data[i] = -1;
        }

        setIconImage(Toolkit.getDefaultToolkit().getImage(GUIFRAME.class.getResource("/resources/FPR.png")));
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setBounds(100, 100, 1023, 580);

        JMenuBar menuBar = new JMenuBar();
        setJMenuBar(menuBar);

        mnEx = new JMenu("File");
        menuBar.add(mnEx);

        mntmExport = new JMenuItem("Export");
        mnEx.add(mntmExport);

        mntmExit = new JMenuItem("Exit");
        mnEx.add(mntmExit);
        contentPane = new JPanel();
        contentPane.setBorder(new EmptyBorder(5, 5, 5, 5));
        contentPane.setLayout(new BorderLayout(0, 0));
        setContentPane(contentPane);

        JTabbedPane tabbedPane = new JTabbedPane(JTabbedPane.TOP);
        contentPane.add(tabbedPane, BorderLayout.CENTER);

        JPanel molecpanel = new JPanel();
        molecpanel.setFocusCycleRoot(true);
        tabbedPane.addTab("Molecules", null, molecpanel, null);

        JScrollPane scrollPaneMolecule = new JScrollPane();

        JScrollPane scrollPaneInterfaceList = new JScrollPane();
        scrollPaneInterfaceList.setBorder(new TitledBorder(UIManager.getBorder("TitledBorder.border"), "Binding Site List", TitledBorder.LEADING, TitledBorder.TOP, null, new Color(0, 0, 0)));

        JPanel panelEditInterfaces = new JPanel();
        panelEditInterfaces.setBorder(new TitledBorder(UIManager.getBorder("TitledBorder.border"), "Binding Sites", TitledBorder.LEADING, TitledBorder.TOP, null, new Color(0, 0, 0)));

        JLabel lblMoleculeName = new JLabel("<html>Molecule <br>Name</html>");

        JLabel lblMoleculeCount = new JLabel("<html>Molecule <br>Count</html>");

        JLabel lblTranslationalDiffusionCoefficient = new JLabel("<html>Translational <br>Diffusion <br>Coeff. (um2/s)</html>");

        JLabel lblTranslationalDiffusionCoefficient_1 = new JLabel("<html>Rotational <br> Diffusion <br>Coeff. (1/us)<html>");

        textFieldMoleculeName = new JTextField();
        textFieldMoleculeName.setText("MolecX");
        textFieldMoleculeName.setColumns(10);

        textFieldMoleculeCount = new JTextField();
        textFieldMoleculeCount.setText("100");
        textFieldMoleculeCount.setColumns(10);

        textFieldTransDifCoeff = new JTextField();
        textFieldTransDifCoeff.setText("0.0");
        textFieldTransDifCoeff.setColumns(10);

        textFieldRotDifCoeff = new JTextField();
        textFieldRotDifCoeff.setText("0.0");
        textFieldRotDifCoeff.setColumns(10);

        chckbxAnchoredToMembrane = new JCheckBox("Is anchored to membrance?");
        chckbxIsLipid = new JCheckBox("Is lipid?");
        chckbxIsImplicitLipid = new JCheckBox("Is implicit lipid?");

        btnAddMolecule = new JButton("Add Molecule");

        btnUpdateMolecule = new JButton("Update Molecule");


        btnDeleteMolecule = new JButton("Delete Molecule");

        lblMoleculeNameEmptyError = new JLabel(" ");
        lblMoleculeNameEmptyError.setForeground(Color.RED);

        separator_1 = new JSeparator();
        separator_1.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));

        lblMoleculeRadiusnm = new JLabel("<html>Molecule <br>Radius (nm)<br>(optional)</html>");

        textFieldRadius = new JTextField();
        textFieldRadius.setText("10.0");
        textFieldRadius.setColumns(10);

        lblcalculatediffusioncoefficients = new JLabel("<html>Calculate <br>Diffusion <br>Coefficients?</html>");

        btnCalcDiff = new JButton("Yes");

        GroupLayout gl_molecpanel = new GroupLayout(molecpanel);
        gl_molecpanel.setHorizontalGroup(
                gl_molecpanel.createParallelGroup(Alignment.TRAILING)
                        .addGroup(gl_molecpanel.createSequentialGroup()
                                .addContainerGap()
                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.TRAILING)
                                        .addComponent(scrollPaneMolecule, Alignment.LEADING, GroupLayout.DEFAULT_SIZE, 972, Short.MAX_VALUE)
                                        .addGroup(Alignment.LEADING, gl_molecpanel.createSequentialGroup()
                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.LEADING, false)
                                                        .addGroup(gl_molecpanel.createSequentialGroup()
                                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.LEADING)
                                                                        .addComponent(textFieldTransDifCoeff, GroupLayout.PREFERRED_SIZE, 73, GroupLayout.PREFERRED_SIZE)
                                                                        .addComponent(lblTranslationalDiffusionCoefficient, GroupLayout.PREFERRED_SIZE, 84, GroupLayout.PREFERRED_SIZE))
                                                                .addGap(2)
                                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.LEADING, false)
                                                                        .addComponent(lblTranslationalDiffusionCoefficient_1, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                                                        .addComponent(textFieldRotDifCoeff, GroupLayout.DEFAULT_SIZE, 74, Short.MAX_VALUE)))
                                                        .addGroup(gl_molecpanel.createSequentialGroup()
                                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.LEADING)
                                                                        .addComponent(textFieldMoleculeName, GroupLayout.PREFERRED_SIZE, 72, GroupLayout.PREFERRED_SIZE)
                                                                        .addComponent(lblMoleculeName)
                                                                        .addComponent(lblMoleculeRadiusnm, GroupLayout.PREFERRED_SIZE, 97, GroupLayout.PREFERRED_SIZE)
                                                                        .addComponent(textFieldRadius, GroupLayout.PREFERRED_SIZE, 73, GroupLayout.PREFERRED_SIZE))
                                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.LEADING)
                                                                        .addComponent(lblMoleculeCount)
                                                                        .addComponent(textFieldMoleculeCount, GroupLayout.PREFERRED_SIZE, 71, GroupLayout.PREFERRED_SIZE)
                                                                        .addGroup(gl_molecpanel.createParallelGroup(Alignment.TRAILING, false)
                                                                                .addComponent(btnCalcDiff, Alignment.LEADING, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                                                                .addComponent(lblcalculatediffusioncoefficients, Alignment.LEADING, GroupLayout.DEFAULT_SIZE, 72, Short.MAX_VALUE))))
                                                        .addGroup(gl_molecpanel.createSequentialGroup()
                                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.LEADING)
                                                                        .addComponent(chckbxIsLipid, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                                                        .addComponent(chckbxIsImplicitLipid, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))))
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addComponent(panelEditInterfaces, GroupLayout.DEFAULT_SIZE, 477, Short.MAX_VALUE)
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addComponent(scrollPaneInterfaceList, GroupLayout.DEFAULT_SIZE, 307, Short.MAX_VALUE))
                                        .addGroup(Alignment.LEADING, gl_molecpanel.createSequentialGroup()
                                                .addComponent(btnAddMolecule)
                                                .addPreferredGap(ComponentPlacement.UNRELATED)
                                                .addComponent(btnUpdateMolecule)
                                                .addPreferredGap(ComponentPlacement.UNRELATED)
                                                .addComponent(btnDeleteMolecule)
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addComponent(lblMoleculeNameEmptyError, GroupLayout.DEFAULT_SIZE, 635, Short.MAX_VALUE))
                                        .addComponent(separator_1, Alignment.LEADING, GroupLayout.DEFAULT_SIZE, 972, Short.MAX_VALUE))
                                .addContainerGap())
        );
        gl_molecpanel.setVerticalGroup(
                gl_molecpanel.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_molecpanel.createSequentialGroup()
                                .addContainerGap()
                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.LEADING)
                                        .addGroup(gl_molecpanel.createSequentialGroup()
                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.BASELINE)
                                                        .addComponent(lblMoleculeName)
                                                        .addComponent(lblMoleculeCount))
                                                .addGap(8)
                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.BASELINE)
                                                        .addComponent(textFieldMoleculeName, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                        .addComponent(textFieldMoleculeCount, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                                                .addGap(18)
                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.BASELINE)
                                                        .addComponent(lblMoleculeRadiusnm)
                                                        .addComponent(lblcalculatediffusioncoefficients))
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.BASELINE)
                                                        .addComponent(textFieldRadius, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                        .addComponent(btnCalcDiff))
                                                .addGap(29)
                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.BASELINE)
                                                        .addComponent(lblTranslationalDiffusionCoefficient)
                                                        .addComponent(lblTranslationalDiffusionCoefficient_1))
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.BASELINE)
                                                        .addComponent(textFieldTransDifCoeff, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                        .addComponent(textFieldRotDifCoeff, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                                                .addGap(18)
                                                .addComponent(chckbxIsLipid)
                                                .addGap(2)
                                                .addComponent(chckbxIsImplicitLipid)
                                                .addGap(1))
                                        .addComponent(scrollPaneInterfaceList, GroupLayout.DEFAULT_SIZE, 284, Short.MAX_VALUE)
                                        .addComponent(panelEditInterfaces, GroupLayout.DEFAULT_SIZE, 284, Short.MAX_VALUE))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(separator_1, GroupLayout.PREFERRED_SIZE, 4, GroupLayout.PREFERRED_SIZE)
                                .addGap(5)
                                .addGroup(gl_molecpanel.createParallelGroup(Alignment.BASELINE)
                                        .addComponent(btnAddMolecule, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(btnUpdateMolecule, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(btnDeleteMolecule, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(lblMoleculeNameEmptyError))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(scrollPaneMolecule, GroupLayout.DEFAULT_SIZE, 132, Short.MAX_VALUE)
                                .addContainerGap())
        );

        JLabel lblInterfaceName = new JLabel("Site Name");

        textFieldInterfaceName = new JTextField();
        textFieldInterfaceName.setText("Site1");
        textFieldInterfaceName.setColumns(10);

        JLabel lblInterfaceCoordinatesRelative = new JLabel("<html>Coordinates <br>Relative<br> to C.o.M. (nm)</html>");
        lblInterfaceCoordinatesRelative.setHorizontalAlignment(SwingConstants.LEFT);

        JLabel lblZ = new JLabel("z:");

        textFieldInterfaceX = new JTextField();
        textFieldInterfaceX.setText("3.0");
        textFieldInterfaceX.setColumns(10);

        textFieldInterfaceY = new JTextField();
        textFieldInterfaceY.setText("4.0");
        textFieldInterfaceY.setColumns(10);

        textFieldInterfaceZ = new JTextField();
        textFieldInterfaceZ.setText("5.0");
        textFieldInterfaceZ.setColumns(10);

        btnAddInterface = new JButton("Add");
        btnUpdateInterface = new JButton("Update");
        btnDeleteInterface = new JButton("Delete");

        panelMolecule = new JPanel(new BorderLayout());
        panelMolecule.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
        panelMolecule.add(gljpanel);///////////////////////////////////////////////////////

        JLabel lblY = new JLabel("y:");

        lblX_1 = new JLabel("x:");

        GroupLayout gl_panelEditInterfaces = new GroupLayout(panelEditInterfaces);
        gl_panelEditInterfaces.setHorizontalGroup(
                gl_panelEditInterfaces.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_panelEditInterfaces.createSequentialGroup()
                                .addGroup(gl_panelEditInterfaces.createParallelGroup(Alignment.LEADING, false)
                                        .addGroup(gl_panelEditInterfaces.createSequentialGroup()
                                                .addGroup(gl_panelEditInterfaces.createParallelGroup(Alignment.LEADING, false)
                                                        .addComponent(lblZ, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                                        .addComponent(lblY, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                                        .addComponent(lblX_1, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addGroup(gl_panelEditInterfaces.createParallelGroup(Alignment.LEADING)
                                                        .addComponent(textFieldInterfaceZ, GroupLayout.PREFERRED_SIZE, 50, GroupLayout.PREFERRED_SIZE)
                                                        .addComponent(textFieldInterfaceY, GroupLayout.PREFERRED_SIZE, 50, GroupLayout.PREFERRED_SIZE)
                                                        .addComponent(textFieldInterfaceX, GroupLayout.PREFERRED_SIZE, 50, GroupLayout.PREFERRED_SIZE)))
                                        .addComponent(lblInterfaceName)
                                        .addComponent(lblInterfaceCoordinatesRelative, GroupLayout.DEFAULT_SIZE, 81, Short.MAX_VALUE)
                                        .addComponent(btnAddInterface, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(btnDeleteInterface, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(btnUpdateInterface, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(textFieldInterfaceName, 0, 0, Short.MAX_VALUE))
                                .addPreferredGap(ComponentPlacement.UNRELATED)
                                .addComponent(panelMolecule, GroupLayout.DEFAULT_SIZE, 188, Short.MAX_VALUE))
        );
        gl_panelEditInterfaces.setVerticalGroup(
                gl_panelEditInterfaces.createParallelGroup(Alignment.TRAILING)
                        .addGroup(gl_panelEditInterfaces.createSequentialGroup()
                                .addComponent(lblInterfaceName)
                                .addGap(1)
                                .addComponent(textFieldInterfaceName, GroupLayout.PREFERRED_SIZE, 19, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(lblInterfaceCoordinatesRelative)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addGroup(gl_panelEditInterfaces.createParallelGroup(Alignment.BASELINE)
                                        .addComponent(lblX_1)
                                        .addComponent(textFieldInterfaceX, GroupLayout.PREFERRED_SIZE, 18, GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addGroup(gl_panelEditInterfaces.createParallelGroup(Alignment.BASELINE)
                                        .addComponent(lblY)
                                        .addComponent(textFieldInterfaceY, GroupLayout.PREFERRED_SIZE, 18, GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addGroup(gl_panelEditInterfaces.createParallelGroup(Alignment.BASELINE)
                                        .addComponent(lblZ)
                                        .addComponent(textFieldInterfaceZ, GroupLayout.PREFERRED_SIZE, 18, GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(ComponentPlacement.RELATED, 11, Short.MAX_VALUE)
                                .addComponent(btnAddInterface)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(btnUpdateInterface)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(btnDeleteInterface))
                        .addComponent(panelMolecule, GroupLayout.DEFAULT_SIZE, 246, Short.MAX_VALUE)
        );
        panelEditInterfaces.setLayout(gl_panelEditInterfaces);

        tableInterface = new JTable();
        tableInterface.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

        tableInterface.setModel(new DefaultTableModel(
                new Object[][] {
                        {"Center of Mass", "0.0", "0.0", "0.0", "0", "-"},
                },
                new String[] {
                        "Name", "<html><center>x<br>(nm)</br></html>", "<html><center>y<br>(nm)</br></html>", "<html><center>z<br>(nm)</br></html>", "<html><center>#Extra<br>States</br></html>", "States"
                }
        ));
        tableInterface.getColumnModel().getColumn(1).setPreferredWidth(15);
        tableInterface.getColumnModel().getColumn(2).setPreferredWidth(15);
        tableInterface.getColumnModel().getColumn(3).setPreferredWidth(15);
        tableInterface.getColumnModel().getColumn(4).setPreferredWidth(40);
        tableInterface.getColumnModel().getColumn(5).setPreferredWidth(40);
        tableInterface.getTableHeader().setPreferredSize(new Dimension(35,35));
        scrollPaneInterfaceList.setViewportView(tableInterface);

        tableMolecule = new JTable();
        tableMolecule.setComponentOrientation(ComponentOrientation.LEFT_TO_RIGHT);
        tableMolecule.setModel(new DefaultTableModel(
                new Object[][] {
                },
                new String[] {
                        "Name", "Molecule Count", "Anchored to Membrane?", "Trans. Diffusion (um2/s)", "Rot. Diffusion (1/us)", "#Binding Sites"
                }
        ) {
            /**
             *
             */
            private static final long serialVersionUID = 1L;
            boolean[] columnEditables = new boolean[] {
                    false, false, false, false, false, false
            };
            public boolean isCellEditable(int row, int column) {
                return columnEditables[column];
            }
        });
        tableMolecule.getColumnModel().getColumn(0).setPreferredWidth(40);
        tableMolecule.getColumnModel().getColumn(1).setPreferredWidth(70);
        tableMolecule.getColumnModel().getColumn(2).setPreferredWidth(120);
        tableMolecule.getColumnModel().getColumn(3).setPreferredWidth(130);
        tableMolecule.getColumnModel().getColumn(4).setPreferredWidth(100);
        tableMolecule.getColumnModel().getColumn(5).setPreferredWidth(45);
        scrollPaneMolecule.setViewportView(tableMolecule);
        molecpanel.setLayout(gl_molecpanel);
        molecpanel.setFocusTraversalPolicy(new FocusTraversalOnArray(new Component[]{textFieldMoleculeName, textFieldMoleculeCount, textFieldTransDifCoeff, textFieldRotDifCoeff, chckbxIsLipid, chckbxIsImplicitLipid, textFieldInterfaceName, textFieldInterfaceX, textFieldInterfaceY, textFieldInterfaceZ, btnAddInterface, btnUpdateInterface, btnDeleteInterface, btnAddMolecule, btnUpdateMolecule, btnDeleteMolecule}));

        JPanel reacpanel = new JPanel();
        tabbedPane.addTab("Reactions", null, reacpanel, null);

        scrollPaneReactions = new JScrollPane();

        btnAddReaction = new JButton("Add Reaction");

        btnUpdateReaction = new JButton("Update Reaction");

        btnDeleteReaction = new JButton("Delete Reaction");

        panelInteractionDesigner = new JPanel();
        panelInteractionDesigner.setEnabled(false);
        panelInteractionDesigner.setBorder(new TitledBorder(UIManager.getBorder("TitledBorder.border"), "Interaction Designer", TitledBorder.LEADING, TitledBorder.TOP, null, new Color(0, 0, 0)));

        panelReactionDesigner = new JPanel();
        panelReactionDesigner.setBorder(new TitledBorder(null, "Reaction Designer", TitledBorder.LEADING, TitledBorder.TOP, null, null));

        lblErrorReaction = new JLabel(" ");
        lblErrorReaction.setAlignmentX(Component.CENTER_ALIGNMENT);
        lblErrorReaction.setForeground(Color.RED);
        lblErrorReaction.setHorizontalAlignment(SwingConstants.LEFT);

        separator_2 = new JSeparator();
        separator_2.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
        GroupLayout gl_reacpanel = new GroupLayout(reacpanel);
        gl_reacpanel.setHorizontalGroup(
                gl_reacpanel.createParallelGroup(Alignment.TRAILING)
                        .addGroup(gl_reacpanel.createSequentialGroup()
                                .addContainerGap()
                                .addGroup(gl_reacpanel.createParallelGroup(Alignment.LEADING)
                                        .addComponent(scrollPaneReactions, GroupLayout.DEFAULT_SIZE, 972, Short.MAX_VALUE)
                                        .addGroup(gl_reacpanel.createSequentialGroup()
                                                .addComponent(panelReactionDesigner, GroupLayout.DEFAULT_SIZE, 456, Short.MAX_VALUE)
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addComponent(panelInteractionDesigner, GroupLayout.DEFAULT_SIZE, 510, Short.MAX_VALUE))
                                        .addGroup(gl_reacpanel.createSequentialGroup()
                                                .addComponent(btnAddReaction)
                                                .addPreferredGap(ComponentPlacement.UNRELATED)
                                                .addComponent(btnUpdateReaction)
                                                .addPreferredGap(ComponentPlacement.UNRELATED)
                                                .addComponent(btnDeleteReaction)
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addComponent(lblErrorReaction, GroupLayout.DEFAULT_SIZE, 629, Short.MAX_VALUE))
                                        .addComponent(separator_2, GroupLayout.DEFAULT_SIZE, 972, Short.MAX_VALUE))
                                .addContainerGap())
        );
        gl_reacpanel.setVerticalGroup(
                gl_reacpanel.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_reacpanel.createSequentialGroup()
                                .addContainerGap()
                                .addGroup(gl_reacpanel.createParallelGroup(Alignment.LEADING)
                                        .addComponent(panelInteractionDesigner, GroupLayout.DEFAULT_SIZE, 279, Short.MAX_VALUE)
                                        .addComponent(panelReactionDesigner, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(separator_2, GroupLayout.PREFERRED_SIZE, 4, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addGroup(gl_reacpanel.createParallelGroup(Alignment.BASELINE, false)
                                        .addComponent(btnAddReaction)
                                        .addComponent(btnUpdateReaction)
                                        .addComponent(btnDeleteReaction)
                                        .addComponent(lblErrorReaction, GroupLayout.PREFERRED_SIZE, 13, GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(scrollPaneReactions, GroupLayout.DEFAULT_SIZE, 136, Short.MAX_VALUE)
                                .addContainerGap())
        );

        lblNumberOfProducts = new JLabel("Number of Products");

        lblNumberOfReactants = new JLabel("Number of Reactants");

        spinnerNumReactant = new JSpinner();

        spinnerNumReactant.setModel(new SpinnerNumberModel(0, 0, 2, 1));
        spinnerNumReactant.setAlignmentY(Component.BOTTOM_ALIGNMENT);
        spinnerNumReactant.setAlignmentX(Component.RIGHT_ALIGNMENT);

        spinnerNumProduct = new JSpinner();
        spinnerNumProduct.setModel(new SpinnerNumberModel(0, 0, 2, 1));
        spinnerNumProduct.setAlignmentY(Component.BOTTOM_ALIGNMENT);
        spinnerNumProduct.setAlignmentX(Component.RIGHT_ALIGNMENT);

        scrollPaneReactant1 = new JScrollPane();

        scrollPaneProduct1 = new JScrollPane();

        scrollPaneReactant2 = new JScrollPane();

        scrollPaneProduct2 = new JScrollPane();

        chckbxReversibleReaction = new JCheckBox("Reversible Reaction?");

        lblForwardRate = new JLabel("On Rate");

        lblBackwardRate = new JLabel("Off Rate");

        textFieldForwardRate = new JTextField();
        textFieldForwardRate.setHorizontalAlignment(SwingConstants.RIGHT);
        textFieldForwardRate.setText("0.0");
        textFieldForwardRate.setMaximumSize(new Dimension(6, 20));
        textFieldForwardRate.setColumns(10);

        textFieldBackwardRate = new JTextField();
        textFieldBackwardRate.setHorizontalAlignment(SwingConstants.RIGHT);
        textFieldBackwardRate.setEnabled(false);
        textFieldBackwardRate.setText("0.0");
        textFieldBackwardRate.setMaximumSize(new Dimension(6, 20));
        textFieldBackwardRate.setColumns(10);

        separator = new JSeparator();

        chckbx2dReaction = new JCheckBox("2D Reaction?");

        GroupLayout gl_panelReactionDesigner = new GroupLayout(panelReactionDesigner);
        gl_panelReactionDesigner.setHorizontalGroup(
                gl_panelReactionDesigner.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_panelReactionDesigner.createSequentialGroup()
                                .addContainerGap()
                                .addGroup(gl_panelReactionDesigner.createParallelGroup(Alignment.TRAILING)
                                        .addComponent(separator, GroupLayout.DEFAULT_SIZE, 500, Short.MAX_VALUE)
                                        .addGroup(gl_panelReactionDesigner.createSequentialGroup()
                                                .addGroup(gl_panelReactionDesigner.createParallelGroup(Alignment.TRAILING)
                                                        .addComponent(scrollPaneProduct1, GroupLayout.DEFAULT_SIZE, 242, Short.MAX_VALUE)
                                                        .addComponent(scrollPaneReactant1, GroupLayout.DEFAULT_SIZE, 242, Short.MAX_VALUE)
                                                        .addComponent(chckbxReversibleReaction, GroupLayout.DEFAULT_SIZE, 242, Short.MAX_VALUE)
                                                        .addGroup(gl_panelReactionDesigner.createSequentialGroup()
                                                                .addComponent(lblNumberOfReactants, GroupLayout.PREFERRED_SIZE, 150, GroupLayout.PREFERRED_SIZE)
                                                                .addPreferredGap(ComponentPlacement.RELATED, 61, Short.MAX_VALUE)
                                                                .addComponent(spinnerNumReactant, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                                                        .addGroup(gl_panelReactionDesigner.createSequentialGroup()
                                                                .addComponent(lblNumberOfProducts, GroupLayout.PREFERRED_SIZE, 150, GroupLayout.PREFERRED_SIZE)
                                                                .addPreferredGap(ComponentPlacement.RELATED, 61, Short.MAX_VALUE)
                                                                .addComponent(spinnerNumProduct, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)))
                                                .addPreferredGap(ComponentPlacement.UNRELATED)
                                                .addGroup(gl_panelReactionDesigner.createParallelGroup(Alignment.TRAILING)
                                                        .addComponent(scrollPaneProduct2, GroupLayout.DEFAULT_SIZE, 248, Short.MAX_VALUE)
                                                        .addComponent(scrollPaneReactant2, GroupLayout.DEFAULT_SIZE, 248, Short.MAX_VALUE)
                                                        .addComponent(chckbx2dReaction, GroupLayout.DEFAULT_SIZE, 248, Short.MAX_VALUE)
                                                        .addGroup(gl_panelReactionDesigner.createSequentialGroup()
                                                                .addComponent(lblBackwardRate, GroupLayout.DEFAULT_SIZE, 175, Short.MAX_VALUE)
                                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                                .addComponent(textFieldBackwardRate, GroupLayout.PREFERRED_SIZE, 69, GroupLayout.PREFERRED_SIZE))
                                                        .addGroup(gl_panelReactionDesigner.createSequentialGroup()
                                                                .addComponent(lblForwardRate, GroupLayout.DEFAULT_SIZE, 177, Short.MAX_VALUE)
                                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                                .addComponent(textFieldForwardRate, GroupLayout.PREFERRED_SIZE, 67, GroupLayout.PREFERRED_SIZE)))))
                                .addContainerGap())
        );
        gl_panelReactionDesigner.setVerticalGroup(
                gl_panelReactionDesigner.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_panelReactionDesigner.createSequentialGroup()
                                .addGroup(gl_panelReactionDesigner.createParallelGroup(Alignment.BASELINE)
                                        .addComponent(lblNumberOfReactants)
                                        .addComponent(textFieldForwardRate, GroupLayout.PREFERRED_SIZE, 18, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(spinnerNumReactant, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(lblForwardRate))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addGroup(gl_panelReactionDesigner.createParallelGroup(Alignment.LEADING)
                                        .addComponent(scrollPaneReactant1, GroupLayout.DEFAULT_SIZE, 81, Short.MAX_VALUE)
                                        .addComponent(scrollPaneReactant2, GroupLayout.DEFAULT_SIZE, 81, Short.MAX_VALUE))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(separator, GroupLayout.PREFERRED_SIZE, 2, GroupLayout.PREFERRED_SIZE)
                                .addGap(8)
                                .addGroup(gl_panelReactionDesigner.createParallelGroup(Alignment.BASELINE)
                                        .addComponent(textFieldBackwardRate, GroupLayout.PREFERRED_SIZE, 18, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(spinnerNumProduct, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(lblBackwardRate)
                                        .addComponent(lblNumberOfProducts))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addGroup(gl_panelReactionDesigner.createParallelGroup(Alignment.LEADING)
                                        .addComponent(scrollPaneProduct1, GroupLayout.DEFAULT_SIZE, 80, Short.MAX_VALUE)
                                        .addComponent(scrollPaneProduct2, GroupLayout.DEFAULT_SIZE, 80, Short.MAX_VALUE))
                                .addGap(4)
                                .addGroup(gl_panelReactionDesigner.createParallelGroup(Alignment.BASELINE)
                                        .addComponent(chckbxReversibleReaction, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(chckbx2dReaction, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
        );

        treeProduct2 = new JTree();
        treeProduct2.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        treeProduct2.setEnabled(false);
        treeProduct2.setModel(new DefaultTreeModel(new DefaultMutableTreeNode("Product2")));
        scrollPaneProduct2.setViewportView(treeProduct2);

        treeProduct1 = new JTree();
        treeProduct1.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        treeProduct1.setEnabled(false);
        treeProduct1.setModel(new DefaultTreeModel(new DefaultMutableTreeNode("Product1")));
        scrollPaneProduct1.setViewportView(treeProduct1);

        treeReactant2 = new JTree();
        treeReactant2.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        treeReactant2.setEnabled(false);
        treeReactant2.setModel(new DefaultTreeModel(new DefaultMutableTreeNode("Reactant2")));
        scrollPaneReactant2.setViewportView(treeReactant2);

        treeReactant1 = new JTree();
        treeReactant1.getSelectionModel().setSelectionMode(TreeSelectionModel.SINGLE_TREE_SELECTION);
        treeReactant1.setEnabled(false);
        treeReactant1.setModel(new DefaultTreeModel(new DefaultMutableTreeNode("Reactant1")));
        scrollPaneReactant1.setViewportView(treeReactant1);

        panelReactionDesigner.setLayout(gl_panelReactionDesigner);

        panelDrawInterface = new JPanel(new BorderLayout());
        panelDrawInterface.setBorder(new EtchedBorder(EtchedBorder.LOWERED, null, null));
        panelDrawInterface.add(gljpanelReaction);/////////////////////////////////////////////////////////////

        //Define rxn angle 1
        lblTheta = new JLabel("Theta Molecule 1 (\u00B0)");
        lblTheta.setEnabled(false);
        lblTheta.setMinimumSize(new Dimension());

        SpinnerTheta = new JSpinner(new SpinnerNumberModel(.5, -179.5, 179.5, .5));
        SpinnerTheta.setAlignmentX(Component.LEFT_ALIGNMENT);
        SpinnerTheta.setEnabled(false);
//		SpinnerTheta.setText("0.0");
        SpinnerTheta.setMaximumSize(new Dimension(6, 20));
//		SpinnerTheta.setColumns(10);

        //Define rxn angle 1b
        lblThetab = new JLabel("Theta Molecule 2 (\u00B0)");
        lblThetab.setEnabled(false);
        lblThetab.setMinimumSize(new Dimension());

        SpinnerThetab = new JSpinner(new SpinnerNumberModel(.5, -179.5, 179.5, .5));
        SpinnerThetab.setAlignmentX(Component.LEFT_ALIGNMENT);
        SpinnerThetab.setEnabled(false);
//		SpinnerThetab.setText("0.0");
        SpinnerThetab.setMaximumSize(new Dimension(6, 20));
//		SpinnerThetab.setColumns(10);

        //Define rxn angle 2
        lblPhi = new JLabel("Phi Molecule 1 (\u00B0)");
        lblPhi.setEnabled(false);

        SpinnerPhi = new JSpinner(new SpinnerNumberModel(.5, -179.5, 179.5, .5));
        SpinnerPhi.setAlignmentX(Component.LEFT_ALIGNMENT);
        SpinnerPhi.setEnabled(false);
//		SpinnerPhi.setText("0.0");
        SpinnerPhi.setMaximumSize(new Dimension(6, 20));

        //Define rxn angle 2b
        lblPhib = new JLabel("Phi Molecule 2 (\u00B0)");
        lblPhib.setEnabled(false);

        SpinnerPhib = new JSpinner(new SpinnerNumberModel(.5, -179.5, 179.5, .5));
        SpinnerPhib.setAlignmentX(Component.LEFT_ALIGNMENT);
        SpinnerPhib.setEnabled(false);
//		SpinnerPhib.setText("0.0");
        SpinnerPhib.setMaximumSize(new Dimension(6, 20));

        //Define rxn angle 3
        lblOmega = new JLabel("Omega (\u00B0)");
        lblOmega.setEnabled(false);

        SpinnerOmega = new JSpinner(new SpinnerNumberModel(.5, -179.5, 179.5, .5));
        SpinnerOmega.setAlignmentX(Component.LEFT_ALIGNMENT);
        SpinnerOmega.setEnabled(false);
//		SpinnerOmega.setText("0.0");
        SpinnerOmega.setMaximumSize(new Dimension(6, 20));
//		SpinnerOmega.setColumns(10);

        //Define rxn x offset
        textFieldSigmaX = new JTextField();
        textFieldSigmaX.setAlignmentX(Component.LEFT_ALIGNMENT);
        textFieldSigmaX.setEnabled(false);
        textFieldSigmaX.setText("1.0");
        textFieldSigmaX.setMaximumSize(new Dimension(6, 20));
        textFieldSigmaX.setColumns(10);

        lblSigmanm = new JLabel("\u03C3 (nm):");
        lblSigmanm.setEnabled(false);

        labelSigmaX = new JLabel("x:");
        labelSigmaX.setSize(new Dimension(14, 14));
        labelSigmaX.setPreferredSize(new Dimension(14, 14));
        labelSigmaX.setMinimumSize(new Dimension(14, 14));
        labelSigmaX.setMaximumSize(new Dimension(14, 14));
        labelSigmaX.setEnabled(false);

        GroupLayout gl_panelInteractionDesigner = new GroupLayout(panelInteractionDesigner);
        gl_panelInteractionDesigner.setHorizontalGroup(
                gl_panelInteractionDesigner.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_panelInteractionDesigner.createSequentialGroup()
                                .addGroup(gl_panelInteractionDesigner.createParallelGroup(Alignment.LEADING, false)
                                        .addComponent(SpinnerOmega, 0, 0, Short.MAX_VALUE)
                                        .addComponent(lblOmega, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(SpinnerPhi, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(SpinnerPhib, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(lblTheta, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(SpinnerTheta, 0, 0, Short.MAX_VALUE)
                                        .addComponent(lblThetab, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(SpinnerThetab, 0, 0, Short.MAX_VALUE)
                                        .addComponent(lblPhi, GroupLayout.DEFAULT_SIZE, 62, Short.MAX_VALUE)
                                        .addComponent(lblPhib, GroupLayout.DEFAULT_SIZE, 62, Short.MAX_VALUE)
                                        .addGroup(gl_panelInteractionDesigner.createSequentialGroup()
                                                .addGroup(gl_panelInteractionDesigner.createParallelGroup(Alignment.LEADING, false)
                                                        .addComponent(labelSigmaX, GroupLayout.PREFERRED_SIZE, 29, GroupLayout.PREFERRED_SIZE))
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addGroup(gl_panelInteractionDesigner.createParallelGroup(Alignment.LEADING)
                                                        .addComponent(textFieldSigmaX, GroupLayout.PREFERRED_SIZE, 40, GroupLayout.PREFERRED_SIZE)))
                                        .addComponent(lblSigmanm, GroupLayout.DEFAULT_SIZE, GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(panelDrawInterface, GroupLayout.DEFAULT_SIZE, 419, Short.MAX_VALUE))
        );
        gl_panelInteractionDesigner.setVerticalGroup(
                gl_panelInteractionDesigner.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_panelInteractionDesigner.createSequentialGroup()
                                .addComponent(lblTheta, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(SpinnerTheta, GroupLayout.PREFERRED_SIZE, 21, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(lblThetab, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(SpinnerThetab, GroupLayout.PREFERRED_SIZE, 21, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(lblPhi)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(SpinnerPhi, GroupLayout.PREFERRED_SIZE, 21, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(lblPhib)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(SpinnerPhib, GroupLayout.PREFERRED_SIZE, 21, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(lblOmega)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(SpinnerOmega, GroupLayout.PREFERRED_SIZE, 21, GroupLayout.PREFERRED_SIZE)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(lblSigmanm)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addGroup(gl_panelInteractionDesigner.createParallelGroup(Alignment.TRAILING)
                                        .addComponent(labelSigmaX, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(textFieldSigmaX, GroupLayout.PREFERRED_SIZE, 21, GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addContainerGap(20, Short.MAX_VALUE))
                        .addComponent(panelDrawInterface, GroupLayout.DEFAULT_SIZE, 256, Short.MAX_VALUE)
        );
        panelInteractionDesigner.setLayout(gl_panelInteractionDesigner);

        tableReaction = new JTable();
        tableReaction.setModel(new DefaultTableModel(
                new Object[][] {
                },
                new String[] {
                        "Reaction", "\u03C3 (nm)",
                        "On Rate", "Units", "Off Rate", "Units", "Theta 1", "Theta 2", "Phi 1",  "Phi 2", "Omega"
                }
        ) {
            /**
             *
             */
            private static final long serialVersionUID = 1L;
            boolean[] columnEditables = new boolean[] {
                    false, false, false, false, false, false, false, false, false
            };
            public boolean isCellEditable(int row, int column) {
                return columnEditables[column];
            }
        });
        tableReaction.getColumnModel().getColumn(0).setPreferredWidth(100);
        tableReaction.getColumnModel().getColumn(1).setPreferredWidth(80);
        tableReaction.getColumnModel().getColumn(1).setMaxWidth(200);
        tableReaction.getColumnModel().getColumn(2).setPreferredWidth(60);
        tableReaction.getColumnModel().getColumn(2).setMaxWidth(200);
        tableReaction.getColumnModel().getColumn(3).setPreferredWidth(45);
        tableReaction.getColumnModel().getColumn(3).setMaxWidth(200);
        tableReaction.getColumnModel().getColumn(4).setPreferredWidth(60);
        tableReaction.getColumnModel().getColumn(4).setMaxWidth(100);
        tableReaction.getColumnModel().getColumn(5).setPreferredWidth(45);
        tableReaction.getColumnModel().getColumn(5).setMaxWidth(200);
        tableReaction.getColumnModel().getColumn(6).setPreferredWidth(55);
        tableReaction.getColumnModel().getColumn(6).setMaxWidth(200);
        tableReaction.getColumnModel().getColumn(7).setPreferredWidth(55);
        tableReaction.getColumnModel().getColumn(7).setMaxWidth(200);
        tableReaction.getColumnModel().getColumn(8).setPreferredWidth(55);
        tableReaction.getColumnModel().getColumn(8).setMaxWidth(200);
        tableReaction.getColumnModel().getColumn(9).setPreferredWidth(55);
        tableReaction.getColumnModel().getColumn(9).setMaxWidth(200);
        tableReaction.getColumnModel().getColumn(10).setPreferredWidth(55);
        tableReaction.getColumnModel().getColumn(10).setMaxWidth(200);
        scrollPaneReactions.setViewportView(tableReaction);
        reacpanel.setLayout(gl_reacpanel);
        reacpanel.setFocusTraversalPolicy(new FocusTraversalOnArray(new Component[]{spinnerNumReactant, spinnerNumProduct, treeProduct1, treeReactant1, textFieldBackwardRate, textFieldForwardRate, treeProduct2, treeReactant2, chckbxReversibleReaction, SpinnerTheta, SpinnerThetab, SpinnerPhi, SpinnerPhib, SpinnerOmega, textFieldSigmaX, tableReaction, btnAddReaction, btnUpdateReaction, btnDeleteReaction}));

        JPanel parmpanel = new JPanel();
        tabbedPane.addTab("Parameters", null, parmpanel, null);

        JLabel lblNUMiterations = new JLabel("<html>Number of <br>Time Steps</html>");

        textFieldNUMiterations = new JTextField();
        textFieldNUMiterations.setColumns(10);

        JLabel lblTimeStep = new JLabel("Time Step (us)");

        textFieldTimeStep = new JTextField();
        textFieldTimeStep.setColumns(10);

        JLabel lblWriteConfigurationEvery = new JLabel("<html>Freq to Print <br>Statistics (in timesteps)</html>"); // come back here to fix labels - Nomo

        textFieldConfOutFreq = new JTextField();
        textFieldConfOutFreq.setColumns(10);

        JLabel lblWriteOutputFrequency = new JLabel("<html> Freq to Print <br>Configuration (in timesteps)</html>");

        textFieldStatOutFreq = new JTextField();
        textFieldStatOutFreq.setColumns(10);

        JLabel lblWriteRestartFrequency = new JLabel("<html> Freq to Print <br>Restart file (in timesteps)</html>");

        textFieldRestartOutFreq = new JTextField();
        textFieldRestartOutFreq.setColumns(10);

        chckbxRestart = new JCheckBox("Is this a restart simulation?");

        panel_1 = new JPanel();
        panel_1.setAlignmentX(Component.LEFT_ALIGNMENT);
        GroupLayout gl_parmpanel = new GroupLayout(parmpanel);
        gl_parmpanel.setHorizontalGroup(
                gl_parmpanel.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_parmpanel.createSequentialGroup()
                                .addContainerGap()
                                .addGroup(gl_parmpanel.createParallelGroup(Alignment.LEADING, false)
                                        .addComponent(textFieldNUMiterations, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(lblNUMiterations)
                                        .addComponent(lblWriteConfigurationEvery, GroupLayout.PREFERRED_SIZE, 155, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(lblWriteOutputFrequency, GroupLayout.PREFERRED_SIZE, 155, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(lblWriteRestartFrequency,GroupLayout.PREFERRED_SIZE, 155, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(textFieldStatOutFreq, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(textFieldConfOutFreq, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(textFieldRestartOutFreq, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(lblTimeStep)
                                        .addComponent(textFieldTimeStep, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                        .addComponent(chckbxRestart, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(panel_1, GroupLayout.PREFERRED_SIZE, 553, GroupLayout.PREFERRED_SIZE)
                                .addGap(272))
        );
        gl_parmpanel.setVerticalGroup(
                gl_parmpanel.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_parmpanel.createSequentialGroup()
                                .addContainerGap()
                                .addGroup(gl_parmpanel.createParallelGroup(Alignment.LEADING)
                                        .addComponent(panel_1, GroupLayout.DEFAULT_SIZE, 460, Short.MAX_VALUE)
                                        .addGroup(gl_parmpanel.createSequentialGroup()
                                                .addComponent(lblTimeStep)
                                                .addGap(4)
                                                .addComponent(textFieldTimeStep, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                .addGap(10)
                                                .addComponent(lblNUMiterations)
                                                .addGap(2)
                                                .addComponent(textFieldNUMiterations, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                .addGap(10)
                                                .addComponent(lblWriteOutputFrequency)
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addComponent(textFieldStatOutFreq, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                .addPreferredGap(ComponentPlacement.UNRELATED)
                                                .addComponent(lblWriteConfigurationEvery)
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addComponent(textFieldConfOutFreq, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                .addGap(10)
                                                .addComponent(lblWriteRestartFrequency)
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addComponent(textFieldRestartOutFreq, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                                                .addPreferredGap(ComponentPlacement.RELATED)
                                                .addComponent(chckbxRestart, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                                )
                                .addContainerGap()
                        )
        );


        parmpanel.setLayout(gl_parmpanel);
        contentPane.setFocusTraversalPolicy(new FocusTraversalOnArray(new Component[]{molecpanel, scrollPaneMolecule, tableMolecule, btnAddMolecule, btnUpdateMolecule, btnDeleteMolecule, lblMoleculeNameEmptyError, lblMoleculeName, tabbedPane, chckbxIsLipid, chckbxIsImplicitLipid,textFieldRotDifCoeff, lblTranslationalDiffusionCoefficient_1, textFieldTransDifCoeff, lblTranslationalDiffusionCoefficient, textFieldMoleculeName, textFieldMoleculeCount, lblMoleculeCount, panelEditInterfaces, textFieldInterfaceName, lblInterfaceName, lblInterfaceCoordinatesRelative, textFieldInterfaceX, textFieldInterfaceY, lblZ, textFieldInterfaceZ, btnAddInterface, btnUpdateInterface, btnDeleteInterface, scrollPaneInterfaceList, tableInterface, separator_1, reacpanel, panelReactionDesigner, separator, lblNumberOfReactants, spinnerNumReactant, lblNumberOfProducts, spinnerNumProduct, scrollPaneProduct1, treeProduct1, scrollPaneReactant1, treeReactant1, lblBackwardRate, textFieldBackwardRate, lblForwardRate, textFieldForwardRate, scrollPaneProduct2, treeProduct2, scrollPaneReactant2, treeReactant2, chckbxReversibleReaction, panelInteractionDesigner, lblSigmanm, lblOmega, lblPhi, lblTheta, SpinnerOmega, SpinnerPhi, SpinnerTheta, SpinnerPhib, SpinnerThetab, textFieldSigmaX, panelDrawInterface, scrollPaneReactions, tableReaction, btnAddReaction, btnUpdateReaction, btnDeleteReaction, parmpanel, textFieldNUMiterations, lblNUMiterations, lblTimeStep, textFieldTimeStep, lblWriteConfigurationEvery, lblWriteOutputFrequency, lblWriteRestartFrequency, textFieldStatOutFreq, textFieldConfOutFreq}));

        JPanel boundpanel = new JPanel(); //added new panel to the screen (added a tab)
        tabbedPane.addTab("Boundary", null, boundpanel, null);
        JLabel lblBoxLengths = new JLabel("Box Dimensions (nm)");

        textFieldBoxX = new JTextField(); //below are defining the textboxes, labels, and buttons that will be on the boundary tab
        textFieldBoxX.setColumns(10);

        textFieldBoxY = new JTextField();
        textFieldBoxY.setColumns(10);

        textFieldBoxZ = new JTextField();
        textFieldBoxZ.setColumns(10);

        JLabel lblX = new JLabel("x:");
        JLabel label = new JLabel("y:");
        JLabel label_1 = new JLabel("z:");

        JLabel lblBorderTypes = new JLabel("Type of boundary");

        JLabel lblbordx = new JLabel("x:");
        JLabel lblbordy = new JLabel("y:");
        JLabel lblbordz = new JLabel("z:");
        reflectiveX = new JCheckBox("Reflective?");
        reflectiveY = new JCheckBox("Reflective?");
        reflectiveZ = new JCheckBox("Reflective?");
        periodicX = new JCheckBox("Periodic?");
        periodicY = new JCheckBox("Periodic?");
        periodicZ = new JCheckBox("Periodic?");

        for (int i=1; i<=4; i++) {
            for (int j=1; j<=4; j++) {
                parminput = new JButton("Create Parameter Input File");
                parminput.setPreferredSize(new Dimension(1000000000, 100000000));
            }
        }


        panel_2 = new JPanel();
        panel_2.setAlignmentX(Component.LEFT_ALIGNMENT);
        GroupLayout gl_boundpanel = new GroupLayout(boundpanel);
        gl_boundpanel.setAutoCreateGaps(true);
        gl_boundpanel.setAutoCreateContainerGaps(true);

        gl_boundpanel.setHorizontalGroup(gl_boundpanel.createSequentialGroup()
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.TRAILING)
                        .addComponent(lblX)
                        .addComponent(label)
                        .addComponent(label_1))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING)
                        .addComponent(lblBoxLengths)
                        .addComponent(textFieldBoxX, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(textFieldBoxY, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(textFieldBoxZ, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.TRAILING)
                        .addComponent(lblbordx)
                        .addComponent(lblbordy)
                        .addComponent(lblbordz))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING)
                        .addComponent(lblBorderTypes)
                        .addComponent(reflectiveX, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(reflectiveY, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(reflectiveZ, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING)
                        .addComponent(periodicY, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(periodicX, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(periodicZ, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addComponent(parminput));

        gl_boundpanel.setVerticalGroup(gl_boundpanel.createSequentialGroup()
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING)
                        .addComponent(lblBoxLengths)
                        .addComponent(lblBorderTypes))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.BASELINE)
                        .addComponent(lblX)
                        .addComponent(textFieldBoxX, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(lblbordx)
                        .addComponent(reflectiveX, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(periodicX, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING)
                        .addComponent(label)
                        .addComponent(textFieldBoxY, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(lblbordy)
                        .addComponent(reflectiveY, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(periodicY, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING)
                        .addComponent(label_1)
                        .addComponent(textFieldBoxZ, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(lblbordz)
                        .addComponent(reflectiveZ, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE)
                        .addComponent(periodicZ, GroupLayout.PREFERRED_SIZE, GroupLayout.DEFAULT_SIZE, GroupLayout.PREFERRED_SIZE))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING))
                .addGroup(gl_boundpanel.createParallelGroup(GroupLayout.Alignment.LEADING))
                .addComponent(parminput));
        //the following line is very necessary for the panel to actually format correctly.

        lblNewLabel = new JLabel("");
        lblNewLabel.setIcon(new ImageIcon(GUIFRAME.class.getResource("/resources/FPRmid.png")));

        lblNewLabel_1 = new JLabel("\u00A9 Osman N Yogurtcu, Spencer Loggia, Nomongo Dorjsuren, Margaret E Johnson, Johns Hopkins University, 2019");
        GroupLayout gl_panel_1 = new GroupLayout(panel_2);
        gl_panel_1.setHorizontalGroup(
                gl_panel_1.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_panel_1.createSequentialGroup()
                                .addContainerGap()
                                .addGroup(gl_panel_1.createParallelGroup(Alignment.LEADING)
                                        .addComponent(lblNewLabel)
                                        .addComponent(lblNewLabel_1))
                                .addContainerGap(276, Short.MAX_VALUE))
        );
        gl_panel_1.setVerticalGroup(
                gl_panel_1.createParallelGroup(Alignment.LEADING)
                        .addGroup(gl_panel_1.createSequentialGroup()
                                .addContainerGap()
                                .addComponent(lblNewLabel)
                                .addPreferredGap(ComponentPlacement.RELATED)
                                .addComponent(lblNewLabel_1)
                                .addContainerGap(262, Short.MAX_VALUE))
        );

        boundpanel.setLayout(gl_boundpanel);
        panel_2.setLayout(gl_panel_1);
        contentPane.setFocusTraversalPolicy(new FocusTraversalOnArray(new Component[]{reflectiveX, reflectiveY, reflectiveZ, periodicX, periodicY, periodicZ, lblBorderTypes, lblbordx, lblbordy, lblbordz, lblX, label, label_1, textFieldBoxX, textFieldBoxY, textFieldBoxZ, lblBoxLengths, scrollPaneMolecule, tableMolecule, btnAddMolecule, btnUpdateMolecule, btnDeleteMolecule, lblMoleculeNameEmptyError, lblMoleculeName, tabbedPane, chckbxIsLipid, chckbxIsImplicitLipid,textFieldRotDifCoeff, lblTranslationalDiffusionCoefficient_1, textFieldTransDifCoeff, lblTranslationalDiffusionCoefficient, textFieldMoleculeName, textFieldMoleculeCount, lblMoleculeCount, panelEditInterfaces, textFieldInterfaceName, lblInterfaceName, lblInterfaceCoordinatesRelative, textFieldInterfaceX, textFieldInterfaceY, lblZ, textFieldInterfaceZ, btnAddInterface, btnUpdateInterface, btnDeleteInterface, scrollPaneInterfaceList, tableInterface, separator_1, reacpanel, panelReactionDesigner, separator, lblNumberOfReactants, spinnerNumReactant, lblNumberOfProducts, spinnerNumProduct, scrollPaneProduct1, treeProduct1, scrollPaneReactant1, treeReactant1, lblBackwardRate, textFieldBackwardRate, lblForwardRate, textFieldForwardRate, scrollPaneProduct2, treeProduct2, scrollPaneReactant2, treeReactant2, chckbxReversibleReaction, panelInteractionDesigner, lblSigmanm, lblOmega, lblPhi, lblTheta, SpinnerOmega, SpinnerPhi, SpinnerTheta, SpinnerPhib, SpinnerThetab, textFieldSigmaX, panelDrawInterface, scrollPaneReactions, tableReaction, btnAddReaction, btnUpdateReaction, btnDeleteReaction, parmpanel, textFieldNUMiterations, lblNUMiterations, lblTimeStep, textFieldTimeStep, lblWriteConfigurationEvery, lblWriteOutputFrequency, lblWriteRestartFrequency, textFieldStatOutFreq, textFieldConfOutFreq}));

        panel_2.setLayout(gl_panel_1);
        boundpanel.setLayout(gl_boundpanel);
        contentPane.setFocusTraversalPolicy(new FocusTraversalOnArray(new Component[]{molecpanel, scrollPaneMolecule, tableMolecule, btnAddMolecule, btnUpdateMolecule, btnDeleteMolecule, lblMoleculeNameEmptyError, lblMoleculeName, tabbedPane, chckbxIsLipid, chckbxIsImplicitLipid, textFieldRotDifCoeff, lblTranslationalDiffusionCoefficient_1, textFieldTransDifCoeff, lblTranslationalDiffusionCoefficient, textFieldMoleculeName, textFieldMoleculeCount, lblMoleculeCount, panelEditInterfaces, textFieldInterfaceName, lblInterfaceName, lblInterfaceCoordinatesRelative, textFieldInterfaceX, textFieldInterfaceY, lblZ, textFieldInterfaceZ, btnAddInterface, btnUpdateInterface, btnDeleteInterface, scrollPaneInterfaceList, tableInterface, separator_1, reacpanel, panelReactionDesigner, separator, lblNumberOfReactants, spinnerNumReactant, lblNumberOfProducts, spinnerNumProduct, scrollPaneProduct1, treeProduct1, scrollPaneReactant1, treeReactant1, lblBackwardRate, textFieldBackwardRate, lblForwardRate, textFieldForwardRate, scrollPaneProduct2, treeProduct2, scrollPaneReactant2, treeReactant2, chckbxReversibleReaction, panelInteractionDesigner, lblSigmanm, lblOmega, lblPhi, lblTheta,SpinnerTheta, SpinnerThetab, SpinnerPhi, SpinnerPhib, SpinnerOmega, textFieldSigmaX, panelDrawInterface, scrollPaneReactions, tableReaction, btnAddReaction, btnUpdateReaction, btnDeleteReaction, parmpanel, lblBoxLengths, lblX, textFieldBoxX, textFieldNUMiterations, lblNUMiterations, lblTimeStep, textFieldTimeStep, label, textFieldBoxY, label_1, textFieldBoxZ, chckbxRestart, lblWriteConfigurationEvery, lblWriteOutputFrequency, textFieldStatOutFreq, textFieldConfOutFreq}));
        ;


    }


    //----------------- for glDrawArrays and Vertex Buffer Object ----------------------------------
    /**
     * Creates the vertex coordinate and normal vectors for a sphere.
     * The data is stored in the FloatBuffers sphereVertexBuffer and
     * sphereNormalBuffer.  In addition, VBOs are created to hold
     * the data and the data is copied from the FloatBuffers into
     * the VBOs.  (Note: The VBOs are used for render mode 4; the
     * FloatBuffers are used for render mode 3.)
     */

    /**
     * Used for drawing rotation handle axis at interface and generating phi rotation axes

     * @return length 3 array (x,y,z) coordinates for rot vector
     */



    private void createSphereArraysAndVBOs(GL2 gl)
    {
        double radius = 0.2;
        int stacks = 16;
        int slices = 32;
        int size = stacks * (slices+1) * 2 * 3;
        sphereVertexBuffer = GLBuffers.newDirectFloatBuffer(size);
        sphereNormalBuffer = GLBuffers.newDirectFloatBuffer(size);
        sphereVertexArray = new float[size];
        sphereNormalArray = new float[size];
        for (int j = 0; j < stacks; j++) {
            double latitude1 = (Math.PI/stacks) * j - Math.PI/2;
            double latitude2 = (Math.PI/stacks) * (j+1) - Math.PI/2;
            double sinLat1 = Math.sin(latitude1);
            double cosLat1 = Math.cos(latitude1);
            double sinLat2 = Math.sin(latitude2);
            double cosLat2 = Math.cos(latitude2);
            for (int i = 0; i <= slices; i++) {
                double longitude = (2*Math.PI/slices) * i;
                double sinLong = Math.sin(longitude);
                double cosLong = Math.cos(longitude);
                double x1 = cosLong * cosLat1;
                double y1 = sinLong * cosLat1;
                double z1 = sinLat1;
                double x2 = cosLong * cosLat2;
                double y2 = sinLong * cosLat2;
                double z2 = sinLat2;
                sphereNormalBuffer.put( (float)x2 );
                sphereNormalBuffer.put( (float)y2 );
                sphereNormalBuffer.put( (float)z2 );
                sphereVertexBuffer.put( (float)(radius*x2) );
                sphereVertexBuffer.put( (float)(radius*y2) );
                sphereVertexBuffer.put( (float)(radius*z2) );
                sphereNormalBuffer.put( (float)x1 );
                sphereNormalBuffer.put( (float)y1 );
                sphereNormalBuffer.put( (float)z1 );
                sphereVertexBuffer.put( (float)(radius*x1) );
                sphereVertexBuffer.put( (float)(radius*y1) );
                sphereVertexBuffer.put( (float)(radius*z1) );
            }
        }
        for (int i = 0; i < size; i++) {
            sphereVertexArray[i] = sphereVertexBuffer.get(i);
            sphereNormalArray[i] = sphereNormalBuffer.get(i);
        }
        sphereVertexBuffer.rewind();
        sphereNormalBuffer.rewind();
        int[] bufferIDs = new int[2];
        gl.glGenBuffers(2, bufferIDs,0);
        vertexVboId = bufferIDs[0];
        normalVboId = bufferIDs[1];
        gl.glBindBuffer(GL2.GL_ARRAY_BUFFER, vertexVboId);
        gl.glBufferData(GL2.GL_ARRAY_BUFFER, size*4, sphereVertexBuffer, GL2.GL_STATIC_DRAW);
        gl.glBindBuffer(GL2.GL_ARRAY_BUFFER, normalVboId);
        gl.glBufferData(GL2.GL_ARRAY_BUFFER, size*4, sphereNormalBuffer, GL2.GL_STATIC_DRAW);
        gl.glBindBuffer(GL.GL_ARRAY_BUFFER, 0);
    }

    private void drawAxis(GL2 gl)
    {
        gl.glPushMatrix();
        glut.glutSolidCylinder(0.2, 4, 6, 10); // Cylinder, radius 0.02, height 1,
        // base at (0,0,0), lying on z-axis.
        gl.glTranslatef(0, 0, 4); // Move the cone to the top of the cylinder
        glut.glutSolidCone(0.45, 0.55, 12, 5); // Cone, radius 0.1, height 0.3,
        // base at (0,0,0), pointing along z-axis.
        gl.glPopMatrix();
    }

    private void axes(GL2 gl)
    {

        /* DRAW COORDINATE AXES */
        /* Draw the x-axis in red. */
        float r = 1.0f;//i/10.0f;
        float g = 0.0f;//j/10.0f;
        float b = 0.0f;//k/10.0f;
        gl.glColor3f(r,g,b);//color
        gl.glPushMatrix();
        gl.glRotatef(90, 0, 1, 0);
        drawAxis(gl);
        gl.glPopMatrix();
        gl.glRasterPos3d(5.5, 0, 0);
        glut.glutBitmapString(GLUT.BITMAP_HELVETICA_18, "X");

        /* Draw the y-axis in green. */
        r = 0.0f;//i/10.0f;
        g = 1.0f;//j/10.0f;
        b = 0.0f;//k/10.0f;
        gl.glColor3f(r,g,b);//color
        gl.glRasterPos3d(0, 5.5, 0);
        glut.glutBitmapString(GLUT.BITMAP_HELVETICA_18, "Y");
        gl.glPushMatrix();
        gl.glRotatef(-90, 1, 0, 0); // drawAxis draws a z-axis; rotate it onto the y-axis.
        drawAxis(gl);
        gl.glPopMatrix();

        /* Draw the z-axis in blue . */
        r = 0.0f;//i/10.0f;
        g = 0.0f;//j/10.0f;
        b = 1.0f;//k/10.0f;
        gl.glColor3f(r,g,b);//color
        gl.glRasterPos3d(0, 0, 5.5);
        glut.glutBitmapString(GLUT.BITMAP_HELVETICA_18, "Z");
        gl.glPushMatrix();
        gl.glRotatef(-90, 0, 0, 1); // drawAxis draws a z-axis; rotate it onto the z-axis.
        drawAxis(gl);
        gl.glPopMatrix();

    }

    private void axesNOTEXT(GL2 gl, rxnCoordSys sys, int moleculeID, int interfaceID)
    {
        Vector3f modPhiNorm = new Vector3f(0,0,0);

        //normRot = phiRotationAxis(moleculeID, interfaceID, zed);

        if (sys == null) {

        } else {
            modPhiNorm = sys.getNormal();
        }
        gl.glPushMatrix();
        drawMolecVertex(gl, 5*modPhiNorm.x, 5*modPhiNorm.y, 5*modPhiNorm.z, 0, 0, 0);
        //gl.glRotatef((float) Math.toDegrees(modAngle), modRot[0], modRot[1], modRot[2]);
        //gl.glRotatef((float) Math.toDegrees(rotAngle), normRot[0], normRot[1], normRot[2]);
        //drawAxis(gl);
        gl.glPopMatrix();

    }

    private void paintSphere(GL2 gl, float x, float y, float z, float r, float g, float b)
    {

        gl.glColor3f(r,g,b);//color
        gl.glPushMatrix();
        gl.glTranslatef(x,y,z); //position
        gl.glCallList(sphereDisplayList);  // Draw by calling a display list.
        gl.glPopMatrix();

    }

    private void drawVertex(GL2 gl, float x0, float y0, float z0, float x, float y, float z, float r, float g, float b)
    {

        gl.glColor3f(r,g,b);//color
        gl.glPushMatrix();
        gl.glEnable(GL2.GL_LINE_STIPPLE);
        int factor = 1;    // Stippling factor
        short pattern = 0x5555;  // Stipple pattern
        gl.glLineStipple(factor,pattern);
        gl.glLineWidth(4.0f);
        gl.glBegin( GL2.GL_LINES );
        gl.glVertex3f(x0, y0, z0);
        gl.glVertex3f(x, y, z);
        gl.glEnd();
        gl.glPopMatrix();

    }

    private void drawMolecVertex(GL2 gl, float x, float y, float z, float r, float g, float b)
    {

        gl.glColor3f(r,g,b);//color

        float angle = (float) (180/Math.PI*Math.acos(z/Math.sqrt(x*x+y*y+z*z)));
        gl.glRotatef(angle,-y,x, 0.0f);
        glut.glutSolidCylinder(0.25, Math.sqrt(x*x+y*y+z*z), 6, 10); // Cylinder, radius 0.02, height 1,

    }

    private void drawSecondaryMolecVertex(GL2 gl, float phioffset, Vector3f phioffsetAxis, float x, float y, float z, float r, float g, float b)
    {

        gl.glColor3f(r,g,b);//color

        float angle = (float) (180/Math.PI*Math.acos(z/Math.sqrt(x*x+y*y+z*z)));
        gl.glRotatef(-1 * (float)Math.toDegrees(phioffset), phioffsetAxis.x, phioffsetAxis.y, phioffsetAxis.z);
        gl.glRotatef(angle,-y,x, 0.0f);

        glut.glutSolidCylinder(0.25, Math.sqrt(x*x+y*y+z*z), 6, 10); // Cylinder, radius 0.02, height 1,

    }

    private void drawMolecule1(GL2 gl, int moleculeID, rxnCoordSys sys, float rVertex, float gVertex, float bVertex, float rInterface, float gInterface, float bInterface)
    {

        float x, y, z;
        //		CENTER OF MASS
        paintSphere(gl, 0.0f, 0.0f, 0.0f, 1.0f, 1.0f, 1.0f);

        int adjust = 0;
        for(int i=0; i < MolCla.get(moleculeID).getnumInterface(); i++){
            if (i==sys.interfaceNum) {
                x = (sys.scale) * sys.curVec.x;
                y = (sys.scale) * sys.curVec.y;
                z = (sys.scale) * sys.curVec.z;

                gl.glPushMatrix();
                paintSphere(gl, x, y, z, rInterface, gInterface, bInterface);
                drawMolecVertex(gl, x, y, z, rVertex, gVertex, bVertex);
                gl.glPopMatrix();

                gl.glPushMatrix();
                //		DISPLAY normal axis
                axesNOTEXT(gl, sys, moleculeID, i);
                //axes(gl);
                gl.glPopMatrix();
                adjust ++;
            } else {

                try {
                    x = sys.secondaryScales.get(i - adjust) * sys.cur_secondaryInterfaces.get(i - adjust).x;
                    y = sys.secondaryScales.get(i - adjust) * sys.cur_secondaryInterfaces.get(i - adjust).y;
                    z = sys.secondaryScales.get(i - adjust) * sys.cur_secondaryInterfaces.get(i - adjust).z;
                    float phi_offset = sys.phiAngle / 2;
                    Vector3f phi_axis = sys.curVec;

                    gl.glPushMatrix();
                    drawSecondaryMolecVertex(gl, phi_offset, phi_axis, x, y, z, rVertex, gVertex, bVertex);
                    gl.glPopMatrix();
                } catch (IndexOutOfBoundsException e){
                    this.sys1 = null;
                    this.sys2 = null;
                }
            }

        }



    }



    public void cleanReactionTab()
    {

        DefaultTreeModel modelreactant1 = (DefaultTreeModel) treeReactant1.getModel();
        DefaultTreeModel modelreactant2 = (DefaultTreeModel) treeReactant2.getModel();

        spinnerNumReactant.setValue(0);
        spinnerNumProduct.setValue(0);
        changerateunits(0,0);
        textFieldBackwardRate.setText("0.0");
        textFieldForwardRate.setText("0.0");
        chckbxReversibleReaction.setSelected(false);
        chckbx2dReaction.setSelected(false);

        //set initial angle values for new reaction
        SpinnerTheta.setValue(.5);
        SpinnerThetab.setValue(.5);
        SpinnerPhi.setValue(.5);
        SpinnerPhib.setValue(.5);
        SpinnerOmega.setValue(.5);
        textFieldSigmaX.setText("1.0");

        fillinNodes(MolCla,modelreactant1);
        fillinNodes(MolCla,modelreactant2);
        spinnerNumReactant.setValue(0);
        spinnerNumProduct.setValue(0);
        treeReactant1.clearSelection();
        treeReactant2.clearSelection();
        treeProduct1.clearSelection();
        treeProduct2.clearSelection();
    }

    private void expandAllNodes(JTree tree, int startingIndex, int rowCount)
    {
        for(int i=startingIndex;i<rowCount;++i){
            tree.expandRow(i);
        }

        if(tree.getRowCount()!=rowCount){
            expandAllNodes(tree, rowCount, tree.getRowCount());
        }
    }

    public int getRowTree(int[] vec)
    {//find the rownumber of the molecule on the reaction tree

        int row=0;
        int molID = vec[0];
        int intID = vec[1];
        int staID = vec[2];

        for(int i = 0; i<molID+1; i++){

            if(i==molID){
                row = row + 1;//for the molecule itself

                for(int j=0;j<MolCla.get(molID).getnumInterface();j++){

                    if(j==intID){
                        row = row + 1;//for the interface
                        row = row + staID + 1;
                        break;
                    }else if(intID==-1){
                        break;
                    }else{
                        row = row + 1;//for the interface
                        row = row + MolCla.get(molID).getnumInterfaceStates(j);
                    }
                }
            }else{
                row = row + MolCla.get(i).getBulkNum();
            }

        }

        return row;
    }

    public int checkConsistentSigma()
    {

        int consflag = 1;

        if(NumberUtils.isCreatable(textFieldSigmaX.getText())) {

        }else{
            consflag = 0;
            lblErrorReaction.setText("Can't add. Sigma missing.");
        }

        return consflag;
    }

    public int checkSphere()
    {

        int sphereflag = 0;
        DefaultTableModel model = (DefaultTableModel) tableInterface.getModel();

        for(int i=0;i<model.getRowCount();i++){
            if((Objects.equals(Double.parseDouble(model.getValueAt(i, 1).toString().trim()), 0.0))
                    &&Objects.equals(Double.parseDouble(model.getValueAt(i, 2).toString().trim()), 0.0)
                    &&Objects.equals(Double.parseDouble(model.getValueAt(i, 3).toString().trim()), 0.0))
            {
                sphereflag = 1;
            }else{
                sphereflag = 0;
            }
        }

        return sphereflag;
    }

    public double calcTransDif()
    {

        double transDif = 0.0, r;

        if(NumberUtils.isCreatable(textFieldRadius.getText()))
        {
            r = Double.parseDouble(textFieldRadius.getText().toString().trim());
            transDif=1.3806E-23*298/(6*3.14*0.00089*r*1E-9)*1E12;
        }else{
            transDif = 0.0;
            lblMoleculeNameEmptyError.setText("Can't calculate. Radius is zero.");
        }

        return transDif;
    }

    public double calcRotDif()
    {

        double rotDif = 0.0, r;

        if(NumberUtils.isCreatable(textFieldRadius.getText()))
        {
            r = Double.parseDouble(textFieldRadius.getText().toString().trim());
            rotDif=1E-6*1.3806E-23*298/(8*3.14*0.00089*r*r*r*1E-27);
        }else{
            rotDif = 0.0;
            lblMoleculeNameEmptyError.setText("Can't calculate. Radius is zero.");
        }

        return rotDif;
    }

    public int checkConsistentAngles()
    {

        int consflag = 1;

        //use some dumb method to check weather user input for angle is a number

        if(NumberUtils.isCreatable(SpinnerTheta.getValue().toString()) && NumberUtils.isCreatable(SpinnerThetab.getValue().toString()))
        {
            if(NumberUtils.isCreatable(SpinnerPhi.getValue().toString()) && NumberUtils.isCreatable(SpinnerPhib.getValue().toString()))
            {
                if(NumberUtils.isCreatable(SpinnerOmega.getValue().toString()))
                {

                }else{
                    consflag = 0;
                    lblErrorReaction.setText("Can't add. Angle3 missing.");
                }
            }else{
                consflag = 0;
                lblErrorReaction.setText("Can't add. Angle2 missing.");
            }
        }else{
            consflag = 0;
            lblErrorReaction.setText("Can't add. Angle1 missing.");
        }

        return consflag;
    }

    public int checkConsistentReactionDimension()
    {

        int consflag = 1;

        if(chckbx2dReaction.isSelected())
        {//this is a 2D reaction, make sure all the molecules selected are in 2D

            int prod1 = productdata[0];
            int prod2 = productdata[3];
            int reac1 = react1data[0];
            int reac2 = react2data[0];

            if(prod1!=-1 && !MolCla.get(prod1).getmolOnMembrane()){
                consflag = 0;
                lblErrorReaction.setText("Make sure product1 is on the membrane");
            }else if(prod2!=-1 && !MolCla.get(prod2).getmolOnMembrane()){
                consflag = 0;
                lblErrorReaction.setText("Make sure product2 is on the membrane");
            }else if(reac1!=-1 && !MolCla.get(reac1).getmolOnMembrane()){
                consflag = 0;
                lblErrorReaction.setText("Make sure reactant1 is on the membrane");
            }else if(reac2!=-1 && !MolCla.get(reac2).getmolOnMembrane()){
                consflag = 0;
                lblErrorReaction.setText("Make sure reactant2 is on the membrane");
            }

        }

        return consflag;
    }

    public int checkConsistentReactionSelection()
    {

        int consflag = 1;
        int npro, nrea;
        npro = (Integer) spinnerNumProduct.getValue(); //spinner product number
        nrea = (Integer) spinnerNumReactant.getValue(); //spinner reactant number

        DefaultTreeModel modelproduct1 = (DefaultTreeModel) treeProduct1.getModel();
        DefaultTreeModel modelproduct2 = (DefaultTreeModel) treeProduct2.getModel();

        int levelreactant1 = digitizeTree(modelproduct1, modelproduct2, 1, false); //0,1,2,3
        int levelreactant2 = digitizeTree(modelproduct1, modelproduct2, 2, false);	//0,1,2,3

        int levelproduct1 = digitizeTreeComplex(1); //0,1,2,3
        int levelproduct2 = digitizeTreeComplex(2); //0,1,2,3

        if(nrea==0){
            if(npro==0){//not possible
                consflag = 0;
                lblErrorReaction.setText("Can't add. No molecule present.");
            }else if(npro==1){//make sure one product selected, birth reaction
                if(levelproduct1<1){//not even one molecule selected
                    lblErrorReaction.setText("Can't add. Please select a product.");
                    consflag = 0;
                }
            }else if(npro==2){//make sure two products selected
                if(levelproduct1<1 || levelproduct2<1){//not even one molecule selected
                    lblErrorReaction.setText("Can't add. Please select two products.");
                    consflag = 0;
                }
            }
        }else if(nrea==1){
            if(levelreactant1<1){
                lblErrorReaction.setText("Can't add. Please select a reactant.");
                consflag = 0;
            }else{
                if(npro==0){//make sure one reactant selected, decay reaction

                }else if(npro==1){//make sure one reactant selected, one product, conversion reaction
                    if(levelproduct1<1){//not even one molecule selected
                        lblErrorReaction.setText("Can't add. Please select a product.");
                        consflag = 0;
                    }
                }else if(npro==2){//collapse reaction A->B+C
                    if(levelproduct1<1 || levelproduct2<1){//not even one molecule selected
                        lblErrorReaction.setText("Can't add. Please select two products.");
                        consflag = 0;
                    }
                    //ALSO!!!! make sure A is the actual "sum" of B and C molecules
                    lblErrorReaction.setText("Can't add. Products can't constitue the reactant.");
                    consflag = 0;
                }
            }
        }else if(nrea==2){
            if(levelreactant1<1 || levelreactant2<1){//not even one molecule selected
                lblErrorReaction.setText("Can't add. Please select two reactants.");
                consflag = 0;
            }else{
                if(npro==0){
                    consflag = 0;
                }else if(npro==1){
                    if(levelreactant1<2 || levelreactant2<2){//not even one molecule selected
                        lblErrorReaction.setText("Can't add. Please select two binding interfaces.");
                        consflag = 0;
                    }
                }else if(npro==2){
                    if(levelproduct1<1 || levelproduct2<1){//not even one molecule selected
                        lblErrorReaction.setText("Can't add. Please select two products.");
                        consflag = 0;
                    } else if (chckbxReversibleReaction.isSelected()){
                        lblErrorReaction.setText("Can't add. This reaction cannot be reversible.");
                        consflag = 0;
                    }
                }
            }
        }

        return consflag;
    }

    public int digitizeTreeComplex(int activeTreeId)
    {

        int level=0, maxlevel=3;
        int[] productdataT = new int[6];
        DefaultMutableTreeNode node = null;
        for(int i=0;i<6;i++){
            productdataT[i]=-1;//-1 if NULL, initialization
        }

        if(activeTreeId==1){
            node = (DefaultMutableTreeNode) treeProduct1.getLastSelectedPathComponent();
        }else if(activeTreeId==2){
            node = (DefaultMutableTreeNode) treeProduct2.getLastSelectedPathComponent();
        }

        /* if nothing is selected */
        if (node != null){
            level = node.getLevel();//0 if root, 1 if molecule, 2 if interface, 3 if interface state (levels 0&1 are illegal for reactions, 1 would be fine if this is a unimolecular reaction)

            switch(level)
            {
                case 0: //root
                    if((Integer) spinnerNumProduct.getValue()>0){
                        lblErrorReaction.setText("Please select at least one molecule.");
                    }
                    break;
                case 1: //up to molecule
                    productdataT[0] = node.getParent().getIndex(node);//index of the molecule
                    break;
                case 2: //up to interface
                    productdataT[1]=node.getParent().getIndex(node);
                    productdataT[0]=(node.getParent()).getParent().getIndex(node.getParent());
                    break;
                case 3: //up to state
                    productdataT[2]=node.getParent().getIndex(node);//index of the state
                    productdataT[1]=(node.getParent()).getParent().getIndex(node.getParent());//index of the interface
                    productdataT[0]=((node.getParent()).getParent()).getParent().getIndex((node.getParent()).getParent());//index of the molecule
                    break;
            }
        }
        if(activeTreeId==1){
            for(int i=0;i<maxlevel;i++){
                productdata[i]=productdataT[i];
            }
        }else if(activeTreeId==2){
            for(int i=0;i<maxlevel;i++){
                productdata[i+maxlevel]=productdataT[i];
            }
        }

        return level;
    }

    public int digitizeTree(DefaultTreeModel p1, DefaultTreeModel p2, int activeTreeId, boolean generateProducts)
    {

        int level=0;
        int[] reactdata = new int[3];
        DefaultMutableTreeNode node = null;
        for(int i=0;i<3;i++){
            reactdata[i]=-1;//-1 if NULL, initialization
        }

        if(activeTreeId==1){
            node = (DefaultMutableTreeNode) treeReactant1.getLastSelectedPathComponent();
        }else if(activeTreeId==2){
            node = (DefaultMutableTreeNode) treeReactant2.getLastSelectedPathComponent();
        }

        /* if nothing is selected */
        if (node != null){

            level = node.getLevel();//0 if root, 1 if molecule, 2 if interface, 3 if interface state (levels 0&1 are illegal for reactions, 1 would be fine if this is a unimolecular reaction)

            switch(level)
            {
                case 0: //root
                    if((Integer) spinnerNumReactant.getValue()>0){
                        lblErrorReaction.setText("Please select at least one molecule.");
                    }
                    flushTree(p1);
                    flushTree(p2);
                    break;
                case 1: //up to molecule
                    reactdata[0] = node.getParent().getIndex(node);//index of the molecule
                    break;
                case 2: //up to interface
                    reactdata[1]=node.getParent().getIndex(node);
                    reactdata[0]=(node.getParent()).getParent().getIndex(node.getParent());
                    break;
                case 3: //up to state
                    reactdata[2]=node.getParent().getIndex(node);//index of the state
                    reactdata[1]=(node.getParent()).getParent().getIndex(node.getParent());//index of the interface
                    reactdata[0]=((node.getParent()).getParent()).getParent().getIndex((node.getParent()).getParent());//index of the molecule
                    break;
            }
        }
        if(activeTreeId==1){
            react1data=reactdata;
        }else if(activeTreeId==2){
            react2data=reactdata;
        }

        if(level>0 && generateProducts==true){
            generateProductforTrees(p1, p2);
        }

        return level;
    }

    public void activateInterfaceDesigner(int activate)
    {
        if(activate == 0){
            panelInteractionDesigner.setEnabled(false);
            SpinnerTheta.setEnabled(false);
            SpinnerThetab.setEnabled(false);
            SpinnerPhi.setEnabled(false);
            SpinnerPhib.setEnabled(false);
            SpinnerOmega.setEnabled(false);
            lblTheta.setEnabled(false);
            lblThetab.setEnabled(false);
            lblPhi.setEnabled(false);
            lblPhib.setEnabled(false);
            lblOmega.setEnabled(false);
            lblSigmanm.setEnabled(false);
            labelSigmaX.setEnabled(false);
//			labelSigmaY.setEnabled(false);
//			labelSigmaZ.setEnabled(false);
            textFieldSigmaX.setEnabled(false);
        }else if(activate == 1){
            panelInteractionDesigner.setEnabled(true);
            SpinnerTheta.setEnabled(true);
            SpinnerThetab.setEnabled(true);
            SpinnerPhi.setEnabled(true);
            SpinnerPhib.setEnabled(true);
            SpinnerOmega.setEnabled(true);
            lblTheta.setEnabled(true);
            lblThetab.setEnabled(true);
            lblPhi.setEnabled(true);
            lblPhib.setEnabled(true);
            lblOmega.setEnabled(true);
            lblSigmanm.setEnabled(true);
            labelSigmaX.setEnabled(true);
            textFieldSigmaX.setEnabled(true);

        }
    }

    public void changerateunits(int spinnerReactantval, int spinnerProductval)
    {

        if(spinnerReactantval==0){
            if(chckbx2dReaction.isSelected()){//for a 2D Reaction
                lblForwardRate.setText("On Rate (nm-2/us)");//conditional
            }else{//for a 3D Reaction
                lblForwardRate.setText("On Rate (uM/s)");//conditional
            }
        }else if(spinnerReactantval==1){
            lblForwardRate.setText("On Rate (1/s)");//conditional
        }else if(spinnerReactantval==2){
            if(chckbx2dReaction.isSelected()){//for a 2D Reaction
                lblForwardRate.setText("On Rate (nm2/us)");//conditional
            }else{//for a 3D Reaction
                lblForwardRate.setText("On Rate (uM-1/s)");//conditional
            }
        }

        if(chckbxReversibleReaction.isSelected()){
            if(spinnerProductval==0){
                if(chckbx2dReaction.isSelected()){//for a 2D Reaction
                    lblBackwardRate.setText("Off Rate (nm-2/us)");//conditional
                }else{//for a 3D Reaction
                    lblBackwardRate.setText("Off Rate (uM/s)");//conditional
                }
            }else if(spinnerProductval==1){
                lblBackwardRate.setText("Off Rate (1/s)");//conditional
            }else if(spinnerProductval==2){
                if(chckbx2dReaction.isSelected()){//for a 2D Reaction
                    lblBackwardRate.setText("Off Rate (nm2/us)");//conditional
                }else{//for a 3D Reaction
                    lblBackwardRate.setText("Off Rate (uM-1/s)");//conditional
                }
            }
        }else{
            lblBackwardRate.setText("Off Rate");
        }
    }

    public void fillinNodes(ArrayList <Molecule> MolCla,DefaultTreeModel treemodel)
    {
        DefaultMutableTreeNode rootreactant1 = (DefaultMutableTreeNode) treemodel.getRoot();
        flushTree(treemodel);
        for(int i = 0; i<MolCla.size();i++)
        {//molecules
            DefaultMutableTreeNode nodemolecule = new DefaultMutableTreeNode(MolCla.get(i).getName());
            rootreactant1.add(nodemolecule);
            for(int j=0; j<MolCla.get(i).getnumInterface();j++)
            {//interfaces
                DefaultMutableTreeNode interfacereactant1 = new DefaultMutableTreeNode(MolCla.get(i).getinterfaceNames(j));
                nodemolecule.add(interfacereactant1);
                for(int k=0; k<MolCla.get(i).getnumInterfaceStates(j);k++)
                {//states
                    //System.out.println(MolCla.get(i).getinterfaceStateNames(j,k));
                    DefaultMutableTreeNode interfacestatereactant1 = new DefaultMutableTreeNode(MolCla.get(i).getinterfaceStateNames(j,k));
                    interfacereactant1.add(interfacestatereactant1);
                }
            }
            treemodel.reload(rootreactant1);
        }
    }

    public int isSameMolecule(ArrayList <Molecule> MolCla, int molecIndex) // check if there is a previous molecule added to the system
    {
        int same = 0, flagintstate = 0, flagint = 0;;

        for(int i=0;i<MolCla.size();i++)
        {
            if(i!=molecIndex)
            {//now check the other molecules
                if(MolCla.get(molecIndex).getmolOnMembrane()==MolCla.get(i).getmolOnMembrane()
                        && MolCla.get(molecIndex).getmolDiffR()==MolCla.get(i).getmolDiffR()
                        && MolCla.get(molecIndex).getmolDiffT()==MolCla.get(i).getmolDiffT()
                        && MolCla.get(molecIndex).getmolCount()==MolCla.get(i).getmolCount()
                        && MolCla.get(molecIndex).getnumInterface()==MolCla.get(i).getnumInterface())
                {//now check the interfaces
                    flagint = 0;
                    for(int j=0;j<MolCla.get(molecIndex).getnumInterface();j++)
                    {
                        if(MolCla.get(molecIndex).getinterfaceNames(j)==MolCla.get(i).getinterfaceNames(j)
                                && MolCla.get(molecIndex).getinterfaceCoords(j, 0)==MolCla.get(i).getinterfaceCoords(j, 0)
                                && MolCla.get(molecIndex).getinterfaceCoords(j, 1)==MolCla.get(i).getinterfaceCoords(j, 1)
                                && MolCla.get(molecIndex).getinterfaceCoords(j, 2)==MolCla.get(i).getinterfaceCoords(j, 2)
                                && MolCla.get(molecIndex).getnumInterfaceStates(j)==MolCla.get(i).getnumInterfaceStates(j))
                        {//now check the interface states
                            flagintstate = 0;
                            for(int k=0;k<MolCla.get(i).getnumInterfaceStates(j);k++)
                            {
                                if(MolCla.get(molecIndex).getinterfaceStateNames(j, k).equals(MolCla.get(i).getinterfaceStateNames(j, k)))
                                {
                                    flagintstate = flagintstate+1;//exactly same states
                                }
                            }
                            if(flagintstate==MolCla.get(i).getnumInterfaceStates(j)){
                                flagint = flagint + 1;//exactly same interfaces
                            }
                        }
                    }
                    if(flagint==MolCla.get(i).getnumInterface()){
                        same = 1;
                        break;
                    }

                }else if(MolCla.get(molecIndex).getName().equals(MolCla.get(i).getName())){ // do they have the same name??
                    same = 1;
                    break;
                }

            }
        }

        return same;
    }

    public int isSameReaction(ArrayList <Reaction> Reax, int reacIndex) // check if there is a previous reaction already added to the system
    {
        int same = 0, flagcheck = 0;
        //	Reax.get(reacIndex).get2D()
        //	Reax.get(reacIndex).getmolOnMembrane()
        //	Reax.get(reacIndex).getNumProduct()
        //	Reax.get(reacIndex).getNumReactant()
        //	Reax.get(reacIndex).getSigma()
        //	Reax.get(reacIndex).getAngle() double[]
        //	Reax.get(reacIndex).getProduct() int[]
        //	Reax.get(reacIndex).getReactant() int[]
        //	Reax.get(reacIndex).getReactantIndiv(j)
        //	Reax.get(reacIndex).getProductIndiv(j)
        //	Reax.get(reacIndex).getAngleIndiv(j)

        for(int i=0;i<Reax.size();i++)
        {
            //Reax.get(i).printReaction();

            if(i!=reacIndex)
            {//now check the other reactions
                if(Reax.get(reacIndex).getRever()==Reax.get(i).getRever()
                        && Reax.get(reacIndex).getmolOnMembrane()==Reax.get(i).getmolOnMembrane()
                        && Reax.get(reacIndex).getNumProduct()==Reax.get(i).getNumProduct()
                        && Reax.get(reacIndex).getNumReactant()==Reax.get(i).getNumReactant()
                        && Reax.get(reacIndex).getSigma()==Reax.get(i).getSigma())
                {//now check the angles and
                    flagcheck = 0;
                    for(int j=0;j<5;j++)
                    {
                        if(Reax.get(reacIndex).getReactantIndiv(j)==Reax.get(i).getReactantIndiv(j)
                                && Reax.get(reacIndex).getProductIndiv(j)==Reax.get(i).getProductIndiv(j)
                                && Reax.get(reacIndex).getAngleIndiv(j)==Reax.get(i).getAngleIndiv(j)
                                && Reax.get(reacIndex).getReactantIndiv(j+3)==Reax.get(i).getReactantIndiv(j+3)
                                && Reax.get(reacIndex).getProductIndiv(j+3)==Reax.get(i).getProductIndiv(j+3))
                        {
                            flagcheck = flagcheck + 1;//exactly same reactions
                        }
                    }
                    if(flagcheck==3){
                        same = 1;
                        break;
                    }
                }
            }
        }

        return same;
    }

    public String stringizeMolecule(int[] reactdata, ArrayList <Molecule> MolCla){

        String moleculeBNGL="";
        int mol1 = reactdata[0];
        int int1 = reactdata[1];
        int ist1 = reactdata[2]; //State changes may be added later

        moleculeBNGL = moleculeBNGL + MolCla.get(mol1).getName()+"(";


        System.out.println(int1);

        for(int j=0;j<MolCla.get(mol1).getnumInterface();j++) {//firstreactant
            System.out.println(j);
            if(j==int1) {//this is the binding interface
                moleculeBNGL = moleculeBNGL + MolCla.get(mol1).getinterfaceNames(j);
                if(ist1 > -1) {
                    moleculeBNGL = moleculeBNGL + "~" + MolCla.get(mol1).getinterfaceStateNames(int1,ist1);
                }
            }
        }
        moleculeBNGL = moleculeBNGL +")";
        return moleculeBNGL;
    }

    public String stringizeComplex(int[] reactdata1, int[] reactdata2, ArrayList <Molecule> MolCla)
    {

        String complexstring="", r1="", r2="";
        int mol1, int1, ist1, mol2, int2, ist2;

        mol1 = reactdata1[0];
        int1 = reactdata1[1];
        ist1 = reactdata1[2];

        mol2 = reactdata2[0];
        int2 = reactdata2[1];
        ist2 = reactdata2[2];

        r1 = MolCla.get(mol1).getName()+"(";
        r2 = MolCla.get(mol2).getName()+"(";

        for(int j=0;j<MolCla.get(mol1).getnumInterface();j++){//firstreactant
            if(j==int1){//this is the binding interface
                r1 = r1 + MolCla.get(mol1).getinterfaceNames(j);
                if(ist1>-1){
                    r1 = r1 + "~" + MolCla.get(mol1).getinterfaceStateNames(int1,ist1);
                }
                r1 = r1 + "!1";
            }
        }

        for(int j=0;j<MolCla.get(mol2).getnumInterface();j++){//secondreactant
            if(j==int2){//this is the binding interface
                r2 = r2 + MolCla.get(mol2).getinterfaceNames(j);
                if(ist2>-1){
                    r2 = r2 + "~" + MolCla.get(mol2).getinterfaceStateNames(int2,ist2);
                }
                r2 = r2 + "!1";
            }
        }

        complexstring = r1+")."+r2+")";

        return complexstring;
    }

    public String stringizeReaction(ArrayList <Molecule> MolCla)
    {

        String reactionstring = "", reactantstring="", productstring="";
        int numproducts = (Integer) spinnerNumProduct.getValue();
        int numreactants = (Integer) spinnerNumReactant.getValue();
        int tempP1[] = new int[3];
        int tempP2[] = new int[3];

        for(int i=0;i<3;i++){
            tempP1[i] = productdata[i];
            tempP2[i] = productdata[i+3];
        }

        if(numreactants==1)
        {
            reactantstring = reactantstring + stringizeMolecule(react1data, MolCla);
        }else if(numreactants==2)
        {
            reactantstring = reactantstring + stringizeMolecule(react1data, MolCla) + " + " + stringizeMolecule(react2data, MolCla);
        }

        if(numproducts==1){	//COMPLEX Formed
            if(numreactants==2){
                productstring = stringizeComplex(react1data, react2data, MolCla);
            }else{
                productstring = productstring + stringizeMolecule(tempP1, MolCla);
            }
        }else if(numproducts==2){
            productstring = productstring + stringizeMolecule(tempP1, MolCla) + " + " + stringizeMolecule(tempP2, MolCla);
        }

        if(chckbxReversibleReaction.isSelected()){
            reactionstring = reactantstring + " <-> " + productstring;
        }else{
            reactionstring = reactantstring + " -> " + productstring;
        }
        return reactionstring;
    }

    public void flushTree(DefaultTreeModel treemodel)
    {
        DefaultMutableTreeNode root = (DefaultMutableTreeNode) treemodel.getRoot();
        root.removeAllChildren();
        treemodel.reload(root);
    }

    public void clearReactionsTable()
    {
        DefaultTableModel model = (DefaultTableModel) tableReaction.getModel();
        int numreactions = model.getRowCount();

        for(int i=0;i<numreactions;i++){
            model.removeRow(i);
        }
    }

    public void fillinNodesReactantSpecific(ArrayList <Molecule> MolCla,DefaultTreeModel treemodel, int reactdata[])
    {//this method fills in product tree for the A+B->C reaction
        DefaultMutableTreeNode rootreactant1 = (DefaultMutableTreeNode) treemodel.getRoot();
        flushTree(treemodel);
        for(int i = 0; i<MolCla.size();i++)
        {//molecules
            if(i==reactdata[0])
            {
                DefaultMutableTreeNode nodemolecule = new DefaultMutableTreeNode(MolCla.get(i).getName());
                rootreactant1.add(nodemolecule);
                for(int j=0; j<MolCla.get(i).getnumInterface();j++)
                {//interfaces
                    DefaultMutableTreeNode interfacereactant1 = new DefaultMutableTreeNode(MolCla.get(i).getinterfaceNames(j));
                    nodemolecule.add(interfacereactant1);
                    for(int k=0; k<MolCla.get(i).getnumInterfaceStates(j);k++)
                    {//states
                        DefaultMutableTreeNode interfacestatereactant1 = new DefaultMutableTreeNode(MolCla.get(i).getinterfaceStateNames(j,k));
                        interfacereactant1.add(interfacestatereactant1);
                    }
                }
                treemodel.reload(rootreactant1);
            }
        }
    }

    public void fillinNodesFormComplex(ArrayList <Molecule> MolCla,DefaultTreeModel treemodel, int reactdata1[], int reactdata2[])
    {

        String complex="";

        DefaultMutableTreeNode rootproduct1 = (DefaultMutableTreeNode) treemodel.getRoot();
        flushTree(treemodel);

        complex = stringizeComplex(reactdata1,reactdata2,MolCla);
        DefaultMutableTreeNode nodemolecule = new DefaultMutableTreeNode(complex);
        rootproduct1.add(nodemolecule);
        treemodel.reload(rootproduct1);

        //select/highlight the path here automatically for a complex
        TreePath path = new TreePath(treemodel.getPathToRoot(nodemolecule));
        treeProduct1.setSelectionPath(path);
        treeProduct1.scrollPathToVisible(path);

    }

    public void generateProductforTrees(DefaultTreeModel product1treemodel, DefaultTreeModel product2treemodel)
    {
        int flaggoodselect=1;
        flushTree(product1treemodel);
        flushTree(product2treemodel);

        //		check for successful selection
        if((react1data[0]==-1 || react2data[0]==-1)
                || (react1data[1]==-1 || react2data[1]==-1))
        {
            flaggoodselect=0;
        }

        lblErrorReaction.setText("");

        try {
            spinnerNumReactant.commitEdit();
        } catch ( java.text.ParseException e ) {}
        int numReactant = (Integer) spinnerNumReactant.getValue();

        try {
            spinnerNumProduct.commitEdit();
        } catch ( java.text.ParseException e ) {}
        int numProduct = (Integer) spinnerNumProduct.getValue();

        switch(numProduct)
        {
            case 0: //No Products set (can be either decay reaction or an error)
                if(numReactant==2){//A+B->0
                    lblErrorReaction.setText("Illegal reaction, A+B->0.");
                }else if(numReactant==1){//decay reaction A->0

                }else if(numReactant==0){//0->0
                    lblErrorReaction.setText("Please select at least one molecule.");
                }
                break;
            case 1: //One product, could be a unimolecular reaction A->B or a combination reaction. Only Product1 tree will be changed
                if(numReactant==2){//A+B->C
                    if(flaggoodselect==1)
                    {
                        fillinNodesFormComplex(MolCla,product1treemodel,react1data,react2data);
                        activateInterfaceDesigner(1);
                        panelDrawInterface.repaint();
                    }else{//nothing got selected
                        lblErrorReaction.setText("Please select two molecules&interfaces.");
                    }
                }else if(numReactant==1){//A->B
                    fillinNodes(MolCla,product1treemodel);
                }else if(numReactant==0){//0->A birth
                    fillinNodes(MolCla,product1treemodel);
                }
                break;
            case 2: //A->B+C or A+B->C+D. Have to change both Product1 and Product2 trees
                if(numReactant==2){//A+B->C+D
                    fillinNodesReactantSpecific(MolCla,product1treemodel,react1data);
                    fillinNodesReactantSpecific(MolCla,product2treemodel,react2data);
                }else if(numReactant==1){//A->B+C

                }else if(numReactant==0){//0->A+B birth
                    lblErrorReaction.setText("Illegal reaction, 0->A+B.");
                }
                break;
        }

    }

    private void createEvents()
    {
        parminput.addActionListener(new ActionListener()
        {

            public void actionPerformed(ActionEvent e)
            {
                if(MolCla.size()>0){

                    DecimalFormat four = new DecimalFormat("#0.0000");
                    int numberofmolecules, numinterfaces, i, j, runninginterfaceindex=0;
                    double DiffT, DiffR, x, y, z;
                    boolean isMem;
                    String filename, interfacename, reactantname, reactantnamefix;

                    JFileChooser fc = new JFileChooser();
                    fc.setCurrentDirectory(new java.io.File(".")); // start at application current directory
                    fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
                    int returnVal = fc.showOpenDialog(getContentPane());
                    if(returnVal == JFileChooser.APPROVE_OPTION)
                    {
                        File yourFolder = fc.getSelectedFile();
                        try {
                            for(j=0;j<MolCla.size();j++){

                                filename = MolCla.get(j).getName().replaceAll("\\s+", "");
                                numberofmolecules = MolCla.get(j).getmolCount();
                                DiffT = MolCla.get(j).getmolDiffT();
                                DiffR = MolCla.get(j).getmolDiffR();
                                isMem = MolCla.get(j).getmolOnMembrane();
                                numinterfaces = MolCla.get(j).getnumInterface();

                                FileWriter stream = new FileWriter(yourFolder + "/" + filename + ".mol");
                                BufferedWriter out = new BufferedWriter(stream);
                                out.write("##");
                                out.newLine();
                                out.write("# " + filename + "molecule information file");
                                out.newLine();
                                out.write("##");
                                out.newLine();
                                out.newLine();
                                out.write("Name	= " +filename);
                                out.newLine();
                                out.write("copies = " + String.valueOf(numberofmolecules));
                                out.newLine();
                                /*out.write("isLipid = ");
                                if (chckbxIsLipid.isSelected())
                                {
                                    out.write("true");
                                }else{
                                    out.write("false");
                                }
                                out.newLine();*/
                                out.write("checkOverlap = true");
                                out.newLine();
                                out.write("mass = 1 #defualt value");
                                out.newLine();
                                out.newLine();

                                out.write("#translational diffusion");
                                out.newLine();
                                out.write("D = [" + four.format(DiffT) + "," + four.format(DiffT) + "," + four.format(0.0));

                                out.newLine();
                                out.newLine();

                                out.write("#rotational diffusion");
                                out.newLine();
                                out.write("Dr = [" + four.format(DiffR) + "," + four.format(DiffR) + "," + four.format(0.0) + "]");
                                out.newLine();
                                out.newLine();

//                                out.write("#interfaces"+ "\t" +"#unique interface id");
//                                out.newLine();
//                                out.write(String.valueOf(numinterfaces));
//                                for(i = 0; i<numinterfaces; i++)
//                                {
//                                    out.write("\t" + runninginterfaceindex);
//                                    runninginterfaceindex++;
//                                }
//                                out.newLine();
//                                out.newLine();
//
//                                out.write("#interface"+ "\t" +"name"+ "\t" +"#extra states"+ "\t" +"states");
//                                out.newLine();
//                                for(i = 0; i<numinterfaces; i++)
//                                {
//                                    interfacename = MolCla.get(j).getinterfaceNames(i).replaceAll("\\s+", "");
//                                    out.write(String.valueOf(i) + "\t" +interfacename + "\t" + String.valueOf(MolCla.get(j).getnumInterfaceStates(i))+ "\t" + MolCla.get(j).getinterfaceStateNamesConcat(i));
//                                    out.newLine();
//                                }
                                out.newLine();

                                out.write("#Coords"+ "\t" +"x"+ "\t" +	"y"+ "\t" +	"z");
                                out.newLine();
                                out.write("COM	0.0000	0.0000	0.0000");
                                out.newLine();
                                for(i = 0; i<numinterfaces; i++)
                                {
                                    interfacename = MolCla.get(j).getinterfaceNames(i).replaceAll("\\s+", "");
                                    x = MolCla.get(j).getinterfaceCoords(i, 0);
                                    y = MolCla.get(j).getinterfaceCoords(i, 1);
                                    z = MolCla.get(j).getinterfaceCoords(i, 2);

                                    out.write(interfacename + "\t" + four.format(x) + "\t" + four.format(y) + "\t" + four.format(z));
                                    out.newLine();
                                }
                                out.newLine();

                                out.close();
                            }
                            JOptionPane.showMessageDialog(null, "FPR molecule .mol files were successfully exported.");

                            FileWriter stream2 = new FileWriter(yourFolder + "/model.inp");
                            BufferedWriter out = new BufferedWriter(stream2);
                            if(NumberUtils.isCreatable(textFieldTimeStep.getText()))
                            {
                                if(NumberUtils.isCreatable(textFieldNUMiterations.getText()))
                                {
                                    if(NumberUtils.isCreatable(textFieldStatOutFreq.getText()))
                                    {
                                        if(NumberUtils.isCreatable(textFieldConfOutFreq.getText()))
                                        {
                                            if(NumberUtils.isCreatable(textFieldBoxX.getText())&&NumberUtils.isCreatable(textFieldBoxY.getText())&&NumberUtils.isCreatable(textFieldBoxZ.getText()))
                                            {

                                                // save button code here - Nomo
                                                out.write("# Input file");
                                                out.newLine();
                                                out.newLine();
                                                out.write("start parameters");
                                                out.newLine();
                                                out.write("	" + "nItr = " + String.valueOf(textFieldNUMiterations.getText()));
                                                out.newLine();
                                                out.write("	" + "timeStep = " + String.valueOf(textFieldTimeStep.getText()));
                                                out.newLine();
                                                out.newLine();
                                                out.write("	" + "timeWrite = " + textFieldConfOutFreq.getText());
                                                out.newLine();
                                                out.write("	" + "trajWrite = " + textFieldStatOutFreq.getText() + "\n");
                                                out.newLine();
                                                out.write("	" + "restartWrite = " + textFieldRestartOutFreq.getText());
                                                out.newLine();
                                                out.write("	" + "overlapSepLimit = .6");
                                                out.newLine();
                                                out.write("end parameters\n");
                                                out.newLine();
                                                out.newLine();
                                                out.write("start boundaries");
                                                out.newLine();
                                                out.write("	" + "WaterBox = [");
                                                out.write(textFieldBoxX.getText() + "," + textFieldBoxY.getText() + "," + textFieldBoxZ.getText() + "]");
                                                out.newLine();
                                                out.write("	" + "implicitLipid = ");
                                                if (chckbxIsImplicitLipid.isSelected())
                                                {
                                                    out.write("true");
                                                }else{
                                                    out.write("false");
                                                }
                                                out.newLine();
                                                out.write("	" + "xBCtype = ");
                                                if (reflectiveX.isSelected()) {
                                                    out.write("reflect");
                                                }else{
                                                    out.write("periodic");
                                                }
                                                out.newLine();
                                                out.write("	" + "yBCtype = ");
                                                if (reflectiveY.isSelected()) {
                                                    out.write("reflect");
                                                }else{
                                                    out.write("periodic");
                                                }
                                                out.newLine();
                                                out.write("	" + "zBCtype = ");
                                                if (reflectiveZ.isSelected()) {
                                                    out.write("reflect");
                                                }else{
                                                    out.write("periodic");
                                                }
                                                out.newLine();
                                                out.write("end boundaries");
                                                out.newLine();
                                                out.newLine();
                                                out.write("start molecules");
                                                out.newLine();
                                                for(j=0;j<MolCla.size();j++)
                                                {
                                                    filename = MolCla.get(j).getName().replaceAll("\\s+", "");
                                                    out.write("	");
                                                    out.write(filename);
                                                    out.newLine();
                                                }
                                                out.write("end molecules");
                                                out.newLine();
                                                out.newLine();
                                                out.write("start reactions");

                                                DefaultTableModel model = (DefaultTableModel) tableReaction.getModel();

                                                if(Reax.size()>0)
                                                {
                                                    out.newLine();
                                                    for(j=0;j<Reax.size();j++)
                                                    {

                                                        out.write("	");
                                                        out.write(model.getValueAt(j, 0).toString());
                                                        System.out.println(model.getValueAt(j, 0).toString());
                                                        out.newLine();
                                                        out.write("	onRate = " + model.getValueAt(j, 2));
                                                        out.newLine();
                                                        out.write("	offRate = " + model.getValueAt(j, 4));
                                                        out.newLine();

                                                        out.write("	sigma = " + model.getValueAt(j, 1));
                                                        out.newLine();
                                                        String out_norm1 = Reax.get(j).norm1.toString();
                                                        String out_norm2 = Reax.get(j).norm2.toString();
                                                        out.write("	norm1 = [" + out_norm1.substring(1, out_norm1.length()-1) + "]");
                                                        out.newLine();
                                                        out.write("	norm2 = [" + out_norm2.substring(1, out_norm2.length()-1) + "]");
                                                        out.newLine();
                                                        out.write("\t bindRadSameCom = 1.5 #scales sigma to define distance");
                                                        out.newLine();
                                                        out.write("\tloopCoopFactor = 0.001");
                                                        out.newLine();
                                                        out.write("\tlength3Dto2D = " + 2 * Float.parseFloat((String)model.getValueAt(j, 1)) + " #default 2*sigma");
                                                        out.newLine();
                                                        out.write("	assocAngles = [");
                                                        String theta1 = Double.toString(Reax.get(j).angles[1]);
                                                        String theta2 = Double.toString(Reax.get(j).angles[2]);
                                                        String phi1s;
                                                        String phi2s;
                                                        float phi1 = (float) Reax.get(j).angles[3];
                                                        if (phi1 == -999) {
                                                            phi1s = "-";
                                                        } else {
                                                            phi1s = Double.toString(phi1);
                                                        }
                                                        float phi2 = (float) Reax.get(j).angles[4];
                                                        if (phi2 == -999) {
                                                            phi2s = "-";
                                                        } else {
                                                            phi2s = Double.toString(phi2);
                                                        }
                                                        float om = (float) Reax.get(j).angles[1];

                                                        out.write(theta1 + "," + theta2 + "," + phi1s +"," + phi2s + "," + om + "]");
                                                        out.newLine();
                                                        out.newLine();
                                                    }
                                                    out.write("end reactions");
                                                }else{
                                                    JOptionPane.showMessageDialog(null, "FPR reactions file were not successfully exported. No reactions present.");
                                                }
                                                out.close();

                                            }else{
                                                JOptionPane.showMessageDialog(null, "FPR parm file was not successfully exported. Please check Parameters>Box size");
                                            }
                                        }else{
                                            JOptionPane.showMessageDialog(null, "FPR parm file was not successfully exported. Please check Parameters>Configuration output frequency.");
                                        }
                                    }else{
                                        JOptionPane.showMessageDialog(null, "FPR parm file was not successfully exported. Please check Parameters>Statistics output frequency.");
                                    }
                                }else{
                                    JOptionPane.showMessageDialog(null, "FPR parm file was not successfully exported. Please check Parameters>Number of iterations.");
                                }
                            }else{
                                JOptionPane.showMessageDialog(null, "FPR parm file was not successfully exported. Please check Parameters>Time step.");
                            }
                            out.close();


                        } catch (IOException e1) {
                            // TODO Auto-generated catch block
                            e1.printStackTrace();
                        }
                    }

                }else{
                    JOptionPane.showMessageDialog(null, "FPR input files were not successfully exported. Please add a molecule.");
                }
            }
        });

        mntmExit.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent arg0)
            {
                int ret = JOptionPane.showConfirmDialog(null, "Are you sure?");
                if (ret==JOptionPane.YES_OPTION)
                    System.exit(0);
            }
        });

        btnAddInterface.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent arg0)
            {
                lblMoleculeNameEmptyError.setText("");
                DefaultTableModel model = (DefaultTableModel) tableInterface.getModel();
                if(!textFieldInterfaceName.getText().trim().equals(""))
                {
                    if(NumberUtils.isCreatable(textFieldInterfaceX.getText())&&NumberUtils.isCreatable(textFieldInterfaceY.getText())&&NumberUtils.isCreatable(textFieldInterfaceZ.getText()))
                    {
                        JSpinner NumberofStates = new JSpinner(new SpinnerNumberModel(2, 2, 10, 1));
                        JTextField NamesofStates = new JTextField();
                        final JComponent[] inputs = new JComponent[] {
                                new JLabel("Number of interface states (Max. 10)"),
                                NumberofStates,
                                new JLabel("List of interface states (e.g. U, P)"),
                                NamesofStates
                        };

                        int dialogButton = JOptionPane.YES_NO_OPTION;
                        int dialogResult = JOptionPane.showConfirmDialog(null, "Can this site exist in additional states?", "Additional states?", dialogButton);
                        if(dialogResult == 0)
                        {
                            int result = JOptionPane.showConfirmDialog(null, inputs, "Can this site exist in additional states?", JOptionPane.PLAIN_MESSAGE);
                            if (result == JOptionPane.OK_OPTION) {

                                int StateDuplicate = 0; //Check for possible interface state duplication
                                String[] strNamesofStates = NamesofStates.getText().split(","); //System.out.println(strNamesofStates.length);
                                for(int i=0;i<strNamesofStates.length;i++){
                                    strNamesofStates[i]=strNamesofStates[i].trim();
                                    if(strNamesofStates.length>1){
                                        for(int j=i+1;j<strNamesofStates.length;j++){
                                            if(strNamesofStates[i].equals(strNamesofStates[j])){
                                                StateDuplicate = 1;
                                            }
                                        }
                                    }
                                }

                                int InterfaceDuplicate = 0;//Check for possible interface duplication
                                for(int i=0;i<model.getRowCount();i++){
                                    if(Objects.equals(model.getValueAt(i, 0).toString().trim(), textFieldInterfaceName.getText().toString().trim())){//Check only for the interfacename duplication

//									if(Objects.equals(model.getValueAt(i, 0).toString().trim(), textFieldInterfaceName.getText().toString().trim()) ||
//											(Objects.equals(Double.parseDouble(model.getValueAt(i, 1).toString().trim()),  Double.parseDouble(textFieldInterfaceX.getText().toString().trim())) &&
//													Objects.equals(Double.parseDouble(model.getValueAt(i, 2).toString().trim()),  Double.parseDouble(textFieldInterfaceY.getText().toString().trim())) &&
//													Objects.equals(Double.parseDouble(model.getValueAt(i, 3).toString().trim()),  Double.parseDouble(textFieldInterfaceZ.getText().toString().trim())) )){
                                        InterfaceDuplicate = 1;
                                    }
                                }

                                if(InterfaceDuplicate == 0){ //duplicate interface state entries
                                    if(StateDuplicate == 0){ //duplicate interface state entries
                                        if(!strNamesofStates[0].equals("") && strNamesofStates.length==(int)NumberofStates.getValue() && (int)NumberofStates.getValue()>0){
                                            model.addRow(new Object[]{textFieldInterfaceName.getText(),Double.parseDouble(textFieldInterfaceX.getText()),Double.parseDouble(textFieldInterfaceY.getText()),Double.parseDouble(textFieldInterfaceZ.getText()),NumberofStates.getValue(),NamesofStates.getText()});
                                            panelMolecule.repaint();
                                        }else if(strNamesofStates[0].equals("") && (strNamesofStates.length-1)==(int)NumberofStates.getValue() && (int)NumberofStates.getValue()==0){
                                            model.addRow(new Object[]{textFieldInterfaceName.getText(),Double.parseDouble(textFieldInterfaceX.getText()),Double.parseDouble(textFieldInterfaceY.getText()),Double.parseDouble(textFieldInterfaceZ.getText()),NumberofStates.getValue(),"-"});
                                            panelMolecule.repaint();
                                        }else{
                                            JOptionPane.showMessageDialog(null, "Mismatch! Number of states is different than the states list provided.", "Error", JOptionPane.ERROR_MESSAGE);
                                        }
                                    }else{
                                        JOptionPane.showMessageDialog(null, "Duplicate state entries.", "Error", JOptionPane.ERROR_MESSAGE);
                                    }
                                }else{
                                    JOptionPane.showMessageDialog(null, "Duplicate interface entries.", "Error", JOptionPane.ERROR_MESSAGE);
                                }
                            }
                        } else
                        {
                            int InterfaceDuplicate = 0;//Check for possible interface duplication
                            for(int i=0;i<model.getRowCount();i++){
                                if(Objects.equals(model.getValueAt(i, 0).toString().trim(), textFieldInterfaceName.getText().toString().trim())){//Check only for the interfacename duplication

//								if(Objects.equals(model.getValueAt(i, 0).toString().trim(), textFieldInterfaceName.getText().toString().trim()) ||
//										(Objects.equals(Double.parseDouble(model.getValueAt(i, 1).toString().trim()),  Double.parseDouble(textFieldInterfaceX.getText().toString().trim())) &&
//												Objects.equals(Double.parseDouble(model.getValueAt(i, 2).toString().trim()),  Double.parseDouble(textFieldInterfaceY.getText().toString().trim())) &&
//												Objects.equals(Double.parseDouble(model.getValueAt(i, 3).toString().trim()),  Double.parseDouble(textFieldInterfaceZ.getText().toString().trim())) )){
                                    InterfaceDuplicate = 1;
                                }
                            }

                            if(InterfaceDuplicate == 0){ //duplicate interface state entries
                                model.addRow(new Object[]{textFieldInterfaceName.getText(),Double.parseDouble(textFieldInterfaceX.getText()),Double.parseDouble(textFieldInterfaceY.getText()),Double.parseDouble(textFieldInterfaceZ.getText()),0,"-"});
                                panelMolecule.repaint();
                            }else{
                                JOptionPane.showMessageDialog(null, "Duplicate interface entries.", "Error", JOptionPane.ERROR_MESSAGE);
                            }
                        }

                    } else{
                        lblMoleculeNameEmptyError.setText("Check the interface coordinates.");
                    }

                } else{
                    lblMoleculeNameEmptyError.setText("Interface name should not be left blank.");
                }
            }
        });

        btnCalcDiff.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {

                DecimalFormat REAL_FORMATTER = new DecimalFormat("0.####");

                textFieldTransDifCoeff.setText(REAL_FORMATTER.format((calcTransDif())));
                textFieldRotDifCoeff.setText(REAL_FORMATTER.format((calcRotDif())));

            }
        });

        btnUpdateInterface.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                lblMoleculeNameEmptyError.setText("");
                DefaultTableModel model = (DefaultTableModel) tableInterface.getModel();

                if(tableInterface.getSelectedRow()!=-1)
                {
                    if(model.getRowCount()>1){

                        if(!textFieldInterfaceName.getText().trim().equals(""))
                        {
                            if(NumberUtils.isCreatable(textFieldInterfaceX.getText())&&NumberUtils.isCreatable(textFieldInterfaceY.getText())&&NumberUtils.isCreatable(textFieldInterfaceX.getText()))
                            {
                                JSpinner NumberofStates = new JSpinner(new SpinnerNumberModel(Integer.parseInt(tableInterface.getValueAt(tableInterface.getSelectedRow(), 4).toString()), 0, 10, 1));
                                JTextField NamesofStates = new JTextField();
                                if(Integer.parseInt(tableInterface.getValueAt(tableInterface.getSelectedRow(), 4).toString())==0){
                                    NamesofStates.setText("");
                                }else{
                                    NamesofStates.setText(tableInterface.getValueAt(tableInterface.getSelectedRow(), 5).toString());
                                }
                                final JComponent[] inputs = new JComponent[] {
                                        new JLabel("Number of interface states (Max. 10)"),
                                        NumberofStates,
                                        new JLabel("List of interface states (e.g. S1,S2)"),
                                        NamesofStates
                                };
                                int result = JOptionPane.showConfirmDialog(null, inputs, "Any States of the Interface?", JOptionPane.PLAIN_MESSAGE);
                                if (result == JOptionPane.OK_OPTION) {

                                    int StateDuplicate = 0; //Check for possible interface state duplication
                                    String[] strNamesofStates = NamesofStates.getText().split(","); //System.out.println(strNamesofStates.length);
                                    for(int i=0;i<strNamesofStates.length;i++){
                                        strNamesofStates[i]=strNamesofStates[i].trim();
                                        if(strNamesofStates.length>1){
                                            for(int j=i+1;j<strNamesofStates.length;j++){
                                                if(strNamesofStates[i].equals(strNamesofStates[j])){
                                                    StateDuplicate = 1;
                                                }
                                            }
                                        }
                                    }

                                    int InterfaceDuplicate = 0;//Check for possible interface duplication
                                    for(int i=0;i<model.getRowCount();i++){

                                        if(Objects.equals(model.getValueAt(i, 0).toString().trim(), textFieldInterfaceName.getText().toString().trim())){//Check only for the interfacename duplication

//										if(Objects.equals(model.getValueAt(i, 0).toString().trim(), textFieldInterfaceName.getText().toString().trim()) ||
//												(Objects.equals(Double.parseDouble(model.getValueAt(i, 1).toString().trim()),  Double.parseDouble(textFieldInterfaceX.getText().toString().trim())) &&
//														Objects.equals(Double.parseDouble(model.getValueAt(i, 2).toString().trim()),  Double.parseDouble(textFieldInterfaceY.getText().toString().trim())) &&
//														Objects.equals(Double.parseDouble(model.getValueAt(i, 3).toString().trim()),  Double.parseDouble(textFieldInterfaceZ.getText().toString().trim())) )){
                                            if(tableInterface.getSelectedRow()!=i){//make sure you're not comparing with your selected interface
                                                InterfaceDuplicate = 1;
                                            }
                                        }
                                    }

                                    if(InterfaceDuplicate == 0){ //duplicate interface state entries
                                        if(StateDuplicate == 0){ //duplicate interface state entries
                                            if(!strNamesofStates[0].equals("") && strNamesofStates.length==(int)NumberofStates.getValue() && (int)NumberofStates.getValue()>0){
                                                model.setValueAt(textFieldInterfaceName.getText(), tableInterface.getSelectedRow(), 0);
                                                model.setValueAt(Double.parseDouble(textFieldInterfaceX.getText()), tableInterface.getSelectedRow(), 1);
                                                model.setValueAt(Double.parseDouble(textFieldInterfaceY.getText()), tableInterface.getSelectedRow(), 2);
                                                model.setValueAt(Double.parseDouble(textFieldInterfaceZ.getText()), tableInterface.getSelectedRow(), 3);
                                                model.setValueAt(NumberofStates.getValue(), tableInterface.getSelectedRow(), 4);
                                                model.setValueAt(NamesofStates.getText(), tableInterface.getSelectedRow(), 5);
                                                panelMolecule.repaint();
                                            }else if(strNamesofStates[0].equals("") && (strNamesofStates.length-1)==(int)NumberofStates.getValue() && (int)NumberofStates.getValue()==0){
                                                model.setValueAt(textFieldInterfaceName.getText(), tableInterface.getSelectedRow(), 0);
                                                model.setValueAt(Double.parseDouble(textFieldInterfaceX.getText()), tableInterface.getSelectedRow(), 1);
                                                model.setValueAt(Double.parseDouble(textFieldInterfaceY.getText()), tableInterface.getSelectedRow(), 2);
                                                model.setValueAt(Double.parseDouble(textFieldInterfaceZ.getText()), tableInterface.getSelectedRow(), 3);
                                                model.setValueAt(NumberofStates.getValue(), tableInterface.getSelectedRow(), 4);
                                                model.setValueAt("-", tableInterface.getSelectedRow(), 5);
                                                panelMolecule.repaint();
                                            }else{
                                                JOptionPane.showMessageDialog(null, "Mismatch! Number of states is different than the states list provided.", "Error", JOptionPane.ERROR_MESSAGE);
                                            }
                                        }else{
                                            JOptionPane.showMessageDialog(null, "Duplicate state entries.", "Error", JOptionPane.ERROR_MESSAGE);
                                        }
                                    }else{
                                        JOptionPane.showMessageDialog(null, "Duplicate interface entries.", "Error", JOptionPane.ERROR_MESSAGE);
                                    }

                                }
                            }
                        }else{
                            lblMoleculeNameEmptyError.setText("Interface update error. Check inputs.");
                        }
                    } else{
                        lblMoleculeNameEmptyError.setText("Add an interface first."); //Can't change center of mass.
                    }

                }else{
                    lblMoleculeNameEmptyError.setText("Select an interface.");
                }
            }
        });

        tableInterface.addMouseListener(new MouseAdapter()
        {
            @Override
            public void mouseClicked(MouseEvent arg0)
            {
                DefaultTableModel model = (DefaultTableModel) tableInterface.getModel();
                if (tableInterface.getSelectedRow()>-1){
                    textFieldInterfaceName.setText(model.getValueAt(tableInterface.getSelectedRow(), 0).toString());
                    textFieldInterfaceX.setText(model.getValueAt(tableInterface.getSelectedRow(), 1).toString());
                    textFieldInterfaceY.setText(model.getValueAt(tableInterface.getSelectedRow(), 2).toString());
                    textFieldInterfaceZ.setText(model.getValueAt(tableInterface.getSelectedRow(), 3).toString());
                }
            }
        });

        btnDeleteInterface.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                lblMoleculeNameEmptyError.setText("");
                DefaultTableModel model = (DefaultTableModel) tableInterface.getModel();

                if(tableInterface.getSelectedRow()!=-1)
                {
                    if(tableInterface.getSelectedRow()>0)
                    {
                        if(model.getRowCount()>1){
                            model.removeRow(tableInterface.getSelectedRow());
                            textFieldInterfaceName.setText(null);
                            textFieldInterfaceX.setText(null);
                            textFieldInterfaceY.setText(null);
                            textFieldInterfaceZ.setText(null);
                            panelMolecule.repaint();
                        }else{
                            lblMoleculeNameEmptyError.setText("Add an interface first."); //Can't change center of mass.
                        }
                    }else{
                        lblMoleculeNameEmptyError.setText("Can't change center of mass.");
                    }
                }else{
                    lblMoleculeNameEmptyError.setText("Select an interface.");
                }
            }
        });

        /**
         * listener for adding new molecule
         */
        btnAddMolecule.addActionListener(new ActionListener()
        {

            /**
             * Action to preform on click
             * @param arg0
             */
            public void actionPerformed(ActionEvent arg0)
            {
                lblMoleculeNameEmptyError.setText("");
                DefaultTableModel model = (DefaultTableModel) tableMolecule.getModel();
                DefaultTableModel modelinterface = (DefaultTableModel) tableInterface.getModel();
                DefaultTreeModel modelreactant1 = (DefaultTreeModel) treeReactant1.getModel();
                DefaultTreeModel modelreactant2 = (DefaultTreeModel) treeReactant2.getModel();

                int isSphere = checkSphere();

                int numinterfaces = tableInterface.getRowCount()-1;

                if(!textFieldMoleculeName.getText().trim().equals(""))
                {
                    if(NumberUtils.isCreatable(textFieldMoleculeCount.getText()))
                    {
                        if(NumberUtils.isCreatable(textFieldTransDifCoeff.getText()))
                        {
                            if(NumberUtils.isCreatable(textFieldRotDifCoeff.getText()))
                            {
                                if(chckbxAnchoredToMembrane.isSelected())
                                {
                                    if(isSphere==1){//Force Drot=0;
                                        MolCla.add(new Molecule(textFieldMoleculeName.getText(), Integer.parseInt(textFieldMoleculeCount.getText()), Double.parseDouble(textFieldTransDifCoeff.getText()),0.0, true, numinterfaces, tableInterface));
                                    }else{
                                        MolCla.add(new Molecule(textFieldMoleculeName.getText(), Integer.parseInt(textFieldMoleculeCount.getText()), Double.parseDouble(textFieldTransDifCoeff.getText()),Double.parseDouble(textFieldRotDifCoeff.getText()), true, numinterfaces, tableInterface));
                                    }
                                    if(isSameMolecule(MolCla, MolCla.size()-1)==0)
                                    {
                                        if(isSphere==1){//Force Drot=0;
                                            model.addRow(new Object[]{textFieldMoleculeName.getText(),textFieldMoleculeCount.getText(),"Yes",textFieldTransDifCoeff.getText(),0.0,numinterfaces});
                                        }else{
                                            model.addRow(new Object[]{textFieldMoleculeName.getText(),textFieldMoleculeCount.getText(),"Yes",textFieldTransDifCoeff.getText(),textFieldRotDifCoeff.getText(),numinterfaces});
                                        }
                                        //MolCla.get(MolCla.size()-1).printMolecule();
                                        textFieldMoleculeName.setText(null);
                                        textFieldMoleculeCount.setText(null);
                                        textFieldTransDifCoeff.setText(null);
                                        textFieldRotDifCoeff.setText(null);
                                        chckbxAnchoredToMembrane.setSelected(false);
                                        textFieldInterfaceName.setText(null);
                                        textFieldInterfaceX.setText(null);
                                        textFieldInterfaceY.setText(null);
                                        textFieldInterfaceZ.setText(null);
                                        while(tableInterface.getRowCount()>1){
                                            modelinterface.removeRow(tableInterface.getRowCount()-1);
                                        }
                                        tableInterface.clearSelection();
                                        tableMolecule.clearSelection();
                                        panelMolecule.repaint();

                                        fillinNodes(MolCla,modelreactant1);
                                        fillinNodes(MolCla,modelreactant2);
                                        spinnerNumReactant.setValue(0);
                                        spinnerNumProduct.setValue(0);
                                        treeReactant1.clearSelection();
                                        treeReactant2.clearSelection();
                                        treeProduct1.clearSelection();
                                        treeProduct2.clearSelection();

                                    }else//same molecule present
                                    {
                                        lblMoleculeNameEmptyError.setText("Can't add. Same molecule already present.");
                                        MolCla.remove(MolCla.size()-1);
                                    }

                                } else{

                                    if(isSphere==1){//Force Drot=0;
                                        MolCla.add(new Molecule(textFieldMoleculeName.getText(), Integer.parseInt(textFieldMoleculeCount.getText()), Double.parseDouble(textFieldTransDifCoeff.getText()),0.0, false, numinterfaces, tableInterface));
                                    }else{
                                        MolCla.add(new Molecule(textFieldMoleculeName.getText(), Integer.parseInt(textFieldMoleculeCount.getText()), Double.parseDouble(textFieldTransDifCoeff.getText()),Double.parseDouble(textFieldRotDifCoeff.getText()), false, numinterfaces, tableInterface));
                                    }
                                    if(isSameMolecule(MolCla, MolCla.size()-1)==0)
                                    {
                                        if(isSphere==1){//Force Drot=0;
                                            model.addRow(new Object[]{textFieldMoleculeName.getText(),textFieldMoleculeCount.getText(),"No",textFieldTransDifCoeff.getText(),0.0,numinterfaces});
                                        }else{
                                            model.addRow(new Object[]{textFieldMoleculeName.getText(),textFieldMoleculeCount.getText(),"No",textFieldTransDifCoeff.getText(),textFieldRotDifCoeff.getText(),numinterfaces});
                                        }
                                        //MolCla.get(MolCla.size()-1).printMolecule();
                                        textFieldMoleculeName.setText(null);
                                        textFieldMoleculeCount.setText(null);
                                        textFieldTransDifCoeff.setText(null);
                                        textFieldRotDifCoeff.setText(null);
                                        chckbxAnchoredToMembrane.setSelected(false);
                                        textFieldInterfaceName.setText(null);
                                        textFieldInterfaceX.setText(null);
                                        textFieldInterfaceY.setText(null);
                                        textFieldInterfaceZ.setText(null);
                                        while(tableInterface.getRowCount()>1){
                                            modelinterface.removeRow(tableInterface.getRowCount()-1);
                                        }
                                        tableInterface.clearSelection();
                                        tableMolecule.clearSelection();
                                        panelMolecule.repaint();

                                        fillinNodes(MolCla,modelreactant1);
                                        fillinNodes(MolCla,modelreactant2);
                                        spinnerNumReactant.setValue(0);
                                        spinnerNumProduct.setValue(0);
                                        treeReactant1.clearSelection();
                                        treeReactant2.clearSelection();
                                        treeProduct1.clearSelection();
                                        treeProduct2.clearSelection();

                                    }else//same molecule present
                                    {
                                        lblMoleculeNameEmptyError.setText("Can't add. Same molecule already present.");
                                        MolCla.remove(MolCla.size()-1);
                                    }
                                }

                            } else{
                                lblMoleculeNameEmptyError.setText("Can't add. Rotational diffusion missing.");
                            }
                        } else{
                            lblMoleculeNameEmptyError.setText("Can't add. Translational diffusion missing.");
                        }
                    } else{
                        lblMoleculeNameEmptyError.setText("Can't add. Molecule count missing.");
                    }

                } else{
                    lblMoleculeNameEmptyError.setText("Can't add. Molecule name should not be left blank.");
                }
            }
        });

        btnUpdateMolecule.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                if(tableMolecule.getSelectedRow()!=-1)
                {

                    lblMoleculeNameEmptyError.setText("");
                    DefaultTableModel model = (DefaultTableModel) tableMolecule.getModel();
                    DefaultTableModel modelinterface = (DefaultTableModel) tableInterface.getModel();
                    DefaultTreeModel modelreactant1 = (DefaultTreeModel) treeReactant1.getModel();
                    DefaultTreeModel modelreactant2 = (DefaultTreeModel) treeReactant2.getModel();

                    int isSphere = checkSphere();
                    int numinterfaces = tableInterface.getRowCount()-1;

                    if(!textFieldMoleculeName.getText().trim().equals(""))
                    {
                        if(NumberUtils.isCreatable(textFieldMoleculeCount.getText()))
                        {
                            if(NumberUtils.isCreatable(textFieldTransDifCoeff.getText()))
                            {
                                if(NumberUtils.isCreatable(textFieldRotDifCoeff.getText()))
                                {
                                    if(chckbxAnchoredToMembrane.isSelected())
                                    {
                                        if(isSphere==1){//Force Drot=0;
                                            MolCla.set(tableMolecule.getSelectedRow(), new Molecule(textFieldMoleculeName.getText(), Integer.parseInt(textFieldMoleculeCount.getText()), Double.parseDouble(textFieldTransDifCoeff.getText()),0.0, true, numinterfaces, tableInterface));
                                        }else{
                                            MolCla.set(tableMolecule.getSelectedRow(), new Molecule(textFieldMoleculeName.getText(), Integer.parseInt(textFieldMoleculeCount.getText()), Double.parseDouble(textFieldTransDifCoeff.getText()),Double.parseDouble(textFieldRotDifCoeff.getText()), true, numinterfaces, tableInterface));
                                        }
                                    } else{
                                        if(isSphere==1){//Force Drot=0;
                                            MolCla.set(tableMolecule.getSelectedRow(), new Molecule(textFieldMoleculeName.getText(), Integer.parseInt(textFieldMoleculeCount.getText()), Double.parseDouble(textFieldTransDifCoeff.getText()),0.0, false, numinterfaces, tableInterface));
                                        }else{
                                            MolCla.set(tableMolecule.getSelectedRow(), new Molecule(textFieldMoleculeName.getText(), Integer.parseInt(textFieldMoleculeCount.getText()), Double.parseDouble(textFieldTransDifCoeff.getText()),Double.parseDouble(textFieldRotDifCoeff.getText()), false, numinterfaces, tableInterface));
                                        }
                                    }

                                    if(isSameMolecule(MolCla, tableMolecule.getSelectedRow())==0)
                                    {
                                        if(chckbxAnchoredToMembrane.isSelected())
                                        {
                                            model.setValueAt("Yes", tableMolecule.getSelectedRow(), 2);
                                        } else{
                                            model.setValueAt("No", tableMolecule.getSelectedRow(), 2);
                                        }
                                        model.setValueAt(textFieldMoleculeName.getText(), tableMolecule.getSelectedRow(), 0);
                                        model.setValueAt(textFieldMoleculeCount.getText(), tableMolecule.getSelectedRow(), 1);
                                        model.setValueAt(textFieldTransDifCoeff.getText(), tableMolecule.getSelectedRow(), 3);
                                        if(isSphere==1){//Force Drot=0;
                                            model.setValueAt("0.0", tableMolecule.getSelectedRow(), 4);
                                        }else{
                                            model.setValueAt(textFieldRotDifCoeff.getText(), tableMolecule.getSelectedRow(), 4);
                                        }
                                        model.setValueAt(numinterfaces, tableMolecule.getSelectedRow(), 5);
                                        //MolCla.get(MolCla.size()-1).printMolecule();
                                        textFieldMoleculeName.setText(null);//clear the form...
                                        textFieldMoleculeCount.setText(null);
                                        textFieldTransDifCoeff.setText(null);
                                        textFieldRotDifCoeff.setText(null);
                                        chckbxAnchoredToMembrane.setSelected(false);
                                        textFieldInterfaceName.setText(null);
                                        textFieldInterfaceX.setText(null);
                                        textFieldInterfaceY.setText(null);
                                        textFieldInterfaceZ.setText(null);
                                        while(tableInterface.getRowCount()>1){
                                            modelinterface.removeRow(tableInterface.getRowCount()-1);
                                        }

                                    }else//same molecule present
                                    {
                                        lblMoleculeNameEmptyError.setText("Can't update. Selected molecule already present and deleted.");
                                        MolCla.remove(tableMolecule.getSelectedRow());
                                        model.removeRow(tableMolecule.getSelectedRow());
                                    }

                                    tableInterface.clearSelection();
                                    tableMolecule.clearSelection();
                                    panelMolecule.repaint();

                                    fillinNodes(MolCla,modelreactant1);
                                    fillinNodes(MolCla,modelreactant2);
                                    spinnerNumReactant.setValue(0);
                                    spinnerNumProduct.setValue(0);
                                    treeReactant1.clearSelection();
                                    treeReactant2.clearSelection();
                                    treeProduct1.clearSelection();
                                    treeProduct2.clearSelection();
                                    clearReactionsTable();

                                } else{
                                    lblMoleculeNameEmptyError.setText("Can't update. Rotational diffusion missing.");
                                }
                            } else{
                                lblMoleculeNameEmptyError.setText("Can't update. Translational diffusion missing.");
                            }
                        } else{
                            lblMoleculeNameEmptyError.setText("Can't update. Molecule coordinates missing.");
                        }

                    } else{
                        lblMoleculeNameEmptyError.setText("Can't update. Molecule name should not be left blank.");
                    }

                }else{
                    lblMoleculeNameEmptyError.setText("Can't update. Select a molecule.");
                }
            }
        });

        tableMolecule.addMouseListener(new MouseAdapter()
        {
            @Override
            public void mouseClicked(MouseEvent arg0)
            {
                Double x, y, z;
                int numstates;
                String intstates;
                DefaultTableModel model = (DefaultTableModel) tableMolecule.getModel();
                DefaultTableModel modelinterface = (DefaultTableModel) tableInterface.getModel();

                if(tableMolecule.getSelectedRow()>-1){
                    textFieldMoleculeName.setText(model.getValueAt(tableMolecule.getSelectedRow(), 0).toString());
                    textFieldMoleculeCount.setText(model.getValueAt(tableMolecule.getSelectedRow(), 1).toString());
                    if(model.getValueAt(tableMolecule.getSelectedRow(), 2).toString()=="Yes"){
                        chckbxAnchoredToMembrane.setSelected(true);
                    }else{
                        chckbxAnchoredToMembrane.setSelected(false);
                    }
                    textFieldTransDifCoeff.setText(model.getValueAt(tableMolecule.getSelectedRow(), 3).toString());
                    textFieldRotDifCoeff.setText(model.getValueAt(tableMolecule.getSelectedRow(), 4).toString());

                    int numinterfaces = MolCla.get(tableMolecule.getSelectedRow()).getnumInterface();
                    while(tableInterface.getRowCount()>1){
                        modelinterface.removeRow(tableInterface.getRowCount()-1);
                    }
                    for(int i=0; i<numinterfaces;i++){
                        x = MolCla.get(tableMolecule.getSelectedRow()).getinterfaceCoords(i, 0);
                        y = MolCla.get(tableMolecule.getSelectedRow()).getinterfaceCoords(i, 1);
                        z = MolCla.get(tableMolecule.getSelectedRow()).getinterfaceCoords(i, 2);
                        numstates = MolCla.get(tableMolecule.getSelectedRow()).getnumInterfaceStates(i);
                        intstates = MolCla.get(tableMolecule.getSelectedRow()).getinterfaceStateNamesConcat(i);
                        modelinterface.addRow(new Object[]{MolCla.get(tableMolecule.getSelectedRow()).getinterfaceNames(i),x,y,z,numstates,intstates});
                    }
                }
            }
        });

        btnDeleteMolecule.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent e)
            {
                lblMoleculeNameEmptyError.setText("");
                DefaultTableModel model = (DefaultTableModel) tableMolecule.getModel();
                DefaultTableModel modelinterface = (DefaultTableModel) tableInterface.getModel();
                DefaultTreeModel modelreactant1 = (DefaultTreeModel) treeReactant1.getModel();
                DefaultTreeModel modelreactant2 = (DefaultTreeModel) treeReactant2.getModel();

                if(tableMolecule.getSelectedRow()!=-1)
                {
                    if(model.getRowCount()>0){
                        MolCla.remove(tableMolecule.getSelectedRow());
                        model.removeRow(tableMolecule.getSelectedRow());
                        textFieldMoleculeName.setText(null);
                        textFieldMoleculeCount.setText(null);
                        textFieldTransDifCoeff.setText(null);
                        textFieldRotDifCoeff.setText(null);
                        chckbxAnchoredToMembrane.setSelected(false);
                        textFieldInterfaceName.setText(null);
                        textFieldInterfaceX.setText(null);
                        textFieldInterfaceY.setText(null);
                        textFieldInterfaceZ.setText(null);
                        while(tableInterface.getRowCount()>1){
                            modelinterface.removeRow(tableInterface.getRowCount()-1);
                        }
                        tableInterface.clearSelection();
                        tableMolecule.clearSelection();

//						for(int i=0;i<MolCla.size();i++){
//							MolCla.get(i).printMolecule();
//						}

                        fillinNodes(MolCla,modelreactant1);
                        fillinNodes(MolCla,modelreactant2);
                        spinnerNumReactant.setValue(0);
                        spinnerNumProduct.setValue(0);
                        treeReactant1.clearSelection();
                        treeReactant2.clearSelection();
                        treeProduct1.clearSelection();
                        treeProduct2.clearSelection();
                        clearReactionsTable();
                        cleanReactionTab();
                    }else{
                        lblMoleculeNameEmptyError.setText("Add a molecule first."); //Can't change center of mass.
                    }

                }else{
                    lblMoleculeNameEmptyError.setText("Select a molecule.");
                }
            }
        });

        chckbxReversibleReaction.addActionListener(new ActionListener()
        {
            public void actionPerformed(ActionEvent arg0) {
                if(chckbxReversibleReaction.isSelected()){
                    textFieldBackwardRate.setEnabled(true);
                }else{
                    textFieldBackwardRate.setEnabled(false);
                }
                changerateunits((Integer) spinnerNumReactant.getValue(), (Integer) spinnerNumProduct.getValue());
            }
        });

        chckbx2dReaction.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                changerateunits((Integer) spinnerNumReactant.getValue(), (Integer) spinnerNumProduct.getValue());
            }
        });

        spinnerNumReactant.addChangeListener(new ChangeListener()
        {
            public void stateChanged(ChangeEvent arg0)
            {

                DefaultTreeModel modelproduct1 = (DefaultTreeModel) treeProduct1.getModel();
                DefaultTreeModel modelproduct2 = (DefaultTreeModel) treeProduct2.getModel();

                if((Integer) spinnerNumReactant.getValue()==0){
                    treeReactant1.setEnabled(false);
                    treeReactant2.setEnabled(false);
                    generateProductforTrees(modelproduct1, modelproduct2);
                    activateInterfaceDesigner(0);
                    changerateunits((Integer) spinnerNumReactant.getValue(), (Integer) spinnerNumProduct.getValue());
                }else if((Integer) spinnerNumReactant.getValue()==1){
                    treeReactant1.setEnabled(true);
                    treeReactant2.setEnabled(false);
                    generateProductforTrees(modelproduct1, modelproduct2);
                    activateInterfaceDesigner(0);
                    changerateunits((Integer) spinnerNumReactant.getValue(), (Integer) spinnerNumProduct.getValue());
                }else{//2
                    treeReactant1.setEnabled(true);
                    treeReactant2.setEnabled(true);
                    generateProductforTrees(modelproduct1, modelproduct2);
                    if((Integer) spinnerNumProduct.getValue()==1){
                        activateInterfaceDesigner(1);
//					}else{
//						lblSigmanm.setEnabled(true);
//						textFieldSigmaX.setEnabled(true);
                    }
                    changerateunits((Integer) spinnerNumReactant.getValue(), (Integer) spinnerNumProduct.getValue());
                }
            }
        });

        spinnerNumProduct.addChangeListener(new ChangeListener()
        {
            public void stateChanged(ChangeEvent arg0)
            {

                DefaultTreeModel modelproduct1 = (DefaultTreeModel) treeProduct1.getModel();
                DefaultTreeModel modelproduct2 = (DefaultTreeModel) treeProduct2.getModel();

                if((Integer) spinnerNumProduct.getValue()==0)
                {
                    flushTree(modelproduct1);
                    flushTree(modelproduct2);
                    treeProduct1.setEnabled(false);
                    treeProduct2.setEnabled(false);
                    activateInterfaceDesigner(0);
//					//not possible A+B-> 0
//					if((Integer) spinnerNumReactant.getValue()==2){
//						lblSigmanm.setEnabled(true);
//						textFieldSigmaX.setEnabled(true);
//					}
//					changerateunits((Integer) spinnerNumReactant.getValue(), (Integer) spinnerNumProduct.getValue());
                }else if((Integer) spinnerNumProduct.getValue()==1)
                {
                    flushTree(modelproduct2);
                    treeProduct1.setEnabled(true);
                    treeProduct2.setEnabled(false);
                    generateProductforTrees(modelproduct1, modelproduct2);
                    if((Integer) spinnerNumReactant.getValue()==2){
//						activateInterfaceDesigner(1);
                    }
                    changerateunits((Integer) spinnerNumReactant.getValue(), (Integer) spinnerNumProduct.getValue());
                }else
                {//2
                    treeProduct1.setEnabled(true);
                    treeProduct2.setEnabled(true);
                    generateProductforTrees(modelproduct1, modelproduct2);
                    activateInterfaceDesigner(0);
                    if((Integer) spinnerNumReactant.getValue()==2){
                        lblSigmanm.setEnabled(true);
                        textFieldSigmaX.setEnabled(true);
                    }
                    changerateunits((Integer) spinnerNumReactant.getValue(), (Integer) spinnerNumProduct.getValue());
                }
            }
        });

        treeReactant1.addTreeSelectionListener(new TreeSelectionListener()
        {
            public void valueChanged(TreeSelectionEvent arg0)
            {

                DefaultTreeModel modelproduct1 = (DefaultTreeModel) treeProduct1.getModel();
                DefaultTreeModel modelproduct2 = (DefaultTreeModel) treeProduct2.getModel();

                for(int i=0;i<3;i++){
                    react1data[i]=-1;//-1 if NULL, initialization
                }
                for(int i=0;i<6;i++){
                    productdata[i]=-1;//-1 if NULL, initialization
                }
                lblErrorReaction.setText("");
                int level = digitizeTree(modelproduct1, modelproduct2, 1, true);
                panelDrawInterface.repaint();
                if(level<2){
                    activateInterfaceDesigner(0);
                }

                //System.out.println(react1data[0]+" "+react1data[1]+" "+react1data[2]);
            }
        });

        treeReactant2.addTreeSelectionListener(new TreeSelectionListener()
        {
            public void valueChanged(TreeSelectionEvent arg0)
            {
                DefaultTreeModel modelproduct1 = (DefaultTreeModel) treeProduct1.getModel();
                DefaultTreeModel modelproduct2 = (DefaultTreeModel) treeProduct2.getModel();

                for(int i=0;i<3;i++){
                    react2data[i]=-1;
                }
                for(int i=0;i<6;i++){
                    productdata[i]=-1;//-1 if NULL, initialization
                }
                lblErrorReaction.setText("");

                int level = digitizeTree(modelproduct1, modelproduct2, 2, true);
                panelDrawInterface.repaint();
                if(level<2){
                    activateInterfaceDesigner(0);
                }
                //System.out.println(level);
            }
        });

        treeProduct1.addTreeSelectionListener(new TreeSelectionListener()
        {
            public void valueChanged(TreeSelectionEvent arg0)
            {
                digitizeTreeComplex(1);
            }
        });

        treeProduct2.addTreeSelectionListener(new TreeSelectionListener()
        {
            public void valueChanged(TreeSelectionEvent arg0)
            {
                digitizeTreeComplex(2);
            }
        });

        btnAddReaction.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {

                int isReactionConsistent = 1, isDimensionConsistent = 1, isAngleConsistent=1, isSigmaConsistent=1;
                int numberofproducts = Integer.parseInt(spinnerNumProduct.getValue().toString());
                int numberofreactants = Integer.parseInt(spinnerNumReactant.getValue().toString());
                int[] reacvec = new int[6];
                int[] prodvec = new int[6];
                double[] anglesvec = new double[5];
                double[] sigmaD = new double[3];
                double onrateD=0.0, offrateD=0.0;
                boolean onMembrane = false, isRevers = false;
                String theta = "-", thetab = "-", phi = "-", phib = "-", omega="-", sigma="-", offrate="0.0", offrateunit="-"; //default strings
                String sigmaX, sigmaY, sigmaZ;
                String onrate = "0.0", onrateunit = "uM-1/s", rxn = "-";
                lblErrorReaction.setText("");
                DefaultTableModel model = (DefaultTableModel) tableReaction.getModel();



//				chckbxReversibleReaction, chckbx2dReaction,spinnerNumReactant,spinnerNumProduct,SpinnerTheta,SpinnerPhi,SpinnerOmega,textFieldSigma

                for(int i=0;i<3;i++){
                    reacvec[i] = react1data[i];
                    reacvec[i+3] = react2data[i];
                    prodvec[i] = productdata[i];
                    prodvec[i+3] = productdata[i+3];
                    sigmaD[i] = 0.0;
                }

                isAngleConsistent = checkConsistentAngles();

                //check if new reaction is allowed and generate array of set angles (anglesvec)
                if(isAngleConsistent==1 && numberofproducts==1 && numberofreactants==2){
                    theta = Double.toString(sys1.thetaOutput());
                    thetab = Double.toString(sys2.thetaOutput());
                    phi = Double.toString(sys1.phiOutput(sys2));
                    phib = Double.toString(sys2.phiOutput(sys1));
                    omega = Double.toString(sys1.omegaOutput(sys2));
                    anglesvec[0]=Double.parseDouble(theta);
                    anglesvec[1]=Double.parseDouble(thetab);
                    anglesvec[2]=Double.parseDouble(phi);
                    anglesvec[3]=Double.parseDouble(phib);
                    anglesvec[4]=Double.parseDouble(omega);
                }
                isSigmaConsistent = checkConsistentSigma();

                isReactionConsistent = checkConsistentReactionSelection();//everything selected correctly?
                if(isReactionConsistent==1){
                    rxn = stringizeReaction(MolCla);
                }
                if(chckbx2dReaction.isSelected()){//everything selected correctly? (if 2D, then everybody's in 2D?)
                    isDimensionConsistent = checkConsistentReactionDimension();
                    onMembrane = true;
                }

                if(isReactionConsistent==1 && isDimensionConsistent==1)
                {
                    String[] splitArrF1 = lblForwardRate.getText().split("\\("); //lblForwardRate.setText("On Rate (nm-2/us)");
                    String[] splitArrF2 = splitArrF1[1].split("\\)");
                    onrateunit = splitArrF2[0];
                    onrate = textFieldForwardRate.getText();
                    if(NumberUtils.isCreatable(onrate))
                    {
                        onrateD = Double.parseDouble(onrate);
                    }else{
                        onrate = "-";
                        lblErrorReaction.setText("On rate non-numeric.");
                    }

                    sigmaX = textFieldSigmaX.getText();
                    sigma = sigmaX;
                    if(isSigmaConsistent==1 && numberofreactants>1)
                    {
                        sigmaD[0] = Double.parseDouble(sigmaX);
                        sigmaD[1] = 0;
                        sigmaD[2] = 0;
                    }else{
                        sigma = "-";
                        if(isSigmaConsistent==0){
                            lblErrorReaction.setText("Sigma non-numeric.");
                        }
                    }

                    if(chckbxReversibleReaction.isSelected())
                    { //Reversible reaction, need off rate & offrate unit
                        String[] splitArrB1 = lblBackwardRate.getText().split("\\("); //lblBackwardRate.setText("Off Rate (uM/s)");//conditional
                        String[] splitArrB2 = splitArrB1[1].split("\\)");
                        offrateunit = splitArrB2[0];
                        offrate = textFieldBackwardRate.getText();
                        if(NumberUtils.isCreatable(offrate))
                        {
                            offrateD = Double.parseDouble(offrate);
                        }else{
                            lblErrorReaction.setText("Off rate non-numeric.");
                            offrate = "0.0";
                        }
                        isRevers = true;
                    }

                    //nproduct,nreactant,sigma,bool membrane,bool revers,int[] reactant,int[]product,double[] angles, onrate, offrate, onrateunit, offrateunit
                    Reax.add(new Reaction(numberofproducts, numberofreactants, sigmaD, onMembrane, isRevers, reacvec, prodvec, anglesvec, sys1.noRotNormal(), sys2.noRotNormal(), onrateD,offrateD,onrateunit,offrateunit));
                    if(isSameReaction(Reax, Reax.size()-1)==0)
                    {
                        model.addRow(new Object[]{rxn, sigma, onrate, onrateunit, offrate, offrateunit, theta, thetab, phi, phib, omega});
                        cleanReactionTab();

                    }else//same reaction present
                    {
                        lblErrorReaction.setText("Can't add. Same reaction already present.");
                        Reax.remove(Reax.size()-1);
                    }

                }

            }
        });

        btnUpdateReaction.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {

                int isReactionConsistent = 1, isDimensionConsistent = 1, isAngleConsistent=1, isSigmaConsistent=1;
                int numberofproducts = Integer.parseInt(spinnerNumProduct.getValue().toString());
                int numberofreactants = Integer.parseInt(spinnerNumReactant.getValue().toString());
                int[] reacvec = new int[6];
                int[] prodvec = new int[6];
                double[] anglesvec = new double[5];
                double[] sigmaD = new double[3];
                double onrateD=0.0, offrateD=0.0;
                boolean onMembrane = false, isRevers = false;
                String sigma = "-", Theta = "-", Thetab = "-", Phi="-", Phib="-", Omega ="-", offrate="0.0", offrateunit="-"; //default strings
                String sigmaX, sigmaY, sigmaZ;
                String onrate = "0.0", onrateunit = "uM-1/s", rxn = "-";
                lblErrorReaction.setText("");
                DefaultTableModel model = (DefaultTableModel) tableReaction.getModel();

//				chckbxReversibleReaction, chckbx2dReaction,spinnerNumReactant,spinnerNumProduct,SpinnerTheta,SpinnerPhi,SpinnerOmega,textFieldSigma

                for(int i=0;i<3;i++){
                    reacvec[i] = react1data[i];
                    reacvec[i+3] = react2data[i];
                    prodvec[i] = productdata[i];
                    prodvec[i+3] = productdata[i+3];
                    sigmaD[i] = 0.0;
                }
                for(int i = 0; i < 5; i++) {
                    anglesvec[i] = 0.0;
                }

                isAngleConsistent = checkConsistentAngles();
                if(isAngleConsistent==1){
                    Theta = Float.toString(sys1.thetaOutput());
                    Thetab = Float.toString(sys2.thetaOutput());
                    Phi = Float.toString(sys1.phiOutput(sys2));
                    Phib = Float.toString(sys2.phiOutput(sys1));
                    Omega = Float.toString(sys1.omegaOutput(sys2));
                    anglesvec[0]= Math.toRadians(Double.parseDouble(Theta));
                    anglesvec[1]= Math.toRadians(Double.parseDouble(Thetab));
                    anglesvec[2]= Math.toRadians(Double.parseDouble(Phi));
                    anglesvec[3]= Math.toRadians(Double.parseDouble(Phib));
                    anglesvec[4]= Math.toRadians(Double.parseDouble(Omega));
                }
                isSigmaConsistent = checkConsistentSigma();
                isReactionConsistent = checkConsistentReactionSelection();//everything selected correctly?
                if(isReactionConsistent==1){
                    rxn = stringizeReaction(MolCla);
                }
                if(chckbx2dReaction.isSelected()){//everything selected correctly? (if 2D, then everybody's in 2D?)
                    isDimensionConsistent = checkConsistentReactionDimension();
                    onMembrane = true;
                }

                if(isReactionConsistent==1 && isDimensionConsistent==1)
                {
                    String[] splitArrF1 = lblForwardRate.getText().split("\\("); //lblForwardRate.setText("On Rate (nm-2/us)");
                    String[] splitArrF2 = splitArrF1[1].split("\\)");
                    onrateunit = splitArrF2[0];
                    onrate = textFieldForwardRate.getText();
                    if(NumberUtils.isCreatable(onrate))
                    {
                        onrateD = Double.parseDouble(onrate);
                    }else{
                        onrate = "-";
                        lblErrorReaction.setText("On rate non-numeric.");
                    }

                    sigmaX = textFieldSigmaX.getText();
                    sigma = sigmaX;
                    if(isSigmaConsistent==1 && numberofreactants>1)
                    {
                        sigmaD[0] = Double.parseDouble(sigmaX);
                        sigmaD[1] = 0;
                        sigmaD[2] = 0;
                    }else{
                        sigma = "-";
                        if(isSigmaConsistent==0){
                            lblErrorReaction.setText("Sigma non-numeric.");
                        }
                    }

                    if(chckbxReversibleReaction.isSelected())
                    { //Reversible reaction, need off rate & offrate unit
                        String[] splitArrB1 = lblBackwardRate.getText().split("\\("); //lblBackwardRate.setText("Off Rate (uM/s)");//conditional
                        String[] splitArrB2 = splitArrB1[1].split("\\)");
                        offrateunit = splitArrB2[0];
                        offrate = textFieldBackwardRate.getText();
                        if(NumberUtils.isCreatable(offrate))
                        {
                            offrateD = Double.parseDouble(offrate);
                        }else{
                            lblErrorReaction.setText("Off rate non-numeric.");
                            offrate = "0.0";
                        }
                        isRevers = true;
                    }

                    if(tableReaction.getSelectedRow()!=-1)
                    {
                        Reax.set(tableReaction.getSelectedRow(), new Reaction(numberofproducts, numberofreactants, sigmaD, onMembrane, isRevers, reacvec, prodvec, anglesvec, sys1.noRotNormal(), sys2.noRotNormal(), onrateD,offrateD,onrateunit,offrateunit));

                        if(isSameReaction(Reax, tableReaction.getSelectedRow())==0)
                        {

                            model.setValueAt(rxn, tableReaction.getSelectedRow(), 0);
                            model.setValueAt(sigma, tableReaction.getSelectedRow(), 1);
                            model.setValueAt(onrate, tableReaction.getSelectedRow(), 2);
                            model.setValueAt(onrateunit, tableReaction.getSelectedRow(), 3);
                            model.setValueAt(offrate, tableReaction.getSelectedRow(), 4);
                            model.setValueAt(offrateunit, tableReaction.getSelectedRow(), 5);
                            model.setValueAt(Theta, tableReaction.getSelectedRow(), 6);
                            model.setValueAt(Thetab, tableReaction.getSelectedRow(), 7);
                            model.setValueAt(Phi, tableReaction.getSelectedRow(), 8);
                            model.setValueAt(Phib, tableReaction.getSelectedRow(), 9);
                            model.setValueAt(Omega, tableReaction.getSelectedRow(), 10);


                        }else//same molecule present
                        {
                            lblErrorReaction.setText("Can't update. Selected reaction already present and will be deleted.");
                            Reax.remove(tableReaction.getSelectedRow());
                            model.removeRow(tableReaction.getSelectedRow());
                        }

                        cleanReactionTab();

                    }else{
                        lblErrorReaction.setText("Can't update. Select a reaction.");
                    }

                }

            }
        });

        btnDeleteReaction.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {

                lblErrorReaction.setText("");
                DefaultTableModel model = (DefaultTableModel) tableReaction.getModel();

                if(tableReaction.getSelectedRow()!=-1)
                {
                    if(model.getRowCount()>0){
                        Reax.remove(tableReaction.getSelectedRow());
                        model.removeRow(tableReaction.getSelectedRow());
                        cleanReactionTab();

                    }else{
                        lblErrorReaction.setText("Add a reaction first."); //Can't change center of mass.
                    }

                }else{
                    lblErrorReaction.setText("Select a reaction.");
                }

            }
        });

        tableReaction.addMouseListener(new MouseAdapter() {
            @Override
            public void mouseClicked(MouseEvent arg0) {

                DefaultTreeModel modelproduct1 = (DefaultTreeModel) treeProduct1.getModel();
                DefaultTreeModel modelproduct2 = (DefaultTreeModel) treeProduct2.getModel();
                int therow = tableReaction.getSelectedRow();

                if(therow>-1){
                    int[] tempvecR = new int[6];
                    int[] tempvecP = new int[6];
                    int[] p1vec = new int[3];
                    int[] p2vec = new int[3];

                    tempvecR = Reax.get(therow).getReactant();
                    tempvecP = Reax.get(therow).getProduct();
                    //				System.out.println(Reax.get(therow).getNumProduct());
                    //				Reax.get(therow).printReaction();
                    spinnerNumProduct.setValue(Reax.get(therow).getNumProduct());
                    spinnerNumReactant.setValue(Reax.get(therow).getNumReactant());

                    for(int i=0;i<3;i++){
                        react1data[i]=  tempvecR[i];
                        react2data[i]= tempvecR[i+3];
                        productdata[i] = tempvecP[i];
                        productdata[i+3] = tempvecP[i+3];
                        p1vec[i] = tempvecP[i];
                        p2vec[i] = tempvecP[i+3];
                    }
                    generateProductforTrees(modelproduct1, modelproduct2);

                    if(Reax.get(therow).getNumProduct()==1){

                        expandAllNodes(treeProduct1, 0, treeProduct1.getRowCount());
                        treeProduct1.setSelectionRow(getRowTree(p1vec));

                    }else if(Reax.get(therow).getNumProduct()==2){

                        expandAllNodes(treeProduct1, 0, treeProduct1.getRowCount());
                        expandAllNodes(treeProduct2, 0, treeProduct2.getRowCount());
                        treeProduct1.setSelectionRow(getRowTree(p1vec));
                        treeProduct2.setSelectionRow(getRowTree(p2vec));
                    }

                    if(Reax.get(therow).getNumReactant()==1){

                        expandAllNodes(treeReactant1, 0, treeReactant1.getRowCount());
                        treeReactant1.setSelectionRow(getRowTree(react1data));

                    }else if(Reax.get(therow).getNumReactant()==2){

                        expandAllNodes(treeReactant1, 0, treeReactant1.getRowCount());
                        expandAllNodes(treeReactant2, 0, treeReactant2.getRowCount());
                        treeReactant1.setSelectionRow(getRowTree(react1data));
                        treeReactant2.setSelectionRow(getRowTree(react2data));
                    }

                    if(Reax.get(therow).getmolOnMembrane()==true){
                        chckbx2dReaction.setSelected(true);
                    }
                    if(Reax.get(therow).getRever()==true){
                        chckbxReversibleReaction.setSelected(true);
                    }

                    textFieldSigmaX.setText(Double.toString(Reax.get(therow).getSigmaIndiv(0)));
                    SpinnerTheta.setValue(Reax.get(therow).getAngleIndiv(0));
                    SpinnerThetab.setValue(Reax.get(therow).getAngleIndiv(1));
                    SpinnerPhi.setValue(Reax.get(therow).getAngleIndiv(2));
                    SpinnerPhib.setValue(Reax.get(therow).getAngleIndiv(3));
                    SpinnerOmega.setValue(Reax.get(therow).getAngleIndiv(4));
                }
            }
        });

        SpinnerTheta.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent arg0) {
                panelDrawInterface.repaint();
            }
        });

        SpinnerThetab.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent arg0) {
                panelDrawInterface.repaint();
            }
        });

        SpinnerPhi.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent arg0) {
                panelDrawInterface.repaint();
            }
        });

        SpinnerPhib.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent arg0) {
                panelDrawInterface.repaint();
            }
        });

        SpinnerOmega.addChangeListener(new ChangeListener() {
            public void stateChanged(ChangeEvent arg0) {
                panelDrawInterface.repaint();
            }
        });

        textFieldSigmaX.addCaretListener(new CaretListener() {
            public void caretUpdate(CaretEvent arg0) {
                panelDrawInterface.repaint();
            }
        });

        panelMolecule.addMouseWheelListener(new MouseWheelListener() {
            public void mouseWheelMoved(MouseWheelEvent arg0) {
                double[] limits = new double[4];
                double delta = 0.2f * arg0.getPreciseWheelRotation();
                limits = camera.getLimits();
                camera.setScale(limits[1]+delta);
                panelMolecule.repaint();
            }
        });

        panelDrawInterface.addMouseWheelListener(new MouseWheelListener() {
            public void mouseWheelMoved(MouseWheelEvent arg0) {
                double[] limits = new double[4];
                double delta = 0.2f * arg0.getPreciseWheelRotation();
                limits = cameraReaction.getLimits();
                cameraReaction.setScale(limits[1]+delta);
                panelDrawInterface.repaint();
            }
        });

    }
}
