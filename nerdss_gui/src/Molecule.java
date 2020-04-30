import java.io.*;
import javax.swing.JTable;
import com.sun.javafx.geom.Vec3f;

public class Molecule 
{

	String name; //Molecule name
	int count; //Molecule count
	int numberofInterfaces;
	double diffT; //Molecule diffusion coeff
	double diffR; //Molecule diffusion coeff
	boolean membrane; //true if anchored to the membrane
	//vector interfaces and interface names
	String[] interfacenames;
	String[] interfaceStatesConcat;
	String[][] interfaceStates;
	int[] numberofInterfaceStates;
	double[][] interfacecoords;// numberofInterfaces-by-3 array
	float [][] interfacenormals; //normal for phi calculations

	// This is the constructor of the class Molecule
	public Molecule(String name, int count, double diffT, double diffR, boolean membrane, int numberofInterfaces, JTable tableInterface) //String[] interfacenames, double[][] interfacecoords
	{
		this.name = name;
		this.count = count;
		this.diffT = diffT;
		this.diffR = diffR;
		this.membrane = membrane;
		this.numberofInterfaces = numberofInterfaces;
        this.interfacecoords = new double[numberofInterfaces][3];
		this.interfacenormals = new float[numberofInterfaces][3];

        this.numberofInterfaceStates = new int[numberofInterfaces];
        this.interfaceStates = new String[numberofInterfaces][10];
        this.interfaceStatesConcat = new String[numberofInterfaces];
        
        for (int i = 0; i < numberofInterfaces; i++){
            for (int j = 0; j < 3; j++){
            		try {
						this.interfacecoords[i][j] = (Double) tableInterface.getValueAt(i+1, j+1);
					} catch (ClassCastException e) {
						this.interfacecoords[i][j] = Double.parseDouble( (String) tableInterface.getValueAt(i+1, j+1));
					}
            }
			this.interfacenormals[i] = getNormal(this.interfacecoords[i]);
            this.numberofInterfaceStates[i] = (int) tableInterface.getValueAt(i+1, 4);
            
            String[] strNamesofStates = tableInterface.getValueAt(i+1, 5).toString().split(",");
            this.interfaceStatesConcat[i] = tableInterface.getValueAt(i+1, 5).toString();
            for(int k=0;k<this.numberofInterfaceStates[i];k++){
//            	System.out.println(strNamesofStates[k].trim());
            	this.interfaceStates[i][k] = strNamesofStates[k].trim();
            }
        }
        this.interfacenames = new String[numberofInterfaces];
        for (int k = 0; k < numberofInterfaces; k++)
            this.interfacenames[k] = (String) tableInterface.getValueAt(k+1, 0);
	}

	private float [] getNormal(double[] coords) throws IndexOutOfBoundsException {
		if (coords.length != 3) {
			throw new IndexOutOfBoundsException("invalid coordinate length");
		}

		Vec3f coord = new Vec3f((float) coords[0], (float) coords[1], (float) coords[2]);
		Vec3f z = new Vec3f(0 , 0, 1);
		Vec3f norm = new Vec3f();
		norm.cross(coord,z);
		float[] temp = {norm.x, norm.y, norm.z};
		return temp;
	}
	
	/* Assign the molecule name.*/
	public void setName(String molname) 
	{
		name = molname;
	}
	public String getName() 
	{
		return name;
	}   

	// Assign the number of molecules of to the variable count.
	public void setmolCount(int molCount) 
	{
		count = molCount;
	}
	public int getmolCount() 
	{
		return count;
	}
	
	// Assign the number of interfaces of to the variable count.
	public void setnumInterface(int numInterface) 
	{
		numberofInterfaces = numInterface;
	}
	public int getnumInterface() 
	{
		return numberofInterfaces;
	}	
	
	public int getnumInterfaceStates(int i) 
	{
		return numberofInterfaceStates[i];
	}
	public String getinterfaceStateNames(int i,int j) 
	{
		return interfaceStates[i][j];
	}
	public String getinterfaceStateNamesConcat(int i) 
	{
		return interfaceStatesConcat[i];
	}

	/* Assign the TransDiff to the variable diffT.*/
	public void setmolDiffT(double molDiffT) 
	{
		diffT = molDiffT;
	}
	public double getmolDiffT() 
	{
		return diffT;
	}

	/* Assign the RotDiff to the variable diffR.*/
	public void setmolDiffR(double molDiffR) 
	{
		diffR = molDiffR;
	}
	public double getmolDiffR() 
	{
		return diffR;
	}

	/* Assign the molOnMembrane to the variable membrane.*/
	public void setmolOnMembrane(boolean molOnMembrane) 
	{
		membrane = molOnMembrane;
	}
	public boolean getmolOnMembrane() 
	{
		return membrane;
	}
	
	public String getinterfaceNames(int i) 
	{
		return interfacenames[i];
	}
	public double getinterfaceCoords(int i, int j) 
	{
		return interfacecoords[i][j];
	}

	public float[] getnormalCoords(int i)
	{
		return interfacenormals[i];
	}
	
	//Total number of interfaces + states + molecule(1)
	public int getBulkNum(){
		
		int bulk = 1;
		
		for(int i = 0; i<numberofInterfaces;i++){
			
			bulk = bulk + 1;
			
			bulk = bulk + numberofInterfaceStates[i]; 
			
		}		
		
		return bulk;
		
	}
	
	/* Print the Molecule details */
	public void printMolecule() 
	{
		System.out.println("Name: "+ name );
		System.out.println("Count: " + count );
		System.out.println("Trans Diffusion coeff: " + diffT );
		System.out.println("Rot Diffusion coeff: " + diffR );
		System.out.println("Anchored to membrane: " + membrane);
		System.out.println("#Interfaces: " + numberofInterfaces );
		System.out.println("Center of mass 0.0 0.0 0.0");
        for (int k = 0; k < numberofInterfaces; k++){
        	System.out.println(interfacenames[k] + " "+ interfacecoords[k][0] +" "+ interfacecoords[k][1] +" "+ interfacecoords[k][2] );
        }
	}
}