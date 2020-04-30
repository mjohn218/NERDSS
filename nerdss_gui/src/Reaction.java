import javax.vecmath.*;

public class Reaction
{
	String onrateunit,offrateunit; //Molecule name
	int nreactant; //Number of reactants
	int nproduct;//Number of products
	double onrate, offrate; //Reaction sigma
	boolean membrane; //true if molecule(s) anchored to the membrane
	boolean reversible;//true if reversible
	//reaction vectors
	double[] angles, sigma;
	int[] reactant; //6x1 vec [012][345] reactant1 reactant2
	int[] product; //6x1 vec [012][345] product1+product2 summary
	Vector3f norm1;
	Vector3f norm2;

	// This is the constructor of the class Molecule
	public Reaction(int nproduct, int nreactant, double[] sigma, boolean membrane, boolean reversible, int[] reactant, int[] product, double[] angles, Vector3f norm1, Vector3f norm2, double onrate, double offrate,String onrateunit, String offrateunit)
	{
		this.nproduct = nproduct;
		this.nreactant = nreactant;
		this.onrateunit = onrateunit;
		this.offrateunit = offrateunit;
		this.membrane = membrane;
		this.reversible = reversible;
		this.onrate = onrate;
		this.offrate = offrate;
        this.reactant = new int[6];
        this.product  = new int[6];
        this.angles = new double[5];
        this.sigma = new double[3];
        this.norm1 = norm1;
        this.norm2 = norm2;
        
        for (int i = 0; i < 5; i++)
        {
        	this.angles[i] = angles[i];
        }
        for (int i = 0; i < 3; i++) {
			this.sigma[i] = sigma[i];
		}
        for (int i = 0; i < 6; i++)
        {
        	this.reactant[i] = reactant[i];
        	this.product[i] = product[i];
        }
	}

	// Assign/Retrieve the number of reactants
	public void setNumReactant(int nreactants) 
	{
		nreactant = nreactants;
	}
	public int getNumReactant() 
	{
		return nreactant;
	}
	
	// Assign/Retrieve the number of products
	public void setNumProduct(int nproducts) 
	{
		nproduct = nproducts;
	}
	public int getNumProduct() 
	{
		return nproduct;
	}
	
	// Assign/Retrieve sigma
	public double[] getSigma() 
	{
		return sigma;
	}
	public double getSigmaIndiv(int i) 
	{
		return sigma[i];
	}	
	public void setSigma(double[] sigmaz) 
	{
        for (int i = 0; i < 3; i++)
        {
        	this.sigma[i] = sigmaz[i];
        }
	}
	
	// Assign/Retrieve the molOnMembrane to the variable membrane.*/
	public void setmolOnMembrane(boolean molOnMembrane) 
	{
		membrane = molOnMembrane;
	}
	public boolean getmolOnMembrane() 
	{
		return membrane;
	}
	
	// Assign/Retrieve the reversiblilty.*/
	public void setRever(boolean isrever) 
	{
		reversible = isrever;
	}
	public boolean getRever() 
	{
		return reversible;
	}
		
	// Assign/Retrieve the reactant 6x1 vec
	public int[] getReactant() 
	{
		return reactant;
	}
	public int getReactantIndiv(int i) 
	{
		return reactant[i];
	}	
	public void setReactant(int[] reactants) 
	{
        for (int i = 0; i < 6; i++)
        {
        	this.reactant[i] = reactants[i];
        }
	}	
	
	// Assign/Retrieve the product 6x1 vec
	public int[] getProduct() 
	{
		return product;
	}
	public int getProductIndiv(int i) 
	{
		return product[i];
	}	
	public void setProduct(int[] products) 
	{
        for (int i = 0; i < 6; i++)
        {
        	this.product[i] = products[i];
        }
	}

	// Assign/Retrieve the product 6x1 vec
	public double[] getAngle() 
	{
		return angles;
	}
	public double getAngleIndiv(int i) 
	{
		return angles[i];
	}	
	public void setAngle(double[] anglesz) 
	{
        for (int i = 0; i < 5; i++)
        {
        	this.angles[i] = anglesz[i];
        }
	}	
	
	/* Print the Molecule details */
	public void printReaction() 
	{
		System.out.println("Number of Products: " + nproduct);
		System.out.println("Number of Reactants: " + nreactant);
		System.out.println("Anchored to membrane: " + membrane);
		System.out.println("Reversible?: " + reversible);
		System.out.println("Sigma: " + sigma);
        
		System.out.println("Angles:" + angles[0] +" "+ angles[1] +" "+ angles[2]);
		for (int k = 0; k < nproduct; k++){
        	System.out.println("Products:" + product[0+3*k] +" "+ product[1+3*k] +" "+ product[2+3*k]);
        }
		for (int k = 0; k < nreactant; k++){
        	System.out.println("Reactants:" + reactant[0+3*k] +" "+ reactant[1+3*k] +" "+ reactant[2+3*k]);
        }
	}
}