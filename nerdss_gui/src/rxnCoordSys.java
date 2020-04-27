import sun.font.TrueTypeFont;

import javax.vecmath.*;
import java.util.ArrayList;

/**
 * A class to compute dependent rotations. Uses a quarterion to transform the interface vector.
 * Needed because openGL is not flexible enough to allow for this sort of rotation easily.
 */
public class rxnCoordSys {
    public Vector3f intVec;
    public ArrayList<Vector3f> init_secondaryInterfaces;
    public ArrayList<Vector3f> cur_secondaryInterfaces;
    public ArrayList<Float> secondaryScales;
    public float scale;
    public Vector3f curVec;
    public Vector3f curPhi;
    public float phiAngle;
    public int interfaceNum;
    public float intLength;
    private Vector3f sigma;
    private Matrix3f rotMatrix;
    private Quat4f rotQuat;
    private Quat4f phiQuat;
    private float oldTheta;
    public float oldPhi;
    private float oldOmega;
    private boolean primary;
    public boolean Tneg;
    public boolean Oneg;
    public boolean Pneg;



    /**
     * construct a new coordinate system.
     * @param bindSite the interface vector
     * @param sigma the spacing between interfaces
     * @param interfaceNum for interface situations
     */
    public rxnCoordSys(Vector3f bindSite, Vector3f sigma, double[][] secondaryInt, int interfaceNum, boolean primary) {
        this.rotMatrix = new Matrix3f();
        this.rotMatrix.setIdentity();
        this.rotQuat = new Quat4f(0,0,0,1);
        this.phiQuat = new Quat4f(0,0,0,1);
        this.sigma = new Vector3f(sigma);
        this.scale = bindSite.length();
        bindSite.normalize();
        this.intVec = new Vector3f(bindSite);
        this.primary = primary;
        if (!primary) {
            this.intVec.scale(-1);
        }
        this.curVec = new Vector3f(this.intVec);
        this.intLength = this.intVec.length();
        this.curPhi = new Vector3f(this.getOrigNormal());
        this.oldOmega = -999;
        this.oldPhi = -999;
        this.oldTheta = -999;
        this.interfaceNum = interfaceNum;

        this.init_secondaryInterfaces = new ArrayList<>();
        this.cur_secondaryInterfaces = new ArrayList<>();
        this.secondaryScales = new ArrayList<>();
        this.Tneg = false;
        this.Oneg = false;
        this.Pneg = false;

        // Below is to add all secondary interfaces so we can properly transform them
        if (secondaryInt.length > 1) {
            for (int i = 0; i < secondaryInt.length; i++) {
                if (i != this.interfaceNum) {
                    Vector3f temp = new Vector3f((float) secondaryInt[i][0], (float) secondaryInt[i][1], (float) secondaryInt[i][2]);
                    this.secondaryScales.add(temp.length());
                    temp.normalize();
                    if (!this.primary) {
                        temp.scale(-1);
                    }
                    this.init_secondaryInterfaces.add(temp);
                }
            }
            this.cur_secondaryInterfaces = (ArrayList<Vector3f>) this.init_secondaryInterfaces.clone();
        }
    }

    /**
     * Generate rotation quatarion for requested rotation.
     * @param a the degrees of rotation (RAD)
     * @param v the axis of rotation.
     * @return a rotation quat
     */
    private Quat4f rotationAboutAxis(float a, Vector3f v) {
        Quat4f R = new Quat4f();
        v.normalize();

        AxisAngle4f axisAngle = new AxisAngle4f(v.x, v.y, v.z, a);

        R.set(axisAngle);

        return R;
    }

    public Vector3f noRotNormal() {
        Vector3f temp = new Vector3f(0,0,1);
        return temp;
    }

    /**
     * get the original normal vector prior to any rotations
     * @returnVector3f normal
     */
    public Vector3f getOrigNormal() {
        Vector3f temp = this.projectToPlane(this.sigma, this.curVec);
        temp.normalize();
        return temp;
    }

    /**
     * get the current normal vector
     * @return Vector3f normal
     */
    public Vector3f getNormal(){
      return this.curPhi;
    }

    /**
     * project a vector onto plane defined by normal
     * @param vector1 the vector to project.
     * @param normal vector orthogonal to the plane to project on
     * @return the projection vector
     */
    public Vector3f projectToPlane(Vector3f vector1, Vector3f normal) {
        Vector3f vec = new Vector3f(vector1);
        Vector3f norm = new Vector3f(normal);
        float dot = vec.dot(norm);
        dot = dot/norm.lengthSquared();
        norm.scale(-1*dot);
        vec.add(norm);
        return vec;
    }

    /**
     * get the current theta angle after rotations
     * @return
     */
    public float curThetaAngle() {
        Vector3f BindSite = curVec;
        BindSite.normalize();
        Vector3f sigma = new Vector3f(1,0,0);

        float theta = BindSite.angle(sigma);
        if (this.Tneg) {
            theta *= -1;
        }
        return theta;
    }

    /**
     * return the proper omega for output files
     * @return
     */
    public float thetaOutput() {
        float trueTheta= curThetaAngle();
        if (this.almostEqual(trueTheta, 0, (float) .017)) {
            trueTheta = 0;
        } else if (this.almostEqual(trueTheta, (float) Math.PI, 1)) {
            trueTheta = 180;
        }
        return trueTheta;
    }

    /**
     * current phi angle after rotations
     * @return
     */
   public float curPhiAngle() {
        Vector3f norm = getNormal();
        norm.normalize();
        Vector3f sig_proj = this.projectToPlane(this.sigma, this.curVec);
        float phi = norm.angle(sig_proj);
        if (this.Pneg) {
            phi *= -1;
        }
       return phi;
    }

    public float phiOutput(rxnCoordSys sys) {
        float truePhi= curPhiAngle();
        if (this.almostEqual(truePhi, -1*(float) Math.PI, (float) .017)) {
            truePhi = -180;
        } else if (this.almostEqual(truePhi, (float) Math.PI, (float) .017)) {
            truePhi = 180;
        } else if (this.thetaOutput() == 0 || this.thetaOutput() == 180 || sys.thetaOutput() == 0 || sys.thetaOutput() == 180) {
            //illogical value -999 used to indicate DNE
            truePhi = -999;
        }
        //compensate for different angle conventions in nerdss
        if (truePhi < 0) {
            truePhi *= -1;
            truePhi = 2* (float) Math.PI - truePhi;
        }
        return truePhi;
    }


    /**
     * current omega angle after rotations
     * @return
     */
    public float curOmegaAngle(rxnCoordSys sys) {
        Vector3f BindSite1 = this.curVec;
        float omega = (projectToPlane(BindSite1, this.sigma).angle(projectToPlane(sys.curVec, this.sigma)));
        if (this.Oneg) {
            omega *= -1;
        }
        return omega;
    }

    /**
     * return the proper omega for output files. if theta is 0 or 180 for either molecule define omega as the angle
     * between the components of the normal vectors orthogonal to sigma in this case omega change is locked, only phi
     * change is allowed to produce effect (since here phi and omega have essentially the same meaning.
     * @param sys
     * @return
     */
    public float omegaOutput(rxnCoordSys sys) {
        float trueOmega = curOmegaAngle(sys);
        if (this.almostEqual(trueOmega, -1*(float) Math.PI, (float) .017)) {
            trueOmega = 0;
        } else if (this.almostEqual(trueOmega, (float) Math.PI, (float) .017)) {
            trueOmega = 180;
        }
        if (this.thetaOutput() == 0 || this.thetaOutput() == 180 || sys.thetaOutput() == 0 || sys.thetaOutput() == 180) {
            trueOmega = (projectToPlane(this.curPhi, this.sigma).angle(projectToPlane(sys.curPhi, this.sigma)));
        }
        //compensate for different angle conventions in nerdss
        if (trueOmega < 0) {
            trueOmega *= -1;
            trueOmega = 2* (float) Math.PI - trueOmega;
        }
        return trueOmega;
    }

    /**
     * return current phi axis.
     * @return
     */
    public Vector3f phiRotationAxis() {
        return this.curVec;
    }


    /**
     * Calculate the current theta axis
     * @return the vector for the axis
     */
    public Vector3f thetaRotationAxis() {
        Vector3f thetaAxis = new Vector3f();
        thetaAxis.cross(this.curVec, this.sigma);

        return thetaAxis;
    }

    /**
     * return omega axis of rotation equal to sigma
     * @return
     */
    public Vector3f omegaAxis() {
        return this.sigma;
    }

    public AxisAngle4f getSingleRot() {
        AxisAngle4f totalRot = new AxisAngle4f();
        totalRot.set(rotQuat);
        return totalRot;
    }

    public Vector3f modifyByGlobalQuat(Vector3f in) {
        Quat4f orig1 = new Quat4f(in.x, in.y, in.z, (float) 0);
        orig1.normalize();
        orig1.mul(this.rotQuat);
        Vector3f out = new Vector3f(orig1.x, orig1.y, orig1.z);
        return out;
    }

    /**
     * update the current interface vector following a requested change.
     * @param axis
     * @param dAngle
     */
    public void updateVector(Vector3f axis, float dAngle) {
        Quat4f orig1 = new Quat4f(this.intVec.x, this.intVec.y, this.intVec.z, (float) 0);
        orig1.normalize();
        this.rotQuat.mul(rotationAboutAxis(dAngle, axis));
        this.rotQuat.normalize();
        orig1.mul(this.rotQuat);
        orig1.normalize();
        Vector3f newVec = new Vector3f(orig1.x, orig1.y, orig1.z);
        this.curVec.set(newVec);
        if (this.init_secondaryInterfaces != null) {
            for (int i = 0; i < this.init_secondaryInterfaces.size(); i++) {
                orig1 = new Quat4f(this.init_secondaryInterfaces.get(i).x,
                        this.init_secondaryInterfaces.get(i).y,
                        this.init_secondaryInterfaces.get(i).z,
                        (float) 0);
                orig1.normalize();
                orig1.mul(this.rotQuat);
                orig1.normalize();
                this.cur_secondaryInterfaces.set(i, new Vector3f(orig1.x, orig1.y, orig1.z));
            }
        }
    }

    /**
     * return the phi modified normal vector.
     * @param Angle angle to rotate by.
     * @return the updated normal vector3f.
     */
    public void phiMod(float Angle) {
        this.phiAngle = Angle;
        Vector3f orig = this.getOrigNormal();
        Quat4f orig1 = new Quat4f(orig.x, orig.y, orig.z, 0);
        orig1.normalize();
        orig1.mul(rotationAboutAxis(Angle, this.curVec));
        orig1.normalize();
        Vector3f orig2 = new Vector3f(orig1.x, orig1.y, orig1.z);
        this.curPhi.set(orig2);
    }

    /**
     * determine if two floats are close enough
     * @param a
     * @param b
     * @param eps
     * @return
     */
    private static boolean almostEqual(float a, float b, float eps){
        return Math.abs(a-b)<eps;
    }

    /**
     * determine if two vectors are close enough.
     * @param a
     * @param b
     * @param eps
     * @return
     */
    private static boolean almostEqual(Vector3f a, Vector3f b, float eps) {
        a.normalize();
        b.normalize();
        return (almostEqual(a.x, b.x, eps) &&
                almostEqual(a.y, b.y, eps) &&
                almostEqual(a.y, b.y, eps));
    }

    /**
     * for testing.
     * @param args
     */
    public static void main(String args[] ){
        /*Vector3f vec1 = new Vector3f(1,1,1);
        Vector3f sigma = new Vector3f(1,0,0);
        rxnCoordSys test = new rxnCoordSys(vec1, sigma, 0);

        float angles[] = {test.curPhiAngle(), test.curThetaAngle()};

        Vector3f normal = test.getNormal();

        assert (almostEqual(angles[0],(float) 2.35, (float).01));
        assert (almostEqual(angles[1],(float) .955, (float).01));

        Vector3f thetaAxis = test.thetaRotationAxis();
        Vector3f phiAxis = test.phiRotationAxis();

        assert (almostEqual(test.phiRotationAxis(), new Vector3f(0,0,-1), (float).01));
        test.updateVector(thetaAxis, (float)(Math.PI / 4));

        //assert the bind site vector is expected after theta transform
        assert (almostEqual(test.curVec, new Vector3f((float).9855985, (float).11957326, (float).11957326), (float) .01));
        phiAxis = test.phiRotationAxis();

        //assert normal is expected after theta rot
        assert (almostEqual(test.getNormal(), new Vector3f((float)-0.207107, (float)1.353553, (float)0.353553), (float) .01));

        test.updateVector(thetaAxis, (float)(Math.PI / 4));*/


    }
}

