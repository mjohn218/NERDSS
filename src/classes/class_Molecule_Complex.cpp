/*! @file class_subs.cpp
 * \brief Class member functions
 *
 * ### Created on 3/29/18 by Matthew Varga
 * ### Purpose
 * ***
 * Class member function definitions for all molecule and associated classes
 *
 * ### TODO List
 * ***
 */

#include "classes/class_Molecule_Complex.hpp"
#include "classes/class_Rxns.hpp"
#include "classes/class_bngl_parser.hpp"
#include "io/io.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "parser/parser_functions.hpp"

#include <cmath>
#include <iomanip>
#include <numeric>

int Molecule::numberOfMolecules = 0;
int Complex::numberOfComplexes = 0;
int Complex::currNumberComTypes = 0;
int Complex::currNumberMolTypes = 0;
std::vector<int> Complex::emptyComList {};
std::vector<int> Molecule::emptyMolList {};
std::vector<int> Complex::obs {};

int propCalled = 0;

bool skipLine(std::string line)
{
    // check if the line is a comment or empty. I can't get regex to work with this in particular, without also
    // including needed lines that have comments at the end
    return line.empty() || line[0] == '#';
}

Complex::Complex(const Molecule& mol, const MolTemplate& oneTemp)
    : comCoord(mol.comCoord)
    , D(oneTemp.D)
    , Dr(oneTemp.Dr)
{
    comCoord = mol.comCoord;
    mass = mol.mass;
    memberList.push_back(mol.index);
    index = mol.index;
    radius = oneTemp.radius;
    // Will elements of this array below be initialized to zero??
    numEachMol = std::vector<int>(MolTemplate::numMolTypes);
    ++numEachMol[oneTemp.molTypeIndex];
}

Complex::Complex(int _index, const Molecule& _memMol, const MolTemplate& _molTemp)
    : comCoord(_memMol.comCoord)
    , D(_molTemp.D)
    , Dr(_molTemp.Dr)
{
    comCoord = _memMol.comCoord;
    mass = _memMol.mass;
    memberList.push_back(_memMol.index);
    index = _index;
    radius = _molTemp.radius;
    // Will elements of this array below be initialized to zero??
    numEachMol = std::vector<int>(MolTemplate::numMolTypes);
    ++numEachMol[_molTemp.molTypeIndex];
}

Complex::Complex(Coord comcoords, Coord D, Coord Dr)
    : comCoord(comcoords)
    , D(D)
    , Dr(Dr)
{
}

// Overloaded operators //
std::ostream& operator<<(std::ostream& os, const std::array<double, 3>& arr)
{
    os << arr[0] << ' ' << arr[1] << ' ' << arr[2];
    return os;
}

void operator+(Coord& c, const double scal)
{
    c.x += scal;
    c.y += scal;
    c.z += scal;
}

void operator*(Coord& c, const double scal)
{
    c.x *= scal;
    c.y *= scal;
    c.z *= scal;
}

std::ostream& operator<<(std::ostream& os, const Molecule& mol)
{
    os << mol.tmpComCoord << std::endl;
    for (auto& iface : mol.tmpICoords)
        os << iface << std::endl;
    return os;
}

// bool operator==(const IntCoordCont::Mol& mol1, const IntCoordCont::Mol& mol2)
//{
//    return (bool{ round(mol1.comCoord) == round(mol2.comCoord) } * bool{ mol1.proType == mol2.proType }
//        * bool{ mol1.reactantIface == mol2.reactantIface } * bool{ mol1.ifaceCoords == mol2.ifaceCoords });
//}

/* MOLECULE::IFACE */
void Molecule::Iface::change_state(int newRelIndex, int newAbsIndex, char newIden)
{
    stateIndex = newRelIndex;
    index = newAbsIndex;
    stateIden = newIden;
}

/* MOLECULE */
void Molecule::display(const MolTemplate& molTemplate) const
{
    std::cout << "Index: " << index << '\n';
    std::cout << "Is empty: " << std::boolalpha << isEmpty << '\n';
    if (!isEmpty) {
        std::cout << "Type: " << molTemplate.molName << '\n';
        std::cout << "Parent complex index: " << myComIndex << '\n';
        std::cout << "Sub volume index: " << mySubVolIndex << '\n';
        std::cout << "Is a lipid: " << std::boolalpha << isLipid << '\n';
        std::cout << "Center of mass coordinate: " << comCoord << '\n';
        std::cout << "Interfaces:\n";
        for (const auto& iface : interfaceList) {
            std::cout << "\t---\n";
            std::cout << "\tRelative index: " << iface.relIndex << '\n';
            std::cout << "\tAbsolute index: " << iface.index << '\n';
            std::cout << "\tInterface name: "
                      << molTemplate.interfaceList[iface.relIndex].name << '\n';
            std::cout << "\tCoordinate: " << iface.coord << '\n';
            if (iface.stateIden != '\0')
                std::cout << "\tCurrent state: " << iface.stateIden << '\n';
            if (iface.isBound) {
                std::cout << "\tInteraction:\n";
                std::cout << "\t\tPartner index: " << iface.interaction.partnerIndex << '\n';
                std::cout << "\t\tPartner interface index " << iface.interaction.partnerIfaceIndex << '\n';
            }
        }
    }
}
void Molecule::display_all() const
{
    std::cout << "Index: " << index << '\n';
    std::cout << "Is empty: " << std::boolalpha << isEmpty << '\n';
    if (!isEmpty) {

        std::cout << "Parent complex index: " << myComIndex << '\n';
        std::cout << "Sub volume index: " << mySubVolIndex << '\n';
        std::cout << "Is a lipid: " << std::boolalpha << isLipid << '\n';
        std::cout << "Center of mass coordinate: " << comCoord << '\n';
        std::cout << "Interfaces:\n";
        for (const auto& iface : interfaceList) {
            std::cout << "\t---\n";
            std::cout << "\tRelative index: " << iface.relIndex << '\n';
            std::cout << "\tAbsolute index: " << iface.index << '\n';

            std::cout << "\tCoordinate: " << iface.coord << '\n';
            if (iface.stateIden != '\0')
                std::cout << "\tCurrent state: " << iface.stateIden << '\n';
            if (iface.isBound) {
                std::cout << "\tInteraction:\n";
                std::cout << "\t\tPartner index: " << iface.interaction.partnerIndex << '\n';
                std::cout << "\t\tPartner interface index " << iface.interaction.partnerIfaceIndex << '\n';
            }
        }
    }
}

void Molecule::display_assoc_icoords(const std::string& name)
{
    std::cout << name << ':' << std::endl;
    std::cout << std::setw(8) << std::setprecision(12) << std::right << tmpComCoord << std::endl;
    for (auto& icoord : tmpICoords)
        std::cout << std::setw(8) << std::setprecision(12) << std::right << icoord << std::endl;
    std::cout << std::endl;
}
void Molecule::display_my_coords(const std::string& name)
{
    std::cout << name << ':' << std::endl;
    std::cout << comCoord << '\n';
    for (const auto& iface : interfaceList) {
        std::cout << iface.coord << '\n';
    }
}

void Molecule::write_crd_file(std::ofstream& os) const
{
    os << tmpComCoord << std::endl;
    for (auto& icoord : tmpICoords)
        os << icoord << std::endl;
}
void Molecule::write_crd_file_cout() const
{
    std::cout << comCoord << std::endl;
    for (auto& iface : interfaceList)
        std::cout << iface.coord << std::endl;
}

void Molecule::set_tmp_association_coords()
{
    tmpComCoord = comCoord;
    for (auto& iface : interfaceList) {
        tmpICoords.emplace_back(iface.coord);
    }
}

void Molecule::clear_tmp_association_coords()
{
    tmpComCoord.zero_crds();
    tmpICoords.erase(tmpICoords.begin(), tmpICoords.end());
}

void Molecule::destroy(std::vector<int>& emptyMolList)
{
    /*! \ingroup Reactions
     * \brief Destroys the parent Molecule.
     *
     * Invoked by Complex::destroy(). Keeps the blank Molecule to fill upon the creation of a new Molecule.
     */

    // add to the list of empty Molecules
    Molecule::emptyMolList.push_back(index);

    // keep track of molecule types
    --MolTemplate::numEachMolType[molTypeIndex];

    myComIndex = -1;
    molTypeIndex = -1;
    mass = -1;
    trajStatus = TrajStatus::empty;

    // clear coordinates
    comCoord.zero_crds();
    interfaceList.clear();

    // clear association lists
    freelist.clear();
    bndpartner.clear();
    bndlist.clear();
    bndRxnList.clear();
    interfaceList.clear();

    // iterate the total number of molecules
    --numberOfMolecules;

    // set to void
    isEmpty = true;
}

void Molecule::create_random_coords(const MolTemplate& molTemplate, const Membrane& membraneObject)
{
    /*!
     * \brief Create random coordinates for a Molecule
     * rotation. Saves time, but could be changed easily.
     */

    comCoord.x = (membraneObject.waterBox.x * rand_gsl()) - (membraneObject.waterBox.x / 2.0);
    comCoord.y = (membraneObject.waterBox.y * rand_gsl()) - (membraneObject.waterBox.y / 2.0);

    bool outOfBox { false }; // TODO: commented out for testing purposes only
    // if the molecule is a lipid, place it along the bottom of the box and don't give it a rotation
    if (molTemplate.isLipid) {
        comCoord.z = -membraneObject.waterBox.z / 2.0;
        for (unsigned int ifaceItr { 0 }; ifaceItr < molTemplate.interfaceList.size(); ++ifaceItr) {
            interfaceList[ifaceItr].coord = Coord { comCoord + molTemplate.interfaceList[ifaceItr].iCoord };

            /* TODO: commented out for testing only */
            if (interfaceList[ifaceItr].coord.isOutOfBox(membraneObject)) {
                outOfBox = true;
                break;
            }

            // initialize the Interface::State to the default state (first listed)
            /*For each physical molecule, initialize the interfaces using
              the interfaces defined for the molTemplate, which were defined in parse_molFile.cpp
             at lines 1038.  Since each interface can exist in multiple states, the
             default choice here is to use the first state, usually an unbound state.
            */
            interfaceList[ifaceItr].index = molTemplate.interfaceList[ifaceItr].stateList[0].index;
            interfaceList[ifaceItr].relIndex = ifaceItr;
            interfaceList[ifaceItr].stateIden = molTemplate.interfaceList[ifaceItr].stateList[0].iden;
            interfaceList[ifaceItr].stateIndex = 0; // because by default we picked the first state
            interfaceList[ifaceItr].molTypeIndex = molTemplate.molTypeIndex;
        }
    } else {
        comCoord.z = (membraneObject.waterBox.z * rand_gsl()) - (membraneObject.waterBox.z / 2.0);

        // set interface coordinates, with a random rotation on the entire molecule
        // TODO: Commented this out for testing against old version
        Quat rotQuat { rand_gsl() * 2 - 1, rand_gsl() * 2 - 1, rand_gsl() * 2 - 1, rand_gsl() * 2 - 1 };
        rotQuat = rotQuat.unit();
        for (unsigned int ifaceItr { 0 }; ifaceItr < molTemplate.interfaceList.size(); ++ifaceItr) {
            Vector ifaceVec { Coord { comCoord + molTemplate.interfaceList[ifaceItr].iCoord } - comCoord };
            rotQuat.rotate(ifaceVec);
            interfaceList[ifaceItr].coord = Coord { comCoord + ifaceVec };

            // TODO: Commented out for testing
            if (interfaceList[ifaceItr].coord.isOutOfBox(membraneObject)) {
                outOfBox = true;
                break;
            }

            interfaceList[ifaceItr].index = molTemplate.interfaceList[ifaceItr].stateList[0].index;
            interfaceList[ifaceItr].relIndex = ifaceItr;
            interfaceList[ifaceItr].stateIden = molTemplate.interfaceList[ifaceItr].stateList[0].iden;
            interfaceList[ifaceItr].stateIndex = 0; // because by default we picked the first state
            interfaceList[ifaceItr].molTypeIndex = molTemplate.molTypeIndex;
        }
    }

    // TODO: Commented out for testing
    if (outOfBox)
        this->create_random_coords(molTemplate, membraneObject);
}

bool Molecule::operator==(const Molecule& rhs) const
{
    return std::tie(myComIndex, molTypeIndex, index, comCoord)
        == std::tie(rhs.myComIndex, rhs.molTypeIndex, rhs.index, rhs.comCoord);
}

bool Molecule::operator!=(const Molecule& rhs) const
{
    return std::tie(myComIndex, molTypeIndex, index, comCoord)
        != std::tie(rhs.myComIndex, rhs.molTypeIndex, rhs.index, rhs.comCoord);
}

void Molecule::Interaction::clear()
{
    partnerIfaceIndex = 0;
    partnerIndex = 0;
    conjBackRxn = 0;
}

void Molecule::update_association_coords(const Vector& vec)
{
    if (tmpICoords.empty()) {
        tmpComCoord = (vec + comCoord);
        for (auto& iCoord : interfaceList)
            tmpICoords.push_back(vec + iCoord.coord);
    } else {
        tmpComCoord = (vec + tmpComCoord);
        for (auto& iCoord : tmpICoords)
            iCoord = (vec + iCoord);
    }
}

// bool Molecule::overlapsWith(const Molecule& otherMol, const std::vector<ForwardRxn>& forwardRxns, const
// std::vector<BackRxn>& backRxns, const std::vector<MolTemplate>& molTemplateList)
//{
//    // check bounding sphere first
//    {
//        Vector sphereRad {comCoord - otherMol.comCoord};
//        sphereRad.calc_magnitude();
//        if (sphereRad.magnitude > molTemplateList[molTypeIndex].radius +
//        molTemplateList[otherMol.molTypeIndex].radius)
//            return false;
//    }
//
//}

/* COMPLEX */
void Complex::update_properties(
    const std::vector<Molecule>& moleculeList, const std::vector<MolTemplate>& molTemplateList)
{
    // update center of mass
    //update links to Surface
    linksToSurface = 0;
    double totMass { 0 };
    comCoord.zero_crds();
    for (auto& memMol : memberList) {
        totMass += moleculeList[memMol].mass;
        comCoord.x += moleculeList[memMol].comCoord.x * moleculeList[memMol].mass;
        comCoord.y += moleculeList[memMol].comCoord.y * moleculeList[memMol].mass;
        comCoord.z += moleculeList[memMol].comCoord.z * moleculeList[memMol].mass;
        linksToSurface += moleculeList[memMol].linksToSurface;
    }
    comCoord /= totMass;
    mass = totMass;

    // update radius (bounding sphere)
    /*For each protein, set radius to be non-zero, even for particles, so Diffusion calc
      is finite.
     */
    radius = 0;
    for (auto& memMol : memberList) {
        Vector distVec { moleculeList[memMol].comCoord - this->comCoord };
        distVec.calc_magnitude();
        if ((distVec.magnitude + molTemplateList[moleculeList[memMol].molTypeIndex].radius) > radius)
            radius = distVec.magnitude + molTemplateList[moleculeList[memMol].molTypeIndex].radius;
    }

    // update diffusion constants
    /*doesn't work for lipids: use Dt=c_avg/radius, where c_avg is a weighted average based on ci=radi*Di for each protein, weighted by 1/radi
      for Dr=c_avg/radius^3, where c_avg is a weighted average based on ci=radi^3*Dri for each protein, weighted by 1/radi^3
     
    */
    double currRad = radius;
    //std::cout <<" curr radius of complex: "<<index<<" radius: "<<radius<<std::endl;
    //    Coord sumD {};
    //Coord sumDr {};
    /*    double normDr=0;
    double normD=0;
    for (int memMol : memberList) {
        const MolTemplate& oneTemp { molTemplateList[moleculeList[memMol].molTypeIndex] };
        // Rotational diffusion constants
	//std::cout<<"curr c2 value: "<<pow(oneTemp.radius,  3.0) * oneTemp.Dr.x<<" curr c1 value: "<< oneTemp.D.x * oneTemp.radius<<'\n';
        sumDr.x += pow(oneTemp.radius,  6.0) * oneTemp.Dr.x;//6 is (r^3)*(r^3)
        sumDr.y += pow(oneTemp.radius, 6.0) * oneTemp.Dr.y;
        sumDr.z += pow(oneTemp.radius,  6.0) * oneTemp.Dr.z;

	normDr  +=pow(oneTemp.radius,  3.0);
	
        // Translational diffusion constants
        sumD.x += oneTemp.D.x * oneTemp.radius*oneTemp.radius;//r*r is r^2
        sumD.y += oneTemp.D.y * oneTemp.radius*oneTemp.radius;
        sumD.z += oneTemp.D.z * oneTemp.radius*oneTemp.radius;

	normD += oneTemp.radius;
	//std::cout<<" curr normDr: "<<normDr<<" currnormD: "<<normD<<std::endl;
    }
    
    double invR3=1.0/(currRad*currRad*currRad*normDr);
    Dr.x = sumDr.x*invR3;
    Dr.y = sumDr.x*invR3;
    Dr.z = sumDr.z*invR3;
    // upate translational diffusion constants
    double invDenom=1.0/(normD*currRad);
    D.x = sumD.x*invDenom;
    D.y = sumD.y*invDenom;
    D.z = sumD.z*invDenom;

    */
    // update rotational diffusion constants
    double inf = 1E300;
    Coord sumD {};
    Coord sumDr {};
    for (int memMol : memberList) {
        const MolTemplate& oneTemp { molTemplateList[moleculeList[memMol].molTypeIndex] };
        // Rotational diffusion constants
        sumDr.x += (oneTemp.Dr.x != 0) ? (1.0 / pow(oneTemp.Dr.x, 1.0 / 3.0)) : inf;
        sumDr.y += (oneTemp.Dr.y != 0) ? (1.0 / pow(oneTemp.Dr.y, 1.0 / 3.0)) : inf;
        sumDr.z += (oneTemp.Dr.z != 0) ? (1.0 / pow(oneTemp.Dr.z, 1.0 / 3.0)) : inf;

        // Translational diffusion constants
        sumD.x += (oneTemp.D.x != 0) ? (1.0 / oneTemp.D.x) : inf;
        sumD.y += (oneTemp.D.y != 0) ? (1.0 / oneTemp.D.y) : inf;
        sumD.z += (oneTemp.D.z != 0) ? (1.0 / oneTemp.D.z) : inf;
    }
    /*If linksToSurface>0, add in diffusion due to implicit lipids.*/
    for (int i = 0; i < linksToSurface; i++) {
        //	std::cout<<" complex:" <<index<<" has links to surface: "<<linksToSurface<<std::endl;
        //find the protein that is the Implicit Lipid.
        const MolTemplate& oneTemp { molTemplateList[moleculeList[iLipidIndex].molTypeIndex] };
        sumDr.x += (oneTemp.Dr.x != 0) ? (1.0 / pow(oneTemp.Dr.x, 1.0 / 3.0)) : inf;
        sumDr.y += (oneTemp.Dr.y != 0) ? (1.0 / pow(oneTemp.Dr.y, 1.0 / 3.0)) : inf;
        sumDr.z += (oneTemp.Dr.z != 0) ? (1.0 / pow(oneTemp.Dr.z, 1.0 / 3.0)) : inf;

        // Translational diffusion constants
        sumD.x += (oneTemp.D.x != 0) ? (1.0 / oneTemp.D.x) : inf;
        sumD.y += (oneTemp.D.y != 0) ? (1.0 / oneTemp.D.y) : inf;
        sumD.z += (oneTemp.D.z != 0) ? (1.0 / oneTemp.D.z) : inf;
    }
    Dr.x = 1.0 / (sumDr.x * sumDr.x * sumDr.x);
    Dr.y = 1.0 / (sumDr.y * sumDr.y * sumDr.y);
    Dr.z = 1.0 / (sumDr.z * sumDr.z * sumDr.z);

    if (Dr.x < 1E-50)
        Dr.x = 0;
    if (Dr.y < 1E-50)
        Dr.y = 0;
    if (Dr.z < 1E-50)
        Dr.z = 0;

    // upate translational diffusion constants
    D.x = 1.0 / sumD.x;
    D.y = 1.0 / sumD.y;
    D.z = 1.0 / sumD.z;

    //std::cout<<" Drx: "<<Dr.x<<" Dx: "<<D.x<<std::endl;
    if (D.z < 1E-50) {
        /*On the membrane. Only allow 2D. diffusion at certain intervals, to avoid generating too many 2D. Tables.*/
        double dtmp;
        if (D.x < 0.0001)
            dtmp = D.x * 100000;
        else if (D.x < 0.001)
            dtmp = D.x * 10000;
        else if (D.x < 0.01)
            dtmp = D.x * 1000;
        else if (D.x < 0.1)
            dtmp = D.x * 100;
        else
            dtmp = D.x * 100;

        /*Keep only one sig fig for <0.1, 2 for 0.1<d<10, 3 for 10<d<100, etc*/
        int d_ones = int(round(dtmp));
        //        std::cout << "FORCED. D._2D. to fewer sig-figs, starting value: " << D.x << " final value: ";
        /*Now put back in correct size*/
        if (D.x < 0.0001)
            D.x = d_ones * 0.00001;
        else if (D.x < 0.001)
            D.x = d_ones * 0.0001;
        else if (D.x < 0.01)
            D.x = d_ones * 0.001;
        else if (D.x < 0.1)
            D.x = d_ones * 0.01;
        else
            D.x = d_ones * 0.01;
        D.y = D.x; // set D.y equal to D.x.
    }

    if (D.x < 1E-50)
        D.x = 0;
    if (D.y < 1E-50)
        D.y = 0;
    if (D.z < 1E-50)
        D.z = 0;

    // update number of each constituent molecule
    numEachMol.clear();
    numEachMol.resize(molTemplateList.size());
    for (auto& memMol : memberList)
        ++numEachMol[moleculeList[memMol].molTypeIndex];
}

void Complex::display()
{
    std::cout << "Comcoords: " << comCoord << std::endl;
    std::cout << "Member list:";
    for (auto& mp : memberList)
        std::cout << ' ' << mp;
    std::cout << std::endl;
    std::cout << "D: " << D << '\n';
    std::cout << "Dr: " << Dr << '\n';
    std::cout << std::endl;
}

void Complex::display(const std::string& name)
{
    std::cout << name << std::endl;
    std::cout << "Comcoords: " << comCoord << std::endl;
    std::cout << "Member list:";
    for (auto& mp : memberList)
        std::cout << ' ' << mp;
    std::cout << std::endl;
    std::cout << "D: " << D << '\n';
    std::cout << "Dr: " << Dr << '\n';
    std::cout << std::endl;
}

void Complex::destroy(
    std::vector<Molecule>& moleculeList, std::vector<int>& emptyMolList, std::vector<Complex>& complexList, std::vector<int>& emptyComList)
{
    /*! \ingroup Reactions
     * \brief This function destroys its parent Complex.
     *
     * Invoked when a destruction reaction involving a member Molecule occurs. Destroys molecule but leaves shell behind
     * for later use if a new complex forms.
     *
     * \param [in] emptyMolList list of indices of empty Molecules in moleculeList
     * \param [in] emptyComList list of indices of empty Complexes in complexList
     */

    // add to the empty complex list
    Complex::emptyComList.push_back(index);

    // zero all the coordinates and diffusion coefficients
    comCoord.zero_crds();
    D.zero_crds();
    Dr.zero_crds();

    for (auto& memMol : memberList)
        moleculeList[memMol].destroy(emptyMolList);

    memberList.clear();
    numEachMol.clear();
    isEmpty = true;

    // iterate down the number of complexes in the system.
    trajStatus = TrajStatus::empty;
    --Complex::numberOfComplexes;
    /*
    // put the last non-empty complex in the list to the non-last empty slot
    int slotIndex { Complex::emptyComList.back() };
    int previousIndex { complexList.back().index };
    // if the empty one is the last, just pop it
    if (slotIndex == previousIndex) {
        complexList.pop_back();
        Complex::emptyComList.pop_back();
    } else {
        complexList[slotIndex] = complexList.back();
        complexList[slotIndex].index = slotIndex;
        complexList.pop_back();
        Complex::emptyComList.pop_back();

        // change the mol.myComIndex with previousIndex to slotIndex
        for (auto& mp : complexList[slotIndex].memberList) {
            moleculeList[mp].myComIndex = slotIndex;
        }
    }*/
}

void Complex::put_back_into_SimulVolume(
    int& itr, Molecule& errantMol, const Membrane& membraneObject, std::vector<Molecule>& moleculeList)
{
    std::cout << "Attempting to put complex " << index << " back into simulation volume...\n";
    display();

    Vector transVec { 0, 0, 0 };

    // check x dimension
    double xDiff = errantMol.comCoord.x - (membraneObject.waterBox.x / 2);
    if (xDiff > 0)
        transVec.x = -xDiff - 1E-6;
    else if (xDiff < -membraneObject.waterBox.x)
        transVec.x = -(xDiff + membraneObject.waterBox.x) + 1E-6;

    // check y dimension
    double yDiff = errantMol.comCoord.y - (membraneObject.waterBox.y / 2);
    if (yDiff > 0)
        transVec.y = -yDiff - 1E-6;
    else if (yDiff < -membraneObject.waterBox.y)
        transVec.y = -(yDiff + membraneObject.waterBox.y) + 1E-6;

    // check z dimension
    double zDiff = errantMol.comCoord.z - (membraneObject.waterBox.z / 2);
    if (zDiff > 0)
        transVec.z = -zDiff - 1E-6;
    else if (zDiff < -membraneObject.waterBox.z)
        transVec.z = -(zDiff + membraneObject.waterBox.z) + 1E-6;

    // TRANSLATE //
    comCoord += transVec;

    for (auto& memMol : memberList) {
        moleculeList.at(memMol).comCoord += transVec;
        for (auto& iface : moleculeList.at(memMol).interfaceList)
            iface.coord += transVec;
    }

    ++itr;
    if (itr == 1000) {
        std::cout << "Cannot fit complex " << index << " into simulation volume. Exiting...\n";
        exit(1);
    }
}

void Complex::propagate(std::vector<Molecule>& moleculeList)
{
    ++propCalled;
    // Create the quaternion
    double cosZ { cos(trajRot.z * 0.5) };
    double sinZ { sin(trajRot.z * 0.5) };
    double cosY { cos(trajRot.y * 0.5) };
    double sinY { sin(trajRot.y * 0.5) };
    double cosX { cos(trajRot.x * 0.5) };
    double sinX { sin(trajRot.x * 0.5) };

    Quat rotQuat {};
    rotQuat.x = (sinX * cosY * cosZ) - (cosX * sinY * sinZ);
    rotQuat.y = (cosX * sinY * cosZ) + (sinX * cosY * sinZ);
    rotQuat.z = (cosX * cosY * sinZ) - (sinX * sinY * cosZ);
    rotQuat.w = (cosX * cosY * cosZ) + (sinX * sinY * sinZ);
    rotQuat.unit();

    // update the member proteins
    for (auto mol : memberList) {
        //        if (moleculeList[mol].comCoord != this->comCoord) {
        Vector comVec { moleculeList[mol].comCoord - this->comCoord };
        rotQuat.rotate(comVec);
        moleculeList[mol].comCoord = Coord { comVec.x, comVec.y, comVec.z } + this->comCoord + trajTrans;
        //        } else {
        //            moleculeList[mol].comCoord += trajTrans;
        //        }

        // now rotate each member molecule of the complex
        for (auto& iface : moleculeList[mol].interfaceList) {
            // get the vector from the interface to the target interface
            Vector ifaceVec { iface.coord - comCoord };
            // rotate
            rotQuat.rotate(ifaceVec);
            iface.coord = Coord { ifaceVec.x, ifaceVec.y, ifaceVec.z } + comCoord + trajTrans;
        }
        moleculeList[mol].trajStatus = TrajStatus::propagated;
    }

    // propagate the complex's center of mass
    comCoord += trajTrans;

    // zero the propagation values
    trajTrans.zero_crds();
    trajRot.zero_crds();
}

void Complex::translate(Vector transVec, std::vector<Molecule>& moleculeList)
{
    for (auto memMol : memberList) {
        moleculeList[memMol].comCoord += transVec;
        for (auto& iface : moleculeList[memMol].interfaceList)
            iface.coord += transVec;
    }
}
