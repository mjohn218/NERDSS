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

#include <cmath>
#include <iomanip>
#include <numeric>

#include "classes/class_Rxns.hpp"
#include "classes/class_bngl_parser.hpp"
#include "io/io.hpp"
#include "math/matrix.hpp"
#include "math/rand_gsl.hpp"
#include "parser/parser_functions.hpp"
#include "reactions/association/functions_for_spherical_system.hpp"
#include "split.cpp"

int Molecule::numberOfMolecules = 0;
int Molecule::maxID = 0;
int Complex::numberOfComplexes = 0;
int Complex::currNumberComTypes = 0;
int Complex::currNumberMolTypes = 0;
int Complex::maxID = 0;
std::vector<int> Complex::emptyComList{};
std::vector<int> Molecule::emptyMolList{};
std::vector<int> Complex::obs{};

int propCalled = 0;

bool skipLine(std::string line) {
  // check if the line is a comment or empty. I can't get regex to work with
  // this in particular, without also including needed lines that have comments
  // at the end
  return line.empty() || line[0] == '#';
}

Complex::Complex(const Molecule& mol, const MolTemplate& oneTemp)
    : Complex(mol.comCoord, oneTemp.D, oneTemp.Dr, Complex::maxID++) {
  mass = mol.mass;
  memberList.push_back(mol.index);
  index = mol.index;
  radius = oneTemp.radius;
  // Will elements of this array below be initialized to zero??
  numEachMol = std::vector<int>(MolTemplate::numMolTypes);
  ++numEachMol[oneTemp.molTypeIndex];

  lastNumberUpdateItrEachMol.resize(MolTemplate::numMolTypes);
  receivedFromNeighborRank = true;
}

Complex::Complex(int _index, const Molecule& _memMol,
                 const MolTemplate& _molTemp)
    : Complex(_memMol.comCoord, _molTemp.D, _molTemp.Dr, Complex::maxID++) {
  mass = _memMol.mass;
  memberList.push_back(_memMol.index);
  index = _index;
  radius = _molTemp.radius;
  // Will elements of this array below be initialized to zero??
  numEachMol = std::vector<int>(MolTemplate::numMolTypes);
  ++numEachMol[_molTemp.molTypeIndex];

  lastNumberUpdateItrEachMol.resize(MolTemplate::numMolTypes);
  receivedFromNeighborRank = true;
}

Complex::Complex(Coord comcoords, Coord D, Coord Dr, int idComplex)
    : comCoord(comcoords), D(D), Dr(Dr), id(idComplex) {
  receivedFromNeighborRank = true;
}

// Overloaded operators //
std::ostream& operator<<(std::ostream& os, const std::array<double, 3>& arr) {
  os << arr[0] << ' ' << arr[1] << ' ' << arr[2];
  return os;
}

void operator+(Coord& c, const double scal) {
  c.x += scal;
  c.y += scal;
  c.z += scal;
}

void operator*(Coord& c, const double scal) {
  c.x *= scal;
  c.y *= scal;
  c.z *= scal;
}

std::ostream& operator<<(std::ostream& os, const Molecule& mol) {
  os << mol.tmpComCoord << std::endl;
  for (auto& iface : mol.tmpICoords) os << iface << std::endl;
  return os;
}

// bool operator==(const IntCoordCont::Mol& mol1, const IntCoordCont::Mol& mol2)
//{
//    return (bool{ round(mol1.comCoord) == round(mol2.comCoord) } * bool{
//    mol1.proType == mol2.proType }
//        * bool{ mol1.reactantIface == mol2.reactantIface } * bool{
//        mol1.ifaceCoords == mol2.ifaceCoords });
//}

/* MOLECULE::IFACE */
void Molecule::Iface::change_state(int newRelIndex, int newAbsIndex,
                                   char newIden) {
  stateIndex = newRelIndex;
  index = newAbsIndex;
  stateIden = newIden;
}

/* MOLECULE */
void Molecule::display(const MolTemplate& molTemplate) const {
  std::cout << "Index: " << index << '\n';
  std::cout << "Id: " << id << '\n';
  std::cout << "Is empty: " << std::boolalpha << isEmpty << '\n';
  if (!isEmpty) {
    std::cout << "Type: " << molTemplate.molName << '\n';
    std::cout << "Complex index: " << myComIndex << '\n';
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
        std::cout << "\t\tPartner index: " << iface.interaction.partnerIndex
                  << '\n';
        std::cout << "\t\tPartner interface index "
                  << iface.interaction.partnerIfaceIndex << '\n';
      }
    }
  }
}
void Molecule::display_all() const {
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
        std::cout << "\t\tPartner index: " << iface.interaction.partnerIndex
                  << '\n';
        std::cout << "\t\tPartner interface index "
                  << iface.interaction.partnerIfaceIndex << '\n';
      }
    }
  }
}

void Molecule::display_assoc_icoords(const std::string& name) {
  std::cout << name << ':' << std::endl;
  std::cout << std::setw(8) << std::setprecision(12) << std::right
            << tmpComCoord << std::endl;
  for (auto& icoord : tmpICoords)
    std::cout << std::setw(8) << std::setprecision(12) << std::right << icoord
              << std::endl;
  std::cout << std::endl;
}
void Molecule::display_my_coords(const std::string& name) {
  std::cout << name << ':' << std::endl;
  std::cout << comCoord << '\n';
  for (const auto& iface : interfaceList) {
    std::cout << iface.coord << '\n';
  }
}

void Molecule::write_crd_file(std::ofstream& os) const {
  os << tmpComCoord << std::endl;
  for (auto& icoord : tmpICoords) os << icoord << std::endl;
}
void Molecule::write_crd_file_cout() const {
  std::cout << comCoord << std::endl;
  for (auto& iface : interfaceList) std::cout << iface.coord << std::endl;
}

void Molecule::set_tmp_association_coords() {
  tmpComCoord = comCoord;
  for (auto& iface : interfaceList) {
    tmpICoords.emplace_back(iface.coord);
  }
}

void Molecule::clear_tmp_association_coords() {
  tmpComCoord.zero_crds();
  tmpICoords.erase(tmpICoords.begin(), tmpICoords.end());
}

void Molecule::destroy() {
  /*! \ingroup Reactions
   * \brief Destroys the parent Molecule.
   *
   * Invoked by Complex::destroy(). Keeps the blank Molecule to fill upon the
   * creation of a new Molecule.
   */

  if (isEmpty) return;

  // add to the list of empty Molecules
  Molecule::emptyMolList.push_back(index);

  // keep track of molecule types
  --MolTemplate::numEachMolType[molTypeIndex];

  myComIndex = -1;
  molTypeIndex = -1;
  mass = -1;
  id = -1;
  trajStatus = TrajStatus::empty;

  // clear coordinates
  comCoord.zero_crds();
  interfaceList.clear();

  // clear association lists
  freelist.clear();
  bndpartner.clear();
  bndlist.clear();
  interfaceList.clear();

  // iterate the total number of molecules
  --numberOfMolecules;

  // set to void
  isEmpty = true;
}

// MPI_remove_from_one_rank destroys the molecule from one of two ranks
// that both kept information about,
// since it was in their shared-area:
void Molecule::MPI_remove_from_one_rank(std::vector<Molecule>& moleculeList,
                                        std::vector<Complex>& complexList) {
  /*! \ingroup Reactions
   * \brief Destroys the parent Molecule.
   *
   * Invoked by deserialize_molecules(). Keeps the blank Molecule to fill upon
   * the creation of a new Molecule.
   */

  if (isEmpty) return;

  // add to the list of empty Molecules
  Molecule::emptyMolList.push_back(index);

  // keep track of molecule types
  --MolTemplate::numEachMolType[molTypeIndex];

  myComIndex = -1;
  molTypeIndex = -1;
  mass = -1;
  id = -1;
  trajStatus = TrajStatus::empty;

  // clear coordinates
  comCoord.zero_crds();
  interfaceList.clear();

  // clear association lists
  freelist.clear();
  bndpartner.clear();
  bndlist.clear();
  interfaceList.clear();

  // iterate the total number of molecules
  --numberOfMolecules;

  // set to void
  isEmpty = true;
}

void Molecule::create_random_coords(const MolTemplate& molTemplate,
                                    const Membrane& membraneObject) {
  /*!
   * \brief Create random coordinates for a Molecule
   * rotation. Saves time, but could be changed easily.
   */
  if (membraneObject.isSphere) {
    double R = membraneObject.sphereR;

    comCoord.x =
        (membraneObject.sphereR * 2 * rand_gsl() - (membraneObject.sphereR));
    comCoord.y =
        (membraneObject.sphereR * 2 * rand_gsl() - (membraneObject.sphereR));
    comCoord.z =
        (membraneObject.sphereR * 2 * rand_gsl() - (membraneObject.sphereR));

    bool outOfBox{false};  // TODO: commented out for testing purposes only
    // if the molecule is a lipid, place it along the bottom of the box and
    // don't give it a rotation

    double molMag = sqrt(comCoord.x * comCoord.x + comCoord.y * comCoord.y +
                         comCoord.z * comCoord.z);

    if (molMag > R) {
      outOfBox = true;
    }

    if (molTemplate.isLipid) {
      double x = comCoord.x;
      double y = comCoord.y;
      double z = comCoord.z;

      comCoord.x = (R / sqrt(x * x + y * y + z * z)) * x;
      comCoord.y = (R / sqrt(x * x + y * y + z * z)) * y;
      comCoord.z = (R / sqrt(x * x + y * y + z * z)) * z;
      for (unsigned int ifaceItr{0};
           ifaceItr < molTemplate.interfaceList.size(); ++ifaceItr) {
        interfaceList[ifaceItr].coord =
            Coord{comCoord + molTemplate.interfaceList[ifaceItr].iCoord};
        Coord iCoord = molTemplate.interfaceList[ifaceItr].iCoord;
        double iLen = sqrt(iCoord.x * iCoord.x + iCoord.y * iCoord.y +
                           iCoord.z * iCoord.z);

        interfaceList[ifaceItr].coord.x = (R - iLen) * (comCoord.x) / R;
        interfaceList[ifaceItr].coord.y = (R - iLen) * (comCoord.y) / R;
        interfaceList[ifaceItr].coord.z = (R - iLen) * (comCoord.z) / R;
        /*
      //random position for check
      Coord a = find_spherical_coords(comCoord);
      double dangle = M_PI/20.0;
      interfaceList[ifaceItr].coord.x = a.z * sin(a.x) * cos(a.y + dangle);
      interfaceList[ifaceItr].coord.y = a.z * sin(a.x) * sin(a.y + dangle);
      interfaceList[ifaceItr].coord.z = a.z * cos(a.x);
      */
        // initialize the Interface::State to the default state (first listed)
        /*For each physical molecule, initialize the interfaces using
        the interfaces defined for the molTemplate, which were defined in
       parse_molFile.cpp at lines 1038.  Since each interface can exist in
       multiple states, the default choice here is to use the first state,
       usually an unbound state.
      */
        interfaceList[ifaceItr].index =
            molTemplate.interfaceList[ifaceItr].stateList[0].index;
        interfaceList[ifaceItr].relIndex = ifaceItr;
        interfaceList[ifaceItr].stateIden =
            molTemplate.interfaceList[ifaceItr].stateList[0].iden;
        interfaceList[ifaceItr].stateIndex =
            0;  // because by default we picked the first state
        interfaceList[ifaceItr].molTypeIndex = molTemplate.molTypeIndex;
      }
    } else {
      // set interface coordinates, with a random rotation on the entire
      // molecule
      // TODO: Commented this out for testing against old version
      Quat rotQuat{rand_gsl() * 2 - 1, rand_gsl() * 2 - 1, rand_gsl() * 2 - 1,
                   rand_gsl() * 2 - 1};
      rotQuat = rotQuat.unit();
      for (unsigned int ifaceItr{0};
           ifaceItr < molTemplate.interfaceList.size(); ++ifaceItr) {
        Vector ifaceVec{
            Coord{comCoord + molTemplate.interfaceList[ifaceItr].iCoord} -
            comCoord};
        rotQuat.rotate(ifaceVec);
        interfaceList[ifaceItr].coord = Coord{comCoord + ifaceVec};

        // TODO: Commented out for testing
        Coord iCoord = interfaceList[ifaceItr].coord;
        double iMag = sqrt(iCoord.x * iCoord.x + iCoord.y * iCoord.y +
                           iCoord.z * iCoord.z);
        if (iMag > R) {
          outOfBox = true;
          break;
        }

        interfaceList[ifaceItr].index =
            molTemplate.interfaceList[ifaceItr].stateList[0].index;
        interfaceList[ifaceItr].relIndex = ifaceItr;
        interfaceList[ifaceItr].stateIden =
            molTemplate.interfaceList[ifaceItr].stateList[0].iden;
        interfaceList[ifaceItr].stateIndex = 0;
        interfaceList[ifaceItr].molTypeIndex = molTemplate.molTypeIndex;
      }
    }

    // TODO: Commented out for testing
    if (outOfBox) this->create_random_coords(molTemplate, membraneObject);

  } else {
    // TODO: nenadko: randomize so that it gets in the domain of this rank
    comCoord.x =
        (membraneObject.waterBox.xLeft +
         (membraneObject.waterBox.xRight - membraneObject.waterBox.xLeft) *
             rand_gsl());
    comCoord.y = (membraneObject.waterBox.y * rand_gsl()) -
                 (membraneObject.waterBox.y / 2.0);

    bool outOfBox{false};  // TODO: commented out for testing purposes only
    // if the molecule is a lipid, place it along the bottom of the box and
    // don't give it a rotation
    if (molTemplate.isLipid) {
      comCoord.z = -membraneObject.waterBox.z / 2.0;
      for (unsigned int ifaceItr{0};
           ifaceItr < molTemplate.interfaceList.size(); ++ifaceItr) {
        interfaceList[ifaceItr].coord =
            Coord{comCoord + molTemplate.interfaceList[ifaceItr].iCoord};

        /* TODO: commented out for testing only */
        if (interfaceList[ifaceItr].coord.isOutOfBox(membraneObject)) {
          outOfBox = true;
          break;
        }

        // initialize the Interface::State to the default state (first listed)
        /*For each physical molecule, initialize the interfaces using
      the interfaces defined for the molTemplate, which were defined in
     parse_molFile.cpp at lines 1038.  Since each interface can exist in
     multiple states, the default choice here is to use the first state, usually
     an unbound state.
    */
        interfaceList[ifaceItr].index =
            molTemplate.interfaceList[ifaceItr].stateList[0].index;
        interfaceList[ifaceItr].relIndex = ifaceItr;
        interfaceList[ifaceItr].stateIden =
            molTemplate.interfaceList[ifaceItr].stateList[0].iden;
        interfaceList[ifaceItr].stateIndex =
            0;  // because by default we picked the first state
        interfaceList[ifaceItr].molTypeIndex = molTemplate.molTypeIndex;
      }
    } else {
      comCoord.z = (membraneObject.waterBox.z * rand_gsl()) -
                   (membraneObject.waterBox.z / 2.0);

      // set interface coordinates, with a random rotation on the entire
      // molecule
      // TODO: Commented this out for testing against old version
      Quat rotQuat{rand_gsl() * 2 - 1, rand_gsl() * 2 - 1, rand_gsl() * 2 - 1,
                   rand_gsl() * 2 - 1};
      rotQuat = rotQuat.unit();
      for (unsigned int ifaceItr{0};
           ifaceItr < molTemplate.interfaceList.size(); ++ifaceItr) {
        Vector ifaceVec{
            Coord{comCoord + molTemplate.interfaceList[ifaceItr].iCoord} -
            comCoord};
        rotQuat.rotate(ifaceVec);
        interfaceList[ifaceItr].coord = Coord{comCoord + ifaceVec};

        // TODO: Commented out for testing
        if (interfaceList[ifaceItr].coord.isOutOfBox(membraneObject)) {
          outOfBox = true;
          break;
        }

        interfaceList[ifaceItr].index =
            molTemplate.interfaceList[ifaceItr].stateList[0].index;
        interfaceList[ifaceItr].relIndex = ifaceItr;
        interfaceList[ifaceItr].stateIden =
            molTemplate.interfaceList[ifaceItr].stateList[0].iden;
        interfaceList[ifaceItr].stateIndex =
            0;  // because by default we picked the first state
        interfaceList[ifaceItr].molTypeIndex = molTemplate.molTypeIndex;
      }
    }

    // TODO: Commented out for testing
    if (outOfBox) this->create_random_coords(molTemplate, membraneObject);
  }
}

bool Molecule::operator==(const Molecule& rhs) const {
  return std::tie(myComIndex, molTypeIndex, index, comCoord) ==
         std::tie(rhs.myComIndex, rhs.molTypeIndex, rhs.index, rhs.comCoord);
}

bool Molecule::operator!=(const Molecule& rhs) const {
  return std::tie(myComIndex, molTypeIndex, index, comCoord) !=
         std::tie(rhs.myComIndex, rhs.molTypeIndex, rhs.index, rhs.comCoord);
}

void Molecule::Interaction::clear() {
  partnerIfaceIndex = -1;
  partnerIndex = -1;
  partnerId = -1;
  conjBackRxn = -1;
}

void Molecule::update_association_coords(const Vector& vec) {
  if (tmpICoords.empty()) {
    tmpComCoord = (vec + comCoord);
    for (auto& iCoord : interfaceList) tmpICoords.push_back(vec + iCoord.coord);
  } else {
    tmpComCoord = (vec + tmpComCoord);
    for (auto& iCoord : tmpICoords) iCoord = (vec + iCoord);
  }
}

void Molecule::print(MpiContext &mpiContext) const {
  const char *TrajStatusTypes[] = {"none",       "dissociate",     "associated",
                                   "propagated", "canBeResampled", "empty"};

  std::cout << "Molecule properties:\n";
  std::cout << "Index: " << index << "; ";
  std::cout << "Id: " << id << "; ";
  std::cout << "Parent Complex Index: " << myComIndex << "; ";
  std::cout << "Parent Complex Id: " << complexId << "; ";
  std::cout << "Molecule Type Index: " << molTypeIndex << "; ";
  std::cout << "mySubVolIndex: " << mySubVolIndex << '\n';
  std::cout << "isEmpty: " << isEmpty << '\n';
  std::cout << "receivedFromNeighborRank: " << receivedFromNeighborRank << '\n';
  std::cout << "trajStatus: " << TrajStatusTypes[static_cast<std::underlying_type<TrajStatus>::type>(
                trajStatus)] << '\n';
  std::cout << "isAssociated: " << isAssociated << "; ";
  std::cout << "isDissociated: " << isDissociated << '\n';

  std::cout << "isLeftGhost: " << isLeftGhost << "; ";
  std::cout << "isLeftEdge: " << isLeftEdge << "; ";
  std::cout << "isRightEdge: " << isRightEdge << "; ";
  std::cout << "isRightGhost: " << isRightGhost << '\n';

  std::cout << "(x,y,z): " << comCoord.x << ", " << comCoord.y << ", " << comCoord.z << '\n';

  int xBin = int((comCoord.x + (*(mpiContext.membraneObject)).waterBox.x / 2) /
            (*(mpiContext.simulVolume)).subCellSize.x) - mpiContext.xOffset;

  std::cout << "xBin: " << xBin << '\n';
  mpiContext.print_xBins();

  std::cout << "freelist: ";
  for (auto& it : freelist) std::cout << it << ' ';
  std::cout << '\n';

  std::cout << "bndlist: ";
  for (auto& it : bndlist) std::cout << it << ' ';
  std::cout << '\n';

  std::cout << "bndpartner: ";
  for (auto& it : bndpartner) std::cout << it << ' ';
  std::cout << '\n';
}

// bool Molecule::overlapsWith(const Molecule& otherMol, const
// std::vector<ForwardRxn>& forwardRxns, const std::vector<BackRxn>& backRxns,
// const std::vector<MolTemplate>& molTemplateList)
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
    const std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList) {
  // update center of mass
  // update links to Surface
  linksToSurface = 0;
  double totMass{0};
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
  /*For each protein, set radius to be non-zero, even for particles, so
    Diffusion calc is finite.
   */
  radius = 0;
  for (auto& memMol : memberList) {
    Vector distVec{moleculeList[memMol].comCoord - this->comCoord};
    distVec.calc_magnitude();
    if ((distVec.magnitude +
         molTemplateList[moleculeList[memMol].molTypeIndex].radius) > radius)
      radius = distVec.magnitude +
               molTemplateList[moleculeList[memMol].molTypeIndex].radius;
  }

  // update diffusion constants
  /*doesn't work for lipids: use Dt=c_avg/radius, where c_avg is a weighted
    average based on ci=radi*Di for each protein, weighted by 1/radi for
    Dr=c_avg/radius^3, where c_avg is a weighted average based on ci=radi^3*Dri
    for each protein, weighted by 1/radi^3

  */
  double currRad = radius;
  // std::cout <<" curr radius of complex: "<<index<<" radius:
  // "<<radius<<std::endl;
  //     Coord sumD {};
  // Coord sumDr {};
  /*    double normDr=0;
  double normD=0;
  for (int memMol : memberList) {
      const MolTemplate& oneTemp {
  molTemplateList[moleculeList[memMol].molTypeIndex] };
      // Rotational diffusion constants
      //std::cout<<"curr c2 value: "<<pow(oneTemp.radius,  3.0) *
  oneTemp.Dr.x<<" curr c1 value: "<< oneTemp.D.x * oneTemp.radius<<'\n'; sumDr.x
  += pow(oneTemp.radius,  6.0) * oneTemp.Dr.x;//6 is (r^3)*(r^3) sumDr.y +=
  pow(oneTemp.radius, 6.0) * oneTemp.Dr.y; sumDr.z += pow(oneTemp.radius,  6.0)
  * oneTemp.Dr.z;

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
  Coord sumD{};
  Coord sumDr{};
  for (int memMol : memberList) {
    const MolTemplate& oneTemp{
        molTemplateList[moleculeList[memMol].molTypeIndex]};
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
    //	std::cout<<" complex:" <<index<<" has links to surface:
    //"<<linksToSurface<<std::endl; find the protein that is the Implicit Lipid.
    const MolTemplate& oneTemp{
        molTemplateList[moleculeList[iLipidIndex].molTypeIndex]};
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

  if (Dr.x < 1E-50) Dr.x = 0;
  if (Dr.y < 1E-50) Dr.y = 0;
  if (Dr.z < 1E-50) Dr.z = 0;

  // upate translational diffusion constants
  D.x = 1.0 / sumD.x;
  D.y = 1.0 / sumD.y;
  D.z = 1.0 / sumD.z;

  // std::cout<<" Drx: "<<Dr.x<<" Dx: "<<D.x<<std::endl;
  if (D.z < 1E-50) {
    /*On the membrane. Only allow 2D. diffusion at certain intervals, to avoid
     * generating too many 2D. Tables.*/
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
    //        std::cout << "FORCED. D._2D. to fewer sig-figs, starting value: "
    //        << D.x << " final value: ";
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
    D.y = D.x;  // set D.y equal to D.x.
  }

  if (D.x < 1E-50) D.x = 0;
  if (D.y < 1E-50) D.y = 0;
  if (D.z < 1E-50) D.z = 0;

  // update number of each constituent molecule
  numEachMol.clear();
  numEachMol.resize(molTemplateList.size());
  for (auto& memMol : memberList)
    ++numEachMol[moleculeList[memMol].molTypeIndex];
}

void Complex::display() {
  std::cout << "Comcoords: " << comCoord << std::endl;
  std::cout << "Member list:";
  for (auto& mp : memberList) std::cout << ' ' << mp;
  std::cout << std::endl;
  std::cout << "D: " << D << '\n';
  std::cout << "Dr: " << Dr << '\n';
  std::cout << std::endl;
}

void Complex::display(const std::string& name) {
  std::cout << name << std::endl;
  std::cout << "Comcoords: " << comCoord << std::endl;
  std::cout << "Member list:";
  for (auto& mp : memberList) std::cout << ' ' << mp;
  std::cout << std::endl;
  std::cout << "D: " << D << '\n';
  std::cout << "Dr: " << Dr << '\n';
  std::cout << std::endl;
}

void Complex::print(MpiContext &mpiContext) const {
  const char *TrajStatusTypes[] = {"none",       "dissociate",     "associated",
                                   "propagated", "canBeResampled", "empty"};

  std::cout << "Complex properties:\n";
  std::cout << "Index: " << index << "; ";
  std::cout << "Id: " << id << "; ";
  std::cout << "isEmpty: " << isEmpty << '\n';
  std::cout << "receivedFromNeighborRank: " << receivedFromNeighborRank << '\n';
  std::cout << "trajStatus: " << TrajStatusTypes[static_cast<std::underlying_type<TrajStatus>::type>(
                trajStatus)] << '\n';
  std::cout << "isLeftGhost: " << isLeftGhost << "; ";
  std::cout << "isLeftEdge: " << isLeftEdge << "; ";
  std::cout << "isRightEdge: " << isRightEdge << "; ";
  std::cout << "isRightGhost: " << isRightGhost << '\n';
  std::cout << "(x,y,z): " << comCoord.x << ", " << comCoord.y << ", " << comCoord.z << '\n';

  std::cout << "memberList: ";
  for (auto& it : memberList) std::cout << it << ' ';
  std::cout << '\n';
}

void Complex::destroy(std::vector<Molecule>& moleculeList,
                      std::vector<Complex>& complexList) {
  /*! \ingroup Reactions
   * \brief This function destroys its parent Complex.
   *
   * Invoked when a destruction reaction involving a member Molecule occurs.
   * Destroys molecule but leaves shell behind for later use if a new complex
   * forms.
   *
   */

  if (isEmpty) return;
  // add to the empty complex list
  Complex::emptyComList.push_back(index);

  // zero all the coordinates and diffusion coefficients
  comCoord.zero_crds();
  D.zero_crds();
  Dr.zero_crds();

  for (auto& memMol : memberList) moleculeList[memMol].destroy();

  memberList.clear();
  numEachMol.clear();
  lastNumberUpdateItrEachMol.clear();
  isEmpty = true;
  id = -1;

  // iterate down the number of complexes in the system.
  trajStatus = TrajStatus::empty;
  --Complex::numberOfComplexes;
}

void Complex::put_back_into_SimulVolume(
    int& itr, Molecule& errantMol, const Membrane& membraneObject,
    std::vector<Molecule>& moleculeList,
    const std::vector<MolTemplate>& molTemplateList) {
  std::cout << "Attempting to put complex " << index
            << " back into simulation volume...\n";
  display();

  Vector transVec{0, 0, 0};

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
  // comCoord += transVec;

  for (auto& memMol : memberList) {
    moleculeList.at(memMol).comCoord += transVec;
    for (auto& iface : moleculeList.at(memMol).interfaceList)
      iface.coord += transVec;
  }
  this->update_properties(moleculeList, molTemplateList);

  ++itr;
  if (itr == 1000) {
    std::cout << "Cannot fit complex " << index
              << " into simulation volume. Exiting...\n";
    exit(1);
  }
}

// void Complex::propagate(std::vector<Molecule>& moleculeList)
// {
//     ++propCalled;
//     // Create the quaternion
//     double cosZ { cos(trajRot.z * 0.5) };
//     double sinZ { sin(trajRot.z * 0.5) };
//     double cosY { cos(trajRot.y * 0.5) };
//     double sinY { sin(trajRot.y * 0.5) };
//     double cosX { cos(trajRot.x * 0.5) };
//     double sinX { sin(trajRot.x * 0.5) };

//     Quat rotQuat {};
//     rotQuat.x = (sinX * cosY * cosZ) - (cosX * sinY * sinZ);
//     rotQuat.y = (cosX * sinY * cosZ) + (sinX * cosY * sinZ);
//     rotQuat.z = (cosX * cosY * sinZ) - (sinX * sinY * cosZ);
//     rotQuat.w = (cosX * cosY * cosZ) + (sinX * sinY * sinZ);
//     rotQuat.unit();

//     // update the member proteins
//     for (auto mol : memberList) {
//         //        if (moleculeList[mol].comCoord != this->comCoord) {
//         Vector comVec { moleculeList[mol].comCoord - this->comCoord };
//         rotQuat.rotate(comVec);
//         moleculeList[mol].comCoord = Coord { comVec.x, comVec.y, comVec.z } +
//         this->comCoord + trajTrans;
//         //        } else {
//         //            moleculeList[mol].comCoord += trajTrans;
//         //        }

//         // now rotate each member molecule of the complex
//         for (auto& iface : moleculeList[mol].interfaceList) {
//             // get the vector from the interface to the target interface
//             Vector ifaceVec { iface.coord - comCoord };
//             // rotate
//             rotQuat.rotate(ifaceVec);
//             iface.coord = Coord { ifaceVec.x, ifaceVec.y, ifaceVec.z } +
//             comCoord + trajTrans;
//         }
//         moleculeList[mol].trajStatus = TrajStatus::propagated;
//     }

//     // propagate the complex's center of mass
//     comCoord += trajTrans;

//     // zero the propagation values
//     trajTrans.zero_crds();
//     trajRot.zero_crds();
// }

void Complex::propagate(std::vector<Molecule>& moleculeList,
                        const Membrane membraneObject,
                        const std::vector<MolTemplate>& molTemplateList) {
  ++propCalled;

  /*debug*/
  // std::cout << "propCalled: " << propCalled << std::endl;
  // std::cout << "startComCoord: " << std::fixed << std::setprecision(20) <<
  // comCoord.x << " " << comCoord.y << " " << comCoord.z << std::endl;
  //  std::cout << "trajTrans: " << std::setprecision(20) << trajTrans.x << " "
  //  << trajTrans.y << " " << trajTrans.z << std::endl; std::cout << "trajRot:
  //  " << std::setprecision(20) << trajRot.x << " " << trajRot.y << " " <<
  //  trajRot.z << " " << std::endl;
  /*end debug*/

  // for the complex on the sphere surface, propagation is special
  if (membraneObject.isSphere && this->D.z < 1E-14) {
    Coord COM = comCoord;
    Coord COMnew = comCoord + trajTrans;
    std::array<double, 9> Crdset = inner_coord_set(COM, COMnew);
    std::array<double, 9> Crdsetnew = inner_coord_set_new(COM, COMnew);
    // propagate the com of the complex
    // comCoord += trajTrans;
    // get the rotation angle: dangle
    double Rotangle = trajRot.x;  // we select the rotation angle x as the angle
                                  // that the complex rotate on the sphere
    // update the member proteins
    for (auto mol : memberList) {
      Coord targmol = moleculeList[mol].comCoord;
      targmol = translate_on_sphere(targmol, COM, COMnew, Crdset, Crdsetnew);
      moleculeList[mol].comCoord =
          rotate_on_sphere(targmol, COMnew, Crdsetnew, Rotangle);
      moleculeList[mol].trajStatus = TrajStatus::propagated;
      // now update each interface of the molecule
      for (auto& iface : moleculeList[mol].interfaceList) {
        Coord targiface = iface.coord;
        targiface =
            translate_on_sphere(targiface, COM, COMnew, Crdset, Crdsetnew);
        iface.coord = rotate_on_sphere(targiface, COMnew, Crdsetnew, Rotangle);
      }
    }
    this->update_properties(moleculeList, molTemplateList);
    trajStatus = TrajStatus::propagated;
  } else {  // for the complex in solution or on box surface
    // Create the quaternion
    double cosZ{cos(trajRot.z * 0.5)};
    double sinZ{sin(trajRot.z * 0.5)};
    double cosY{cos(trajRot.y * 0.5)};
    double sinY{sin(trajRot.y * 0.5)};
    double cosX{cos(trajRot.x * 0.5)};
    double sinX{sin(trajRot.x * 0.5)};

    Quat rotQuat{};
    rotQuat.x = (sinX * cosY * cosZ) - (cosX * sinY * sinZ);
    rotQuat.y = (cosX * sinY * cosZ) + (sinX * cosY * sinZ);
    rotQuat.z = (cosX * cosY * sinZ) - (sinX * sinY * cosZ);
    rotQuat.w = (cosX * cosY * cosZ) + (sinX * sinY * sinZ);
    rotQuat.unit();

    // update the member proteins
    for (auto mol : memberList) {
      // TODO: nenadko
      // take xbin from
      // moleculeList[mol].comCoord

      Vector comVec{moleculeList[mol].comCoord - this->comCoord};
      rotQuat.rotate(comVec);
      moleculeList[mol].comCoord =
          Coord{comVec.x, comVec.y, comVec.z} + this->comCoord + trajTrans;

      // take xbin from
      // moleculeList[mol].comCoord
      // if(moved from non ghosted to ghosted)
      //   enter_ghosted_zone();
      //  now rotate each member molecule of the complex
      for (auto& iface : moleculeList[mol].interfaceList) {
        // get the vector from the interface to the target interface
        Vector ifaceVec{iface.coord - comCoord};
        // rotate
        rotQuat.rotate(ifaceVec);
        // std::cout << "=======" << moleculeList[mol].index << " " <<
        // iface.coord.x << " z:" << iface.coord.z << std::endl;
        iface.coord =
            Coord{ifaceVec.x, ifaceVec.y, ifaceVec.z} + comCoord + trajTrans;
        // std::cout << "=======" << moleculeList[mol].index << " " <<
        // iface.coord.x << " z:" << iface.coord.z << std::endl;
      }
      moleculeList[mol].trajStatus = TrajStatus::propagated;
    }

    // propagate the complex's center of mass
    // std::cout << "comCoord: " << std::fixed << std::setprecision(20) <<
    // comCoord.x << " " << comCoord.y << " " << comCoord.z << std::endl;
    // comCoord += trajTrans;
    this->update_properties(moleculeList, molTemplateList);
    trajStatus = TrajStatus::propagated;
    // std::cout << "comCoord: " << std::setprecision(20) << comCoord.x << " "
    // << comCoord.y << " " << comCoord.z << std::endl;
  }
  // zero the propagation values
  trajTrans.zero_crds();
  trajRot.zero_crds();
}

// only used for the temporary movement on sphere
//  all the input coords are cardesian coords
void Complex::update_association_coords_sphere(
    std::vector<Molecule>& moleculeList, Coord iface, Coord ifacenew) {
  Coord COM = iface;
  Coord COMnew = ifacenew;
  std::array<double, 9> Crdset = inner_coord_set(COM, COMnew);
  std::array<double, 9> Crdsetnew = inner_coord_set_new(COM, COMnew);

  // update the member proteins
  for (auto mol : memberList) {
    Coord targ = moleculeList[mol].comCoord;
    moleculeList[mol].tmpComCoord =
        translate_on_sphere(targ, COM, COMnew, Crdset, Crdsetnew);
    // now update each interface of the molecule
    for (int i = 0; i < moleculeList[mol].interfaceList.size(); i++) {
      auto ifacetmp = moleculeList[mol].interfaceList[i];
      targ = ifacetmp.coord;
      Coord ifacecrds =
          translate_on_sphere(targ, COM, COMnew, Crdset, Crdsetnew);
      if (moleculeList[mol].tmpICoords.empty()) {
        moleculeList[mol].tmpICoords.push_back(ifacecrds);
      } else {
        moleculeList[mol].tmpICoords[i] = ifacecrds;
      }
    }
  }
}

void Complex::translate(Vector transVec, std::vector<Molecule>& moleculeList) {
  for (auto memMol : memberList) {
    moleculeList[memMol].comCoord += transVec;
    for (auto& iface : moleculeList[memMol].interfaceList)
      iface.coord += transVec;
  }
}

void Molecule::create_position_implicit_lipid(Molecule& reactMol1,
                                              int ifaceIndex2,
                                              double bindRadius,
                                              const Membrane& membraneObject) {
  if (membraneObject.isSphere) {
    Coord displace = interfaceList[ifaceIndex2].coord - comCoord;
    double mag = displace.get_magnitude();
    /*
double lambda = (membraneObject.sphereR -
mag)/(reactMol1.comCoord.get_magnitude()); interfaceList[ifaceIndex2].coord =
lambda*reactMol1.comCoord; double lambda2 =
(membraneObject.sphereR)/(reactMol1.comCoord.get_magnitude()); comCoord =
lambda2*reactMol1.comCoord;
*/
    // the implicit lipid is set out of sphere, then the protein will bind onto
    // the sphere. and the boundary condition is easy to carry out.
    double lambda1 = (membraneObject.sphereR + mag + bindRadius) /
                     (reactMol1.comCoord.get_magnitude());
    comCoord = lambda1 * reactMol1.comCoord;
    double lambda2 = (membraneObject.sphereR + bindRadius) /
                     (reactMol1.comCoord.get_magnitude());
    interfaceList[ifaceIndex2].coord = lambda2 * reactMol1.comCoord;
  } else {
    Coord displace = interfaceList[ifaceIndex2].coord - comCoord;
    double shift =
        0.1;  // do not put right underneath, so that sigma starts at non-zero.
    comCoord.x = reactMol1.comCoord.x + shift;
    comCoord.y = reactMol1.comCoord.y - shift;
    comCoord.z = -membraneObject.waterBox.z / 2.0;  // + membraneObject.RS3D;
    Coord direction{0.0, 0.0, 1.0};
    interfaceList[ifaceIndex2].coord =
        comCoord + displace.get_magnitude() * direction;
  }
}
