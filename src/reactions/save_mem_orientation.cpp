#include "reactions/association/association.hpp"

/*Calculate rotation matrix for orienting one molecule to itself (at another timepoint, e.g.)
  One of these molecules (the first one) is true coords: baseTarget, the other is temp coords (base1).
 */
Quat save_mem_orientation(Molecule baseTarget, Molecule baseTmp, MolTemplate onePro)
{
    /*Reorient to template to ensure rigid structure is maintained*/
    /*Need to subtract off COM.*/
    // subtract_off_com_position(base1);
    //subtract_off_com_position(baseTarget);
    //  write_base(base1);
    //write_base(baseTarget);
    MolTemplate copy; //it does need to know whether it is a rod or NOT.
    copy.isRod = onePro.isRod;
    /*Copy the coordinates of base1 into another molecule*/

    //  cout <<"Copy Full mole into Protein class "<<endl;
    //cout <<" ninterface: "<<base1.ninterface <<endl;
    //copy.ninterface=base1.ninterface;
    std::string iface = "tmp";
    //std::cout << " IN SAVE MEM ORIENTATION TARGET COORDS: " << std::endl;
    //baseTarget.write_crd_file_cout();//iface coords
    copy.comCoord = Coord { 0.0, 0.0, 0.0 };
    for (int i = 0; i < baseTarget.interfaceList.size(); i++) {
        copy.interfaceList.emplace_back(iface, Coord { baseTarget.interfaceList[i].coord - baseTarget.comCoord }); //subtract off com if non-zero
    }

    for (int i = 0; i < baseTmp.interfaceList.size(); i++) {
        baseTmp.tmpICoords[i] = baseTmp.tmpICoords[i] - baseTmp.tmpComCoord; //subtract off com
    }
    baseTmp.tmpComCoord = Coord { 0.0, 0.0, 0.0 };
    // std::cout <<" IN SAVE MEM ORIENTATION MOVABLE COORDS: "<<std::endl;
    //baseTmp.display_assoc_icoords("temp coords");//tmp Iface coords

    //std::cout <<" Copied base, now find out rotation matrix "<<std::endl;
    //copy.display("Copy_of_Target");
    //Second input here will use the tmpCoords of the protein!

    Quat qRot = orient_crds_to_template(copy, baseTmp); //first input is a MolTemplate, second is a Molecule
    //std::cout<<" QUAT: "<<qRot.x<<' '<<qRot.y<<' '<<qRot.z<<' '<<qRot.w<<std::endl;
    //cout <<"Reoriented to template Proten. "<<endl;
    //copy.display("Copy_after ROTATION");
    //std::cout <<" IN SAVE MEM ORIENTATION COORDS AFTER ROTATION: "<<std::endl;
    // baseTmp.display_assoc_icoords("rotated coords");//tmp Iface coords

    return qRot;
}
