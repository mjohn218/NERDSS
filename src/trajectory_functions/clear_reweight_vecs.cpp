#include "trajectory_functions/trajectory_functions.hpp"

void clear_reweight_vecs(Molecule& oneMol)
{
    oneMol.prevlist = oneMol.currlist;
    oneMol.prevmyface = oneMol.currmyface;
    oneMol.prevpface = oneMol.currpface;
    oneMol.prevnorm = oneMol.currprevnorm;
    oneMol.ps_prev = oneMol.currps_prev;
    oneMol.prevsep = oneMol.currprevsep;
    oneMol.currlist.clear();
    oneMol.currmyface.clear();
    oneMol.currpface.clear();
    oneMol.currprevnorm.clear();
    oneMol.currps_prev.clear();
    oneMol.currprevsep.clear();
}
