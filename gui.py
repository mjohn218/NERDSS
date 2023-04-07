import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import numpy as np

def gui():
    class Mol:
        def __init__(self, name, sites):
            self.name = name
            self.sites = sites
    
    class Site:
        def __init__(self, name, states = []):
            self.name = name
            self.states = states

    class State:
        def __init__(self, name):
            self.name = name
            
    # Helper function to create labeled entry fields
    def create_labeled_entry_with_hint(parent, label_text, hint_text, row, column, width=50):
        label = tk.Label(parent, text=label_text)
        entry = tk.Entry(parent, width=width)
        label.grid(row=row, column=column, sticky='e')
        entry.grid(row=row, column=column + 1)

        entry.insert(0, hint_text)
        entry.bind("<FocusIn>", lambda event, entry=entry, hint_text=hint_text: on_focus_in(event, entry, hint_text))
        entry.bind("<FocusOut>", lambda event, entry=entry, hint_text=hint_text: on_focus_out(event, entry, hint_text))

        return entry

    def on_focus_in(event, entry, hint_text):
        if entry.get() == hint_text:
            entry.delete(0, tk.END)

    def on_focus_out(event, entry, hint_text):
        if not entry.get():
            entry.insert(0, hint_text)

    def create_labeled_combobox_entry(parent, label_text, options_list, row, column, width=30):
        label = tk.Label(parent, text=label_text)
        combobox = ttk.Combobox(parent, values=options_list, width=width)
        label.grid(row=row, column=column, sticky='e')
        combobox.grid(row=row, column=column + 1)

        return combobox
    
    def create_labeled_combobox_entry_with_callback(parent, label_text, options, row, column, callback, width=30):
        label = tk.Label(parent, text=label_text)
        combobox = ttk.Combobox(parent, values=options, width=width)
        combobox.bind("<<ComboboxSelected>>", callback)
        label.grid(row=row, column=column, sticky="e")
        combobox.grid(row=row, column=column + 1)
        return combobox

    def generate_inp_file():
        with open("parms.inp", "w") as f:
            f.write("start parameters\n")
            f.write(f"    nItr = {nItr_entry.get()}\n")
            f.write(f"    timeStep = {timeStep_entry.get()}\n")
            f.write(f"    timeWrite = {timeWrite_entry.get()}\n")
            f.write(f"    trajWrite = {trajWrite_entry.get()}\n")
            f.write(f"    pdbWrite = {pdbWrite_entry.get()}\n")
            f.write(f"    transitionWrite = {transitionWrite_entry.get()}\n")
            f.write(f"    restartWrite = {restartWrite_entry.get()}\n")
            f.write(f"    checkPoint = {checkPoint_entry.get()}\n")
            f.write(f"    scaleMaxDisplace = {scaleMaxDisplace_entry.get()}\n")
            f.write(f"    overlapSepLimit = {overlapSepLimit_entry.get()}\n")
            f.write(f"    clusterOverlapCheck = {clusterOverlapCheck_entry.get()}\n")
            f.write("end parameters\n\n")

            f.write("start boundaries\n")
            isSphere_value = isSphere_entry.get().strip().lower() == "true"
            if not isSphere_value:
                f.write(f"    WaterBox = {waterBox_entry.get()}\n")
                f.write("    xBCtype = reflect\n")
                f.write("    yBCtype = reflect\n")
                f.write("    zBCtype = reflect\n")
            else:
                f.write("    isSphere = true\n")
                f.write(f"    sphereR  = {sphereR_entry.get()}\n")
            f.write("end boundaries\n\n")

            f.write("start molecules\n")
            for mol_data in added_molecules:
                f.write(f"    {mol_data['name']} : {mol_data['copy_number']}\n")
            f.write("end molecules\n\n")

            f.write("start reactions\n")
            for react_data in reactions:
                f.write(f"{react_data}\n")
            f.write("end reactions\n")
        messagebox.showinfo("File Generated", "The .inp file has been successfully generated.")

    def create_mol_file_from_entries():
        mol_name = mol_name_entry.get()
        check_overlap = checkOverlap_entry.get().lower() == "true"
        count_transition = countTransition_entry.get().lower() == "true"
        transition_matrix_size = int(transitionMatrixSize_entry.get())
        is_implicit_lipid = isImplicitLipid_entry.get().lower() == "true"
        if is_implicit_lipid and not added_molecules:
            messagebox.showerror("Error", "Please add the implicit lipid before other mols.")
            return
        is_lipid = isLipid_entry.get().lower() == "true"
        D = D_entry.get()
        Dr = Dr_entry.get()

        input_str = coords_entry.get()
        sites_list = input_str.split(';')
        sites = []
        coords_formatted = ''
        for site_str in sites_list:
            if ',' not in site_str:
                site_name, site_coord = site_str.split(':')
                x, y, z = site_coord.split()
                coords_formatted += site_name + '\t' + x + '\t' + y + '\t' + z + '\n'
                sites.append(Site(site_name))
            else:
                coord_str, state_str = site_str.split(',')
                coord_str = coord_str.strip()
                state_str = state_str.strip()
                site_name, site_coord = coord_str.split(':')
                x, y, z = site_coord.split()
                coords_formatted += site_name + '\t' + x + '\t' + y + '\t' + z + '\n'
                states = []
                states_str = state_str.split('~')
                for state in states_str:
                    states.append(State(state))
                sites.append(Site(site_name, states=states)) 


        bonds = [bond.strip() for bond in bonds_entry.get().split(';')]

        file_path = f"{mol_name}.mol"

        mols.append(Mol(mol_name, sites=sites))
        reactant1_entry["values"] = [mol.name for mol in mols]
        reactant2_entry["values"] = [mol.name for mol in mols]
        product1_entry["values"] = [mol.name for mol in mols]
        product2_entry["values"] = [mol.name for mol in mols]

        with open(file_path, "w") as file:
            file.write(f"Name = {mol_name}\n")
            file.write(f"checkOverlap = {check_overlap}\n")
            file.write(f"countTransition = {count_transition}\n")
            file.write(f"transitionMatrixSize = {transition_matrix_size}\n")
            if is_lipid:
                file.write("isLipid = true\n")
            if is_implicit_lipid:
                file.write("isImplicitLipid = true\n")
            file.write(f"D = {D}\n")
            file.write(f"Dr = {Dr}\n")
            file.write("COM    0.0000    0.0000    0.0000\n")

            file.write(f"{coords_formatted}\n")

            file.write("bonds = {}\n".format(len(bonds)))
            for bond in bonds:
                file.write(f"{bond}\n")

            for one_site in sites:
                if one_site.states:
                    all_states_str = f'state={one_site.name}'
                    for one_state in one_site.states:
                        all_states_str += '~' + one_state.name.upper()
                    file.write(f"{all_states_str}\n")

    def add_reaction():
        react_type = reaction_type_entry.get()
        if react_type == 'Zeroth Creation: null -> A':
            react_str = '    NULL -> '
            prod_name = product1_entry.get()
            prod_site = product1_site.get()
            prod_state = product1_state.get()
            r = rate_entry.get().strip()
            l = label_entry.get().strip()
            if prod_state == '':
                react_str += prod_name + '(' + prod_site + ')\n'
            else:
                react_str += prod_name + '(' + prod_site + '~' + prod_state + ')\n'
            react_str += '    rate = ' + r +'\n'
            if l != '':
                react_str += '    rxnLagbel = ' + l + '\n'
            reactions.append(react_str)

        elif react_type == 'Destruction: A -> null':
            react_str = '    '
            react_name = reactant1_entry.get()
            react_site = reactant1_site.get()
            react_state = reactant1_state.get()
            r = rate_entry.get().strip()
            l = label_entry.get().strip()
            if react_state == '':
                react_str += react_name + '(' + react_site + ') -> NULL\n'
            else:
                react_str += react_name + '(' + react_site + '~' + react_state + ') -> NULL\n'
            react_str += '    rate = ' + r +'\n'
            if l != '':
                react_str += '    rxnLagbel = ' + l + '\n'
            reactions.append(react_str)

        elif react_type == 'Statechange: B(b~U) -> B(b~P)':
            react_str = '    '
            react_name = reactant1_entry.get()
            react_site = reactant1_site.get()
            react_state = reactant1_state.get()
            prod_name = product1_entry.get()
            prod_site = product1_site.get()
            prod_state = product1_state.get()
            r = rate_entry.get().strip()
            l = label_entry.get().strip()
            if react_state == '':
                react_str += react_name + '(' + react_site + ') -> '
            else:
                react_str += react_name + '(' + react_site + '~' + react_state + ') -> '
            if prod_state == '':
                react_str += prod_name + '(' + prod_site + ')\n'
            else:
                react_str += prod_name + '(' + prod_site + '~' + prod_state + ')\n'
            react_str += '    rate = ' + r +'\n'
            if l != '':
                react_str += '    rxnLagbel = ' + l + '\n'
            reactions.append(react_str)

        elif react_type == 'Bimolecular Association: (A + B <-> A.B)':
            react1_name = reactant1_entry.get()
            react1_site = reactant1_site.get()
            react1_state = reactant1_state.get()
            react1_req = reactant1_bound_site_requirment.get()
            prod1_name = product1_entry.get()
            prod1_site = product1_site.get()
            prod1_state = product1_state.get()
            react2_name = reactant2_entry.get()
            react2_site = reactant2_site.get()
            react2_state = reactant2_state.get()
            react2_req = reactant2_bound_site_requirment.get()
            prod2_name = product2_entry.get()
            prod2_site = product2_site.get()
            prod2_state = product2_state.get()
            micro_on_rate_str = onRate3Dka_entry.get().strip()
            micro_off_rate_str = offRatekb_entry.get().strip()
            macro_on_rate_str = onRate3DMacro_entry.get().strip()
            macro_off_rate_str = offRateMacro_entry.get().strip()
            kcat_str = kcat_entry.get().strip()
            bindRadSameCom_str = bindRadSameCom_entry.get().strip()
            loopCoopFactor_str = loopCoopFactor_entry.get().strip()
            length3dto2d_str = length3dto2d_entry.get().strip()
            exclude_str = exclude_entry.get().strip()
            coupledLabel_str = coupledLabel_entry.get().strip()
            l = label_entry.get().strip()
            c1_str = coord_com_reactant1_entry.get().strip()
            c2_str = coord_com_reactant2_entry.get().strip()
            p1_str = coord_site_reactant1_entry.get().strip()
            p2_str = coord_site_reactant2_entry.get().strip()

            theta1 = 'nan'
            theta2 = 'nan'
            phi1 = 'nan'
            phi2 = 'nan'
            omega = 'nan'
            # get c1
            c1_list = c1_str[1:-1].split(',')
            c1 = np.array([float(one) for one in c1_list])
            # get c2
            c2_list = c2_str[1:-1].split(',')
            c2 = np.array([float(one) for one in c2_list])
            # get p1
            p1_list = p1_str[1:-1].split(',')
            p1 = np.array([float(one) for one in p1_list])
            # get p2
            p2_list = p2_str[1:-1].split(',')
            p2 = np.array([float(one) for one in p2_list])
            # get sigma1
            sigma1 = p1 - p2
            # get sigma2
            sigma2 = -sigma1
            # get sigma
            sigma = np.linalg.norm(sigma1)

            v1 = p1 - c1
            v2 = p2 - c2

            # get n1
            n1 = np.cross(v1, sigma1)
            if not np.isclose(np.linalg.norm(n1), 0):
                n1 = n1 / np.linalg.norm(n1)
            else:
                n1 = np.array([0,0,1])

            # get n2
            n2 = np.cross(v2, sigma2)
            if not np.isclose(np.linalg.norm(n2), 0):
                n2 = n2 / np.linalg.norm(n2)
            else:
                n2 = np.array([0,0,1])

            # get theta1, theta2
            if not np.isclose(np.linalg.norm(v1), 0):
                theta1 = np.arccos(np.dot(v1, sigma1) / (np.linalg.norm(v1) * np.linalg.norm(sigma1)))
            else:
                theta1 = 'nan'
            if not np.isclose(np.linalg.norm(v2), 0):
                theta2 = np.arccos(np.dot(v2, sigma2) / (np.linalg.norm(v2) * np.linalg.norm(sigma2)))
            else:
                theta2 = 'nan'

            # get phi1, phi2
            if not np.isclose(np.linalg.norm(np.cross(v1, sigma1)), 0):
                t1 = np.cross(v1, sigma1)
                t2 = np.cross(v1, n1)
                norm_t1 = t1 / np.linalg.norm(t1)
                norm_t2 = t2 / np.linalg.norm(t2)
                phi1 = np.arccos(np.dot(norm_t1, norm_t2))
            else:
                phi1 = 'nan'
            
            if not np.isclose(np.linalg.norm(np.cross(v2, sigma2)), 0):
                t1 = np.cross(v2, sigma2)
                t2 = np.cross(v2, n2)
                norm_t1 = t1 / np.linalg.norm(t1)
                norm_t2 = t2 / np.linalg.norm(t2)
                phi2 = np.arccos(np.dot(norm_t1, norm_t2))
            else:
                phi1 = 'nan'

            # get omega
            if not np.isclose(np.linalg.norm(v1), 0) and not np.isclose(np.linalg.norm(v2), 0):
                if not np.isclose(np.linalg.norm(np.cross(v1, sigma1)), 0) and not np.isclose(np.linalg.norm(np.cross(v2, sigma2)), 0):
                    t1 = np.cross(sigma1, v1)
                    t2 = np.cross(sigma1, v2)
                else:
                    t1 = np.cross(sigma1, n1)
                    t2 = np.cross(sigma1, n2)

                omega = np.arccos(np.dot(t1, t2) / (np.linalg.norm(t1) * np.linalg.norm(t2)))
            else:
                omega = 'nan'

            react_str = '    '
            if react1_state == '':
                if react1_req == '':
                    react1_str = react1_name + '(' + react1_site + ')'
                    prod1_str = react1_name + '(' + react1_site + '!1)'
                else:
                    react1_str = f"{react1_name}({react1_site}, {react1_req}!*)"
                    prod1_str = f"{react1_name}({react1_site}!1, {react1_req}!*)"
            else:
                if react1_req == '':
                    react1_str = react1_name + '(' + react1_site + '~' + react1_state + ')'
                    prod1_str = react1_name + '(' + react1_site + '~' + react1_state + '!1)'
                else:
                    react1_str = f"{react1_name}({react1_site}~{react1_state}, {react1_req}!*)"
                    prod1_str = f"{react1_name}({react1_site}~{react1_state}!1, {react1_req}!*)"
            if react2_state == '':
                if react2_req == '':
                    react2_str = react2_name + '(' + react2_site + ')'
                    prod2_str = react2_name + '(' + react2_site + '!1)'
                else:
                    react2_str = f"{react2_name}({react2_site}, {react2_req}!*)"
                    prod2_str = f"{react2_name}({react2_site}!1, {react2_req}!*)"
            else:
                if react2_req == '':
                    react2_str = react2_name + '(' + react2_site + '~' + react2_state + ')'
                    prod2_str = react2_name + '(' + react2_site + '~' + react2_state + '!1)'
                else:
                    react2_str = f"{react2_name}({react2_site}~{react2_state}, {react2_req}!*)"
                    prod2_str = f"{react2_name}({react2_site}~{react2_state}!1, {react2_req}!*)"

            prod_str = prod1_str + '.' + prod2_str

            react_str += f"{react1_str} + {react2_str} <-> {prod_str}\n" 
            
            if micro_on_rate_str != '':
                if micro_off_rate_str == '':
                    messagebox.showerror("Error", "Please enter micro off rate.")
                    return
                else:
                    react_str += f'    onRate3Dka = {micro_on_rate_str}\n'
                    react_str += f'    offRatekb = {micro_off_rate_str}\n'
            elif macro_on_rate_str != '':
                if macro_off_rate_str == '':
                    messagebox.showerror("Error", "Please enter macro off rate.")
                    return
                else:
                    react_str += f'    onRate3DMacro = {macro_on_rate_str}\n'
                    react_str += f'    offRateMacro = {macro_off_rate_str}\n'
            else:
                messagebox.showerror("Error", "Please enter reaction rates (either micro or macro).")
                return
            
            if kcat_str != '':
                react_str += f'    kcat = {kcat_str}\n'

            react_str += f'    sigma = {sigma}\n'

            if bindRadSameCom_str != '':
                react_str += f'    bindRadSameCom = {bindRadSameCom_str}\n'

            if loopCoopFactor_str != '':
                react_str += f'    loopCoopFactor = {loopCoopFactor_str}\n'

            if length3dto2d_str != '':
                react_str += f'    length3Dto2D = {length3dto2d_str}\n'

            react_str += f'    norm1 = {str(list(n1))}\n'
            react_str += f'    norm2 = {str(list(n2))}\n'
            react_str += '    assocAngles = [' + str(theta1) + ',' + str(theta2) + ',' + str(phi1) + ',' + str(phi2) + ',' + str(omega) + ']\n'

            if l != '':
                react_str += '    rxnLagbel = ' + l + '\n'

            if coupledLabel_str != '':
                react_str += '    coupledRxnLabel = ' + coupledLabel_str + '\n'

            if exclude_str != '':
                exclude_value = exclude_str.lower() == 'true'
                react_str += f'    excludeVolumeBound = {exclude_value}\n'
            
            reactions.append(react_str)

        elif react_type == 'Michaelis-Menten: A + B <-> A.B -> A + C':
            react1_name = reactant1_entry.get()
            react1_site = reactant1_site.get()
            react1_state = reactant1_state.get()
            react1_req = reactant1_bound_site_requirment.get()
            prod1_name = product1_entry.get()
            prod1_site = product1_site.get()
            prod1_state = product1_state.get()
            react2_name = reactant2_entry.get()
            react2_site = reactant2_site.get()
            react2_state = reactant2_state.get()
            react2_req = reactant2_bound_site_requirment.get()
            prod2_name = product2_entry.get()
            prod2_site = product2_site.get()
            prod2_state = product2_state.get()
            micro_on_rate_str = onRate3Dka_entry.get().strip()
            micro_off_rate_str = offRatekb_entry.get().strip()
            macro_on_rate_str = onRate3DMacro_entry.get().strip()
            macro_off_rate_str = offRateMacro_entry.get().strip()
            kcat_str = kcat_entry.get().strip()
            bindRadSameCom_str = bindRadSameCom_entry.get().strip()
            loopCoopFactor_str = loopCoopFactor_entry.get().strip()
            length3dto2d_str = length3dto2d_entry.get().strip()
            exclude_str = exclude_entry.get().strip()
            coupledLabel_str = coupledLabel_entry.get().strip()
            l = label_entry.get().strip()
            c1_str = coord_com_reactant1_entry.get().strip()
            c2_str = coord_com_reactant2_entry.get().strip()
            p1_str = coord_site_reactant1_entry.get().strip()
            p2_str = coord_site_reactant2_entry.get().strip()

            theta1 = 'nan'
            theta2 = 'nan'
            phi1 = 'nan'
            phi2 = 'nan'
            omega = 'nan'
            # get c1
            c1_list = c1_str[1:-1].split(',')
            c1 = np.array([float(one) for one in c1_list])
            # get c2
            c2_list = c2_str[1:-1].split(',')
            c2 = np.array([float(one) for one in c2_list])
            # get p1
            p1_list = p1_str[1:-1].split(',')
            p1 = np.array([float(one) for one in p1_list])
            # get p2
            p2_list = p2_str[1:-1].split(',')
            p2 = np.array([float(one) for one in p2_list])
            # get sigma1
            sigma1 = p1 - p2
            # get sigma2
            sigma2 = -sigma1
            # get sigma
            sigma = np.linalg.norm(sigma1)

            v1 = p1 - c1
            v2 = p2 - c2

            # get n1
            n1 = np.cross(v1, sigma1)
            if not np.isclose(np.linalg.norm(n1), 0):
                n1 = n1 / np.linalg.norm(n1)
            else:
                n1 = np.array([0,0,1])

            # get n2
            n2 = np.cross(v2, sigma2)
            if not np.isclose(np.linalg.norm(n2), 0):
                n2 = n2 / np.linalg.norm(n2)
            else:
                n2 = np.array([0,0,1])

            # get theta1, theta2
            if not np.isclose(np.linalg.norm(v1), 0):
                theta1 = np.arccos(np.dot(v1, sigma1) / (np.linalg.norm(v1) * np.linalg.norm(sigma1)))
            else:
                theta1 = 'nan'
            if not np.isclose(np.linalg.norm(v2), 0):
                theta2 = np.arccos(np.dot(v2, sigma2) / (np.linalg.norm(v2) * np.linalg.norm(sigma2)))
            else:
                theta2 = 'nan'

            # get phi1, phi2
            if not np.isclose(np.linalg.norm(np.cross(v1, sigma1)), 0):
                t1 = np.cross(v1, sigma1)
                t2 = np.cross(v1, n1)
                norm_t1 = t1 / np.linalg.norm(t1)
                norm_t2 = t2 / np.linalg.norm(t2)
                phi1 = np.arccos(np.dot(norm_t1, norm_t2))
            else:
                phi1 = 'nan'
            
            if not np.isclose(np.linalg.norm(np.cross(v2, sigma2)), 0):
                t1 = np.cross(v2, sigma2)
                t2 = np.cross(v2, n2)
                norm_t1 = t1 / np.linalg.norm(t1)
                norm_t2 = t2 / np.linalg.norm(t2)
                phi2 = np.arccos(np.dot(norm_t1, norm_t2))
            else:
                phi1 = 'nan'

            # get omega
            if not np.isclose(np.linalg.norm(v1), 0) and not np.isclose(np.linalg.norm(v2), 0):
                if not np.isclose(np.linalg.norm(np.cross(v1, sigma1)), 0) and not np.isclose(np.linalg.norm(np.cross(v2, sigma2)), 0):
                    t1 = np.cross(sigma1, v1)
                    t2 = np.cross(sigma1, v2)
                else:
                    t1 = np.cross(sigma1, n1)
                    t2 = np.cross(sigma1, n2)

                omega = np.arccos(np.dot(t1, t2) / (np.linalg.norm(t1) * np.linalg.norm(t2)))
            else:
                omega = 'nan'

            react_str = '    '
            if react1_state == '':
                if react1_req == '':
                    react1_str = react1_name + '(' + react1_site + ')'
                    prod1_str = react1_name + '(' + react1_site + '!1)'
                else:
                    react1_str = f"{react1_name}({react1_site}, {react1_req}!*)"
                    prod1_str = f"{react1_name}({react1_site}!1, {react1_req}!*)"
            else:
                if react1_req == '':
                    react1_str = react1_name + '(' + react1_site + '~' + react1_state + ')'
                    prod1_str = react1_name + '(' + react1_site + '~' + react1_state + '!1)'
                else:
                    react1_str = f"{react1_name}({react1_site}~{react1_state}, {react1_req}!*)"
                    prod1_str = f"{react1_name}({react1_site}~{react1_state}!1, {react1_req}!*)"
            if react2_state == '':
                if react2_req == '':
                    react2_str = react2_name + '(' + react2_site + ')'
                    prod2_str = react2_name + '(' + react2_site + '!1)'
                else:
                    react2_str = f"{react2_name}({react2_site}, {react2_req}!*)"
                    prod2_str = f"{react2_name}({react2_site}!1, {react2_req}!*)"
            else:
                if react2_req == '':
                    react2_str = react2_name + '(' + react2_site + '~' + react2_state + ')'
                    prod2_str = react2_name + '(' + react2_site + '~' + react2_state + '!1)'
                else:
                    react2_str = f"{react2_name}({react2_site}~{react2_state}, {react2_req}!*)"
                    prod2_str = f"{react2_name}({react2_site}~{react2_state}!1, {react2_req}!*)"

            prod_str = prod1_str + '.' + prod2_str

            react_str += f"{react1_str} + {react2_str} <-> {prod_str}\n" 
            
            if micro_on_rate_str != '':
                if micro_off_rate_str == '':
                    messagebox.showerror("Error", "Please enter micro off rate.")
                    return
                else:
                    react_str += f'    onRate3Dka = {micro_on_rate_str}\n'
                    react_str += f'    offRatekb = {micro_off_rate_str}\n'
            elif macro_on_rate_str != '':
                if macro_off_rate_str == '':
                    messagebox.showerror("Error", "Please enter macro off rate.")
                    return
                else:
                    react_str += f'    onRate3DMacro = {macro_on_rate_str}\n'
                    react_str += f'    offRateMacro = {macro_off_rate_str}\n'
            else:
                messagebox.showerror("Error", "Please enter reaction rates (either micro or macro).")
                return
            
            if kcat_str != '':
                react_str += f'    kcat = {kcat_str}\n'
            else:
                messagebox.showerror("Error", "Please specify kcat for Michaelis-Menten.")
                return

            react_str += f'    sigma = {sigma}\n'

            if bindRadSameCom_str != '':
                react_str += f'    bindRadSameCom = {bindRadSameCom_str}\n'

            if loopCoopFactor_str != '':
                react_str += f'    loopCoopFactor = {loopCoopFactor_str}\n'

            if length3dto2d_str != '':
                react_str += f'    length3Dto2D = {length3dto2d_str}\n'

            react_str += f'    norm1 = {str(list(n1))}\n'
            react_str += f'    norm2 = {str(list(n2))}\n'
            react_str += '    assocAngles = [' + str(theta1) + ',' + str(theta2) + ',' + str(phi1) + ',' + str(phi2) + ',' + str(omega) + ']\n'

            if l != '':
                react_str += '    rxnLagbel = ' + l + '\n'

            if coupledLabel_str != '':
                react_str += '    coupledRxnLabel = ' + coupledLabel_str + '\n'
            else:
                messagebox.showerror("Error", "Please specify coupled reaction for Michaelis-Menten.")
                return

            if exclude_str != '':
                exclude_value = exclude_str.lower() == 'true'
                react_str += f'    excludeVolumeBound = {exclude_value}\n'
            
            reactions.append(react_str)

        elif react_type == 'Bimolecular Statechange: A + B(b~U) -> A + B(b~P)':
            react1_name = reactant1_entry.get()
            react1_site = reactant1_site.get()
            react1_state = reactant1_state.get()
            react1_req = reactant1_bound_site_requirment.get()
            prod1_name = product1_entry.get()
            prod1_site = product1_site.get()
            prod1_state = product1_state.get()
            prod1_req = product1_bound_site_requirment.get()
            react2_name = reactant2_entry.get()
            react2_site = reactant2_site.get()
            react2_state = reactant2_state.get()
            react2_req = reactant2_bound_site_requirment.get()
            prod2_name = product2_entry.get()
            prod2_site = product2_site.get()
            prod2_state = product2_state.get()
            prod2_req = product2_bound_site_requirment.get()
            micro_on_rate_str = onRate3Dka_entry.get().strip()
            micro_off_rate_str = offRatekb_entry.get().strip()
            macro_on_rate_str = onRate3DMacro_entry.get().strip()
            macro_off_rate_str = offRateMacro_entry.get().strip()
            kcat_str = kcat_entry.get().strip()
            bindRadSameCom_str = bindRadSameCom_entry.get().strip()
            loopCoopFactor_str = loopCoopFactor_entry.get().strip()
            length3dto2d_str = length3dto2d_entry.get().strip()
            exclude_str = exclude_entry.get().strip()
            coupledLabel_str = coupledLabel_entry.get().strip()
            l = label_entry.get().strip()
            c1_str = coord_com_reactant1_entry.get().strip()
            c2_str = coord_com_reactant2_entry.get().strip()
            p1_str = coord_site_reactant1_entry.get().strip()
            p2_str = coord_site_reactant2_entry.get().strip()

            theta1 = 'nan'
            theta2 = 'nan'
            phi1 = 'nan'
            phi2 = 'nan'
            omega = 'nan'
            # get c1
            c1_list = c1_str[1:-1].split(',')
            c1 = np.array([float(one) for one in c1_list])
            # get c2
            c2_list = c2_str[1:-1].split(',')
            c2 = np.array([float(one) for one in c2_list])
            # get p1
            p1_list = p1_str[1:-1].split(',')
            p1 = np.array([float(one) for one in p1_list])
            # get p2
            p2_list = p2_str[1:-1].split(',')
            p2 = np.array([float(one) for one in p2_list])
            # get sigma1
            sigma1 = p1 - p2
            # get sigma2
            sigma2 = -sigma1
            # get sigma
            sigma = np.linalg.norm(sigma1)

            v1 = p1 - c1
            v2 = p2 - c2

            # get n1
            n1 = np.cross(v1, sigma1)
            if not np.isclose(np.linalg.norm(n1), 0):
                n1 = n1 / np.linalg.norm(n1)
            else:
                n1 = np.array([0,0,1])

            # get n2
            n2 = np.cross(v2, sigma2)
            if not np.isclose(np.linalg.norm(n2), 0):
                n2 = n2 / np.linalg.norm(n2)
            else:
                n2 = np.array([0,0,1])

            # get theta1, theta2
            if not np.isclose(np.linalg.norm(v1), 0):
                theta1 = np.arccos(np.dot(v1, sigma1) / (np.linalg.norm(v1) * np.linalg.norm(sigma1)))
            else:
                theta1 = 'nan'
            if not np.isclose(np.linalg.norm(v2), 0):
                theta2 = np.arccos(np.dot(v2, sigma2) / (np.linalg.norm(v2) * np.linalg.norm(sigma2)))
            else:
                theta2 = 'nan'

            # get phi1, phi2
            if not np.isclose(np.linalg.norm(np.cross(v1, sigma1)), 0):
                t1 = np.cross(v1, sigma1)
                t2 = np.cross(v1, n1)
                norm_t1 = t1 / np.linalg.norm(t1)
                norm_t2 = t2 / np.linalg.norm(t2)
                phi1 = np.arccos(np.dot(norm_t1, norm_t2))
            else:
                phi1 = 'nan'
            
            if not np.isclose(np.linalg.norm(np.cross(v2, sigma2)), 0):
                t1 = np.cross(v2, sigma2)
                t2 = np.cross(v2, n2)
                norm_t1 = t1 / np.linalg.norm(t1)
                norm_t2 = t2 / np.linalg.norm(t2)
                phi2 = np.arccos(np.dot(norm_t1, norm_t2))
            else:
                phi1 = 'nan'

            # get omega
            if not np.isclose(np.linalg.norm(v1), 0) and not np.isclose(np.linalg.norm(v2), 0):
                if not np.isclose(np.linalg.norm(np.cross(v1, sigma1)), 0) and not np.isclose(np.linalg.norm(np.cross(v2, sigma2)), 0):
                    t1 = np.cross(sigma1, v1)
                    t2 = np.cross(sigma1, v2)
                else:
                    t1 = np.cross(sigma1, n1)
                    t2 = np.cross(sigma1, n2)

                omega = np.arccos(np.dot(t1, t2) / (np.linalg.norm(t1) * np.linalg.norm(t2)))
            else:
                omega = 'nan'

            react_str = '    '
            if react1_state == '':
                if react1_req == '':
                    react1_str = react1_name + '(' + react1_site + ')'
                else:
                    react1_str = f"{react1_name}({react1_site}, {react1_req}!*)"
            else:
                if react1_req == '':
                    react1_str = react1_name + '(' + react1_site + '~' + react1_state + ')'
                else:
                    react1_str = f"{react1_name}({react1_site}~{react1_state}, {react1_req}!*)"
            if react2_state == '':
                if react2_req == '':
                    react2_str = react2_name + '(' + react2_site + ')'
                else:
                    react2_str = f"{react2_name}({react2_site}, {react2_req}!*)"
            else:
                if react2_req == '':
                    react2_str = react2_name + '(' + react2_site + '~' + react2_state + ')'
                else:
                    react2_str = f"{react2_name}({react2_site}~{react2_state}, {react2_req}!*)"

            if prod1_state == '':
                if prod1_req == '':
                    prod1_str = prod1_name + '(' + prod1_site + ')'
                else:
                    prod1_str = f"{prod1_name}({prod1_site}, {prod1_req}!*)"
            else:
                if prod1_req == '':
                    prod1_str = prod1_name + '(' + prod1_site + '~' + prod1_state + ')'
                else:
                    prod1_str = f"{prod1_name}({prod1_site}~{prod1_state}, {prod1_req}!*)"
            if prod2_state == '':
                if prod2_req == '':
                    prod2_str = prod2_name + '(' + prod2_site + ')'
                else:
                    prod2_str = f"{prod2_name}({prod2_site}, {prod2_req}!*)"
            else:
                if prod2_req == '':
                    prod2_str = prod2_name + '(' + prod2_site + '~' + prod2_state + ')'
                else:
                    prod2_str = f"{prod2_name}({prod2_site}~{prod2_state}, {prod2_req}!*)"

            react_str += f"{react1_str} + {react2_str} <-> {prod1_str} + {prod2_str}\n" 
            
            if micro_on_rate_str != '':
                if micro_off_rate_str == '':
                    messagebox.showerror("Error", "Please enter micro off rate.")
                    return
                else:
                    react_str += f'    onRate3Dka = {micro_on_rate_str}\n'
                    react_str += f'    offRatekb = {micro_off_rate_str}\n'
            elif macro_on_rate_str != '':
                if macro_off_rate_str == '':
                    messagebox.showerror("Error", "Please enter macro off rate.")
                    return
                else:
                    react_str += f'    onRate3DMacro = {macro_on_rate_str}\n'
                    react_str += f'    offRateMacro = {macro_off_rate_str}\n'
            else:
                messagebox.showerror("Error", "Please enter reaction rates (either micro or macro).")
                return
            
            if kcat_str != '':
                react_str += f'    kcat = {kcat_str}\n'

            react_str += f'    sigma = {sigma}\n'

            if bindRadSameCom_str != '':
                react_str += f'    bindRadSameCom = {bindRadSameCom_str}\n'

            if loopCoopFactor_str != '':
                react_str += f'    loopCoopFactor = {loopCoopFactor_str}\n'

            if length3dto2d_str != '':
                react_str += f'    length3Dto2D = {length3dto2d_str}\n'

            react_str += f'    norm1 = {str(list(n1))}\n'
            react_str += f'    norm2 = {str(list(n2))}\n'
            react_str += '    assocAngles = [' + str(theta1) + ',' + str(theta2) + ',' + str(phi1) + ',' + str(phi2) + ',' + str(omega) + ']\n'

            if l != '':
                react_str += '    rxnLagbel = ' + l + '\n'

            if coupledLabel_str != '':
                react_str += '    coupledRxnLabel = ' + coupledLabel_str + '\n'

            if exclude_str != '':
                exclude_value = exclude_str.lower() == 'true'
                react_str += f'    excludeVolumeBound = {exclude_value}\n'
            
            reactions.append(react_str)

        elif react_type == 'Unimolecular Creation: A -> A + B':
            react1_name = reactant1_entry.get()
            react1_site = reactant1_site.get()
            react1_state = reactant1_state.get()
            react1_req = reactant1_bound_site_requirment.get()
            prod1_name = product1_entry.get()
            prod1_site = product1_site.get()
            prod1_state = product1_state.get()
            prod1_req = product1_bound_site_requirment.get()
            react2_name = reactant2_entry.get()
            react2_site = reactant2_site.get()
            react2_state = reactant2_state.get()
            react2_req = reactant2_bound_site_requirment.get()
            prod2_name = product2_entry.get()
            prod2_site = product2_site.get()
            prod2_state = product2_state.get()
            prod2_req = product2_bound_site_requirment.get()
            micro_on_rate_str = onRate3Dka_entry.get().strip()
            micro_off_rate_str = offRatekb_entry.get().strip()
            macro_on_rate_str = onRate3DMacro_entry.get().strip()
            macro_off_rate_str = offRateMacro_entry.get().strip()
            kcat_str = kcat_entry.get().strip()
            bindRadSameCom_str = bindRadSameCom_entry.get().strip()
            loopCoopFactor_str = loopCoopFactor_entry.get().strip()
            length3dto2d_str = length3dto2d_entry.get().strip()
            exclude_str = exclude_entry.get().strip()
            coupledLabel_str = coupledLabel_entry.get().strip()
            l = label_entry.get().strip()
            r = rate_entry.get().strip()

            react_str = '    '
            if react1_state == '':
                if react1_req == '':
                    react1_str = react1_name + '(' + react1_site + ')'
                else:
                    react1_str = f"{react1_name}({react1_site}, {react1_req}!*)"
            else:
                if react1_req == '':
                    react1_str = react1_name + '(' + react1_site + '~' + react1_state + ')'
                else:
                    react1_str = f"{react1_name}({react1_site}~{react1_state}, {react1_req}!*)"
            
            if prod1_state == '':
                if prod1_req == '':
                    prod1_str = prod1_name + '(' + prod1_site + ')'
                else:
                    prod1_str = f"{prod1_name}({prod1_site}, {prod1_req}!*)"
            else:
                if prod1_req == '':
                    prod1_str = prod1_name + '(' + prod1_site + '~' + prod1_state + ')'
                else:
                    prod1_str = f"{prod1_name}({prod1_site}~{prod1_state}, {prod1_req}!*)"
            if prod2_state == '':
                if prod2_req == '':
                    prod2_str = prod2_name + '(' + prod2_site + ')'
                else:
                    prod2_str = f"{prod2_name}({prod2_site}, {prod2_req}!*)"
            else:
                if prod2_req == '':
                    prod2_str = prod2_name + '(' + prod2_site + '~' + prod2_state + ')'
                else:
                    prod2_str = f"{prod2_name}({prod2_site}~{prod2_state}, {prod2_req}!*)"

            react_str += f"{react1_str} <-> {prod1_str} + {prod2_str}\n" 
            
            
            if r != '':
                react_str += f'    rate = {r}\n'


            if l != '':
                react_str += '    rxnLagbel = ' + l + '\n'

            if coupledLabel_str != '':
                react_str += '    coupledRxnLabel = ' + coupledLabel_str + '\n'
            
            reactions.append(react_str)

        messagebox.showinfo("Reaction Added", "Reaction Added.")
        return

    def add_molecule():
        mol_name = mol_name_entry.get().strip()
        copy_num = mol_copy_num_entry.get().strip()
        concentration = mol_concentration_entry.get().strip()

        if not mol_name:
            messagebox.showerror("Error", "Please enter a molecule name.")
            return

        if not copy_num and not concentration:
            messagebox.showerror("Error", "Please enter either copy number or concentration.")
            return

        if copy_num and concentration:
            messagebox.showerror("Error", f"Please enter either copy number {copy_num} OR concentration {concentration}, not both.")
            return

        if concentration:
            # Convert concentration to copy number
            isSphere_value = isSphere_entry.get().strip().lower() == "true"
            volume = 0
            if not isSphere_value:
                dimensions_str = waterBox_entry.get().strip()
                dimensions = [int(dim) for dim in dimensions_str.strip("[]").split(",")]
                x, y, z = dimensions
                volume = x * y * z / 1E9 # um3
            else:
                r = int(sphereR_entry.get().strip())
                volume = 4 / 3 * 3.1415927 * r ** 3 / 1E9 # um3
            copy_num = int(float(concentration) / 1000000 * volume /1000000000000000 * 6.02E+23)

        percentage = states_entry.get().strip()
        if not percentage:
            mol_data = {
                "name": mol_name,
                "copy_number": copy_num
            }
            added_molecules.append(mol_data)
            molecules_listbox.insert(tk.END, f"{mol_name} : {copy_num}")
        else:
            state_str_list = percentage.split(',')
            new_copy_number_str_list = []
            for one_state in state_str_list:
                num, name = one_state.split()
                num = int(float(num) * int(copy_num))
                new_copy_number_str_list.append(str(num) + ' ' + name)
                new_copy_number_str = ','.join(new_copy_number_str_list)
            mol_data = {
                "name": mol_name,
                "copy_number":new_copy_number_str
            }
            added_molecules.append(mol_data)
            molecules_listbox.insert(tk.END, f"{mol_name} : {new_copy_number_str}")

        create_mol_file_from_entries()
        messagebox.showinfo("Molecule Added", f"The {mol_name}.mol file has been successfully generated.")

    def update_combobox_options(combobox, options):
        combobox["values"] = options
        combobox.set('')  # Clear the current selection
            
    def update_reaction_site_options(selected_molecule, sites_combobox):
        sites_options = []
        for mol in mols:
            if mol.name == selected_molecule:
                for site in mol.sites:
                    sites_options.append(site.name)
        if sites_options:
            update_combobox_options(sites_combobox, sites_options)
        else:
            update_combobox_options(sites_combobox, [])

    def update_reaction_state_options(selected_site, states_combobox):
        states_options = []
        for mol in mols:
            for site in mol.sites:
                if site.name == selected_site:
                    for state in site.states:
                        states_options.append(state.name)
        if states_options:
            update_combobox_options(states_combobox, states_options)
        else:
            update_combobox_options(states_combobox, [])

    def on_reactant1_selected(event):
        selected_mol = reactant1_entry.get()
        update_reaction_site_options(selected_mol, reactant1_site)
        update_reaction_site_options(selected_mol, reactant1_bound_site_requirment)

    def on_reactant2_selected(event):
        selected_mol = reactant2_entry.get()
        update_reaction_site_options(selected_mol, reactant2_site)
        update_reaction_site_options(selected_mol, reactant2_bound_site_requirment)

    def on_product1_selected(event):
        selected_mol = product1_entry.get()
        update_reaction_site_options(selected_mol, product1_site)

    def on_product2_selected(event):
        selected_mol = product2_entry.get()
        update_reaction_site_options(selected_mol, product2_site)

    def on_reactant1_site_selected(event):
        selected_site = reactant1_site.get()
        update_reaction_state_options(selected_site, reactant1_state)

    def on_reactant2_site_selected(event):
        selected_site = reactant2_site.get()
        update_reaction_state_options(selected_site, reactant2_state)

    def on_product1_site_selected(event):
        selected_site = product1_site.get()
        update_reaction_state_options(selected_site, product1_state)

    def on_product2_site_selected(event):
        selected_site = product2_site.get()
        update_reaction_state_options(selected_site, product2_state)



    root = tk.Tk()
    root.title("Simulation Parameters")

    # Create a tabbed interface
    tab_parent = ttk.Notebook(root)

    # Parameters tab
    parameters_tab = ttk.Frame(tab_parent)
    tab_parent.add(parameters_tab, text="Parameters")

    nItr_entry = create_labeled_entry_with_hint(parameters_tab, "Enter total iteration (steps):", "1000", 0, 0)
    timeStep_entry = create_labeled_entry_with_hint(parameters_tab, "Enter timeStep (us):", "0.1", 1, 0)
    timeWrite_entry = create_labeled_entry_with_hint(parameters_tab, "Enter how often to write output (steps):", "10", 2, 0)
    trajWrite_entry = create_labeled_entry_with_hint(parameters_tab, "Enter how often to write traj (steps):", "100", 3, 0)
    pdbWrite_entry = create_labeled_entry_with_hint(parameters_tab, "Enter how often to write pdb (steps):", "100", 4, 0)
    transitionWrite_entry = create_labeled_entry_with_hint(parameters_tab, "Enter how often to write transition matrix (steps):", "100", 5, 0)
    restartWrite_entry = create_labeled_entry_with_hint(parameters_tab, "Enter how often to write restart file (steps):", "100", 6, 0)
    checkPoint_entry = create_labeled_entry_with_hint(parameters_tab, "Enter how often to write check points (steps):", "100", 7, 0)
    scaleMaxDisplace_entry = create_labeled_entry_with_hint(parameters_tab, "Enter threshold to reject association (<RMSD>):", "100", 8, 0)
    overlapSepLimit_entry = create_labeled_entry_with_hint(parameters_tab, "Enter threshold for overlap check between COMs (nm):", "0.1", 9, 0)
    clusterOverlapCheck_entry = create_labeled_entry_with_hint(parameters_tab, "Enter whether overlap is checked based on cluster (bool):", "False", 10, 0)


    # Boundaries tab
    boundaries_tab = ttk.Frame(tab_parent)
    tab_parent.add(boundaries_tab, text="Boundaries")

    isSphere_entry = create_labeled_entry_with_hint(boundaries_tab, "Enter if the system is spherical (bool):", "False", 0, 0)
    waterBox_entry = create_labeled_entry_with_hint(boundaries_tab, "Enter the dimensions of the system if it is a box [x,y,z] (nm):", "[100,100,100]", 1, 0)
    sphereR_entry = create_labeled_entry_with_hint(boundaries_tab, "Enter the radius of the sphere if the system is a sphere (nm):", "100", 2, 0)

    # Molecules tab
    molecules_tab = ttk.Frame(tab_parent)
    tab_parent.add(molecules_tab, text="Molecules")
    # Create widgets for Molecules tab
    mol_name_entry = create_labeled_entry_with_hint(molecules_tab, "Enter molecule name:", "", 0, 0)
    mol_copy_num_entry = create_labeled_entry_with_hint(molecules_tab, "Enter copy number:", "", 1, 0)
    mol_concentration_entry = create_labeled_entry_with_hint(molecules_tab, "Enter concentration if copy number is not provided (uM):", "", 2, 0)
    checkOverlap_entry = create_labeled_entry_with_hint(molecules_tab, "Enter boolean value to check overlap for this molecule type:", "True", 3, 0)
    countTransition_entry = create_labeled_entry_with_hint(molecules_tab, "Enter boolean value to track transition matrix for this molecule type::", "True", 4, 0)
    transitionMatrixSize_entry = create_labeled_entry_with_hint(molecules_tab, "Enter transition matrix size for this molecule type, which should be larger than the size of the formed largest complex:", "500", 5, 0)
    isImplicitLipid_entry = create_labeled_entry_with_hint(molecules_tab, "Enter if it is isImplicitLipid (bool):", "False", 6, 0)
    isLipid_entry = create_labeled_entry_with_hint(molecules_tab, "Enter if it is lipid (bool):", "False", 7, 0)
    D_entry = create_labeled_entry_with_hint(molecules_tab, "Enter translational diffusion constants ([Dx,Dy,Dz] (um2/s)):", "[10,10,10]", 8, 0)
    Dr_entry = create_labeled_entry_with_hint(molecules_tab, "Enter rotational diffusion constants ([Drx,Dry,Drz] (rad2/us)):", "[0.1,0,1,0.1]", 9, 0)
    coords_entry = create_labeled_entry_with_hint(molecules_tab, "Enter coords and states of all binding sites, states (one upper letter) are optional (c1: 1 1 1, U~P;c2: 2 2 2; c3: 3 3 3) (nm):", "", 10, 0)
    bonds_entry = create_labeled_entry_with_hint(molecules_tab, "Enter bonds between sites (com c1;com c2;com c3):", "", 11, 0)
    states_entry = create_labeled_entry_with_hint(molecules_tab, "Enter percentage of different initial states (0.8 (sitename~statename1), 0.2 (sitename~statename2)):", "", 12, 0)
    add_mol_button = tk.Button(molecules_tab, text="Add Molecule", command=add_molecule)
    molecules_listbox = tk.Listbox(molecules_tab)

    add_mol_button.grid(row=13, column=0, columnspan=2)
    molecules_listbox.grid(row=14, column=0, columnspan=2)

    # List to store added molecules
    added_molecules = []

    # all the mols in the system
    mols = []

    # Reactions tab
    reactions_tab = ttk.Frame(tab_parent)
    tab_parent.add(reactions_tab, text="Reactions")

    mol_names = [mol.name for mol in mols]
    reactant1_entry = create_labeled_combobox_entry_with_callback(reactions_tab, "Select reactant 1 if it is not creation reaction:", mol_names, 0, 0, on_reactant1_selected, width=10)
    reactant1_site = create_labeled_combobox_entry_with_callback(reactions_tab, "Select reactant 1 site:", [], 0, 2, on_reactant1_site_selected, width=10)
    reactant1_state = create_labeled_combobox_entry(reactions_tab, "Select reactant 1 state (optional):", [], 1, 0, width=10)
    reactant1_bound_site_requirment = create_labeled_combobox_entry(reactions_tab, "Select the site of reactant 1 that needs to be already \nbound for the reaction to occur (optional):", [], 1, 2, width=10)
    reactant2_entry = create_labeled_combobox_entry_with_callback(reactions_tab, "Select reactant 2 if it is not creation reaction:", mol_names, 2, 0, on_reactant2_selected, width=10)
    reactant2_site = create_labeled_combobox_entry_with_callback(reactions_tab, "Select reactant 2 site:", [], 2, 2, on_reactant2_site_selected, width=10)
    reactant2_state = create_labeled_combobox_entry(reactions_tab, "Select reactant 2 state (optional):", [], 3, 0, width=10)
    reactant2_bound_site_requirment = create_labeled_combobox_entry(reactions_tab, "Select the site of reactant2 that needs to be already \nbound for the reaction to occur (optional):", [], 3, 2, width=10)
    product1_entry = create_labeled_combobox_entry_with_callback(reactions_tab, "Select product 1 if it is not bimolecular reaction or destruction reaction:", mol_names, 4, 0, on_product1_selected, width=10)
    product1_site = create_labeled_combobox_entry_with_callback(reactions_tab, "Select product 1 site:", [], 4, 2, on_product1_site_selected, width=10)
    product1_bound_site_requirment = create_labeled_combobox_entry(reactions_tab, "Select the site of product 1 that needs to be already \nbound for the reaction to occur (optional):", [], 5, 0, width=10)
    product1_state = create_labeled_combobox_entry(reactions_tab, "Select product 1 state (optional):", [], 5, 2, width=10)
    product2_entry = create_labeled_combobox_entry_with_callback(reactions_tab, "Select product 2 if it is not bimolecular reaction or destruction reaction:", mol_names, 6, 0, on_product2_selected, width=10)
    product2_site = create_labeled_combobox_entry_with_callback(reactions_tab, "Select product 2 site:", [], 6, 2, on_product2_site_selected, width=10)
    product2_state = create_labeled_combobox_entry(reactions_tab, "Select product 2 state (optional):", [], 7, 0, width=10)
    product2_bound_site_requirment = create_labeled_combobox_entry(reactions_tab, "Select the site of product 2 that needs to be already \nbound for the reaction to occur (optional):", [], 7, 2, width=10)
    coord_com_reactant1_entry = create_labeled_entry_with_hint(reactions_tab, "Enter reactant 1's com coords \nfor bimolecular reaction ([x,y,z](nm)):", "[0,0,0]", 8, 0, width=10)
    coord_site_reactant1_entry = create_labeled_entry_with_hint(reactions_tab, "Enter reactant 1's reaction site coords \nfor bimolecular reaction ([x,y,z](nm)):", "", 8, 2, width=10)
    coord_com_reactant2_entry = create_labeled_entry_with_hint(reactions_tab, "Enter reactant 2's com coords \n(reactant 1 and 2 are in the same coordinate system) ([x,y,z](nm)):", "", 9, 0, width=10)
    coord_site_reactant2_entry = create_labeled_entry_with_hint(reactions_tab, "Enter reactant 2's reaction site coords: ([x,y,z](nm)):", "", 9, 2, width=10)
    onRate3Dka_entry = create_labeled_entry_with_hint(reactions_tab, "Enter microscopic on rate (nm^3 us^-1):", "", 10, 0, width=10)
    offRatekb_entry = create_labeled_entry_with_hint(reactions_tab, "Enter microscopic off rate (s^-1):", "", 10, 2, width=10)
    onRate3DMacro_entry = create_labeled_entry_with_hint(reactions_tab, "Enter macroscopic on rate if micro not provided (uM^-1 s^-1):", "", 11, 0, width=10)
    offRateMacro_entry = create_labeled_entry_with_hint(reactions_tab, "Enter macroscopic off rate if micro not provided (s^-1):", "", 11, 2, width=10)
    rate_entry = create_labeled_entry_with_hint(reactions_tab, "Enter rate for creation / destruction / uni statechange \n/ uni creation (M s^-1 / s^-1 / s^-1 / s^-1):", "", 12, 0, width=10)
    length3dto2d_entry = create_labeled_entry_with_hint(reactions_tab, "Enter length scale to convert 3D rate to 2D rate\n for bimolecular association (nm):", "10", 13, 0, width=10)
    bindRadSameCom_entry = create_labeled_entry_with_hint(reactions_tab, "Enter distance betweeen two reactants to force reaction \nwithinin the same complex for bimolecular association (sigma):", "1.1", 13, 2, width=10)
    loopCoopFactor_entry = create_labeled_entry_with_hint(reactions_tab, "Enter scale factor of rate when closing loops \nfor bimolecular association (ka):", "1", 14, 0, width=10)
    label_entry = create_labeled_entry_with_hint(reactions_tab, "Enter reaction label, used if you want to couple a \ndifferent reaction to this one:", "", 14, 2, width=10)
    coupledLabel_entry = create_labeled_entry_with_hint(reactions_tab, "Enter coupled reaction label, used if you allow \nthe completion of this reaction to cause another reaction:", "", 15, 0, width=10)
    kcat_entry = create_labeled_entry_with_hint(reactions_tab, "Enter rate of the coupled reaction (s^-1):", "", 15, 2, width=10)
    exclude_entry = create_labeled_entry_with_hint(reactions_tab, "Enter if check exclude volume for bound sites:", "False", 16, 0, width=10)
    reaction_type = ['Bimolecular Association: (A + B <-> A.B)', 'Zeroth Creation: null -> A', 'Destruction: A -> null', 'Michaelis-Menten: A + B <-> A.B -> A + C', 'Unimolecular Creation: A -> A + B', 'Bimolecular Statechange: A + B(b~U) -> A + B(b~P)', 'Statechange: B(b~U) -> B(b~P)']
    reaction_type_entry = create_labeled_combobox_entry(reactions_tab, "Select reaction type:", reaction_type, 17, 0)
    add_reaction_button = tk.Button(reactions_tab, text="Add Reaction", command=add_reaction)
    add_reaction_button.grid(row=18, column=1, pady=10)

    # all reactions added
    reactions = []

    tab_parent.pack(expand=True, fill=tk.BOTH)

    generate_button = tk.Button(root, text="Generate .inp File", command=generate_inp_file)
    generate_button.pack()



    root.mainloop()

def main():
    gui()

if __name__ == "__main__":
    main()
