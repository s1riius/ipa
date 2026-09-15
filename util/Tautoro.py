# ------------------------------------------------------------
# Tautomerisation engine
# ------------------------------------------------------------

"""Tautoro (Version 1.3)
Generate 3D tautomers from a 3D sdf-structure provided as input_sdf
Author: Dr. Tony Froehlich"""

import os
from rdkit import Chem
from rdkit.Chem import AllChem, QED, Crippen, rdmolops, Descriptors, rdMolDescriptors, Draw
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Geometry import Point3D
import numpy as np
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image,ImageDraw,ImageFont
import io
import math

#load input molecule from input_sdf
def load_template_mol(input_file):
    supplier = Chem.SDMolSupplier(input_file)
    for mol in supplier:
        if mol is not None:
            return mol
    raise ValueError("Could not load any molecule from the input SDF file.")

#generate all tautomers from input molecule
def generate_tautomers(input_file):
    supplier = Chem.SDMolSupplier(input_file)
    for mol in supplier:
        if mol is not None:
            enumerator = rdMolStandardize.TautomerEnumerator()
            tautomers = enumerator.Enumerate(mol)
            #for i, tautomer in enumerate(tautomers):
            #    print(f"Unsorted RD-Tautomer {i+1}: {Chem.MolToSmiles(tautomer)}")
            return tautomers
    raise ValueError("Could not load any molecule from the input SDF file.")

#write molecules to an SDF file, saving SMILES and name as properties
def write_sdf(mols,output_sdf,filename = None):
    with Chem.SDWriter(output_sdf) as writer:
        for mol in mols:
            if filename:
                mol.SetProp("_Name",filename)
            mol.SetProp("SMILES",Chem.MolToSmiles(mol))
            writer.SetForceV3000(False)                     #set to True to write SDF files in V3000 format instead of V2000
            writer.write(mol)

#add molecular properties such as MW, logP, TPSA, QED, etc.
def add_properties(mol, output_sdf):
    props = {'Formal_Charge': rdmolops.GetFormalCharge(mol),
             'Nr_Heavy_Atoms': mol.GetNumHeavyAtoms(),
             'MW': Descriptors.ExactMolWt(mol),
             'logP': Crippen.MolLogP(mol),
             'TPSA': rdMolDescriptors.CalcTPSA(mol),
             'QED': QED.qed(mol)}
    for prop_name,prop_value in props.items():
        mol.SetProp(prop_name,str(prop_value))
    writer = Chem.SDWriter(output_sdf)
    writer.SetForceV3000(False)                       #set to True to write SDF files in V3000 format instead of V2000
    writer.write(mol)
    writer.close()

#map atom indices between two molecules using substructure matching
def map_atoms(template_mol,tautomer_mol):
    match = tautomer_mol.GetSubstructMatch(template_mol)
    if not match:
        raise ValueError("No substructure match found between template and tautomer.")
    return {i: j for i,j in enumerate(match)}

#copy bond orders from original tautomer to new molecule
def copy_bond_orders(source_mol, target_mol, atom_map):
    for bond in source_mol.GetBonds():
        begin,end = bond.GetBeginAtomIdx(),bond.GetEndAtomIdx()
        mapped_begin,mapped_end = atom_map.get(begin),atom_map.get(end)
        if mapped_begin is not None and mapped_end is not None:
            tbond = target_mol.GetBondBetweenAtoms(mapped_begin,mapped_end)
            if tbond:
                tbond.SetBondType(bond.GetBondType())

#rebuild tautomer from SMILES and assign 3D coordinates from template
def process_tautomer(template_mol, tautomer_mol, raw_output_file, final_output_file):
    tautomer_smiles = Chem.MolToSmiles(tautomer_mol)
    tautomer_rebuilt = Chem.AddHs(Chem.MolFromSmiles(tautomer_smiles))
    AllChem.EmbedMolecule(tautomer_rebuilt,randomSeed = 0xf00d)
    tautomer_rebuilt.SetProp("Tautomer_SMILES",tautomer_smiles)
    atom_map = map_atoms(tautomer_mol,tautomer_rebuilt)
    copy_bond_orders(tautomer_mol,tautomer_rebuilt,atom_map)
    template_conf = template_mol.GetConformer()
    conf = Chem.Conformer(tautomer_rebuilt.GetNumAtoms())
    for tmpl_idx,taut_idx in atom_map.items():
        pos = template_conf.GetAtomPosition(tmpl_idx)
        conf.SetAtomPosition(taut_idx,pos)
    tautomer_rebuilt.RemoveAllConformers()
    tautomer_rebuilt.AddConformer(conf)
    write_sdf([tautomer_rebuilt],raw_output_file)
    finalize_tautomer_with_matching(template_mol,raw_output_file,final_output_file)
    return tautomer_rebuilt

#get coordinates for heavy atoms only
def get_heavy_atom_coords(mol):
    conf = mol.GetConformer()
    coords = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() > 1:
            idx = atom.GetIdx()
            pos = conf.GetAtomPosition(idx)
            coords.append((pos.x,pos.y,pos.z))
    return np.array(coords)

# For each heavy atom in template, find the matching heavy atom in target by 3D coordinates
def heavy_atom_coordinate_mapping(template_mol, target_mol, tolerance = 1e-3):
    template_coords = get_heavy_atom_coords(template_mol)
    target_coords = get_heavy_atom_coords(target_mol)
    permutation = []
    used = set()
    for t_coord in template_coords:
        found = False
        for idx, g_coord in enumerate(target_coords):
            if idx in used:
                continue
            if np.linalg.norm(np.array(t_coord) - np.array(g_coord)) < tolerance:
                permutation.append(idx)
                used.add(idx)
                found = True
                break
        if not found:
            raise ValueError("No coordinate match found for a heavy atom.")
    return permutation

#finalize molecule by restoring coordinates and adding hydrogens + properties
def finalize_tautomer_with_matching(template_mol, raw_sdf_file, final_output_sdf):
    supplier = Chem.SDMolSupplier(raw_sdf_file)
    raw_mol = next((mol for mol in supplier if mol is not None), None)
    if raw_mol is None:
        raise ValueError(f"Could not load molecule from {raw_sdf_file}")
    template_coords = [
        (template_mol.GetConformer().GetAtomPosition(i).x,
         template_mol.GetConformer().GetAtomPosition(i).y,
         template_mol.GetConformer().GetAtomPosition(i).z)
        for i in range(template_mol.GetNumAtoms())
    ]
    matched_indices = [i for i in range(template_mol.GetNumAtoms()) if template_mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
    editable_raw_mol = Chem.EditableMol(raw_mol)
    for j,(x,y,z) in enumerate(template_coords):
        if j < len(matched_indices):
            editable_raw_mol.GetMol().GetConformer().SetAtomPosition(matched_indices[j],Point3D(x,y,z))
    reordered_mol = editable_raw_mol.GetMol()
    mol_with_hydrogens = Chem.AddHs(reordered_mol,addCoords = True)
    #reorder heavy atoms based on coordinates
    try:
        perm = heavy_atom_coordinate_mapping(template_mol,mol_with_hydrogens)
        #RDKit RenumberAtoms only permutes all atoms, so let's build a full-atom permutation.
        #for heavy atoms: put them in template order. Hydrogens: leave in original order (append after heavy atoms).
        heavy_idxs = [atom.GetIdx() for atom in mol_with_hydrogens.GetAtoms() if atom.GetAtomicNum() > 1]
        hydro_idxs = [atom.GetIdx() for atom in mol_with_hydrogens.GetAtoms() if atom.GetAtomicNum() == 1]
        full_permutation = [heavy_idxs[i] for i in perm] + hydro_idxs
        mol_with_hydrogens = Chem.RenumberAtoms(mol_with_hydrogens,full_permutation)
    except Exception as e:
        print(f"Coordinate-based renumbering failed: {e}")
    #add_properties(mol_with_hydrogens, final_output_sdf)
    write_sdf([mol_with_hydrogens],final_output_sdf,filename = os.path.basename(final_output_sdf))
    return mol_with_hydrogens

#calculate MMFF94 energy for a molecule
def calculate_mmff_energy(mol):
#    mol = Chem.AddHs(mol, addCoords=True)                            #removed
#    if AllChem.EmbedMolecule(mol, randomSeed=0xf00d) != 0:           #removed
#        raise ValueError("MMFF embedding failed.")                   #removed
    props = AllChem.MMFFGetMoleculeProperties(mol)
    if props is None:
        raise ValueError("MMFF parameters missing.")
    ff = AllChem.MMFFGetMoleculeForceField(mol,props)
#    ff.Minimize()               #removed
    return ff.CalcEnergy()

#calculate energy using MMFF94, fallback to UFF if needed
def calculate_energy(mol,amide_filter = True,missing_amides = 0):
    try:
        energy = calculate_mmff_energy(mol)
    except Exception as e:
        print(f"MMFF failed, falling back to UFF: {e}")
#        mol = Chem.AddHs(mol, addCoords=True)                        #removed
#        AllChem.EmbedMolecule(mol, randomSeed=0xf00d)                #removed
        ff = AllChem.UFFGetMoleculeForceField(mol)
#        ff.Minimize()                                                #removed
        energy = ff.CalcEnergy()
    #apply penalty only if explicitly marked for penalization
    if amide_filter and missing_amides > 0:
        penalty = missing_amides*50.0
#       print(f"Penalty applied (amide-filter set to True). {missing_amides} amide group(s) detected in alternative tautomeric form (e.g. imide).") # removed
        energy += penalty
    return energy

#draw a single molecule with a legend and optional property block below
def draw_molecule_with_props(mol,legend,properties = None,mol_size = (600,600)):
    if mol is None:
        return Image.new("RGB",mol_size,"white")
    drawer = rdMolDraw2D.MolDraw2DCairo(mol_size[0],mol_size[1])
    opts = drawer.drawOptions()
    opts.bondLineWidth = 3
    opts.legendFontSize = 40
    rdMolDraw2D.PrepareAndDrawMolecule(drawer,mol,legend = legend)
    drawer.FinishDrawing()
    img_bytes = drawer.GetDrawingText()
    mol_img = Image.open(io.BytesIO(img_bytes)).convert("RGB")
    #draw property block below if any properties given
    if properties:
        line_height = 35                      #increase value if font size was increased
        text_block_height = len(properties)*line_height + 20
        total_height = mol_size[1] + text_block_height
        final_img = Image.new("RGB",(mol_size[0],total_height),"white")
        final_img.paste(mol_img,(0,0))
        draw = ImageDraw.Draw(final_img)
        try:
            font_regular = ImageFont.truetype("arial.ttf",30)         #change value to adjust font size
        except:
            font_regular = ImageFont.load_default()
        y_text = mol_size[1] + 10
        for key, val in properties:
            draw.text((10,y_text),f"{key}: {val}",fill = "black",font = font_regular)
            y_text += line_height
        return final_img
    else:
        return mol_img

#draws molecules in a grid with legends and (optional) properties below each molecule
def draw_molecule_grid(mols, legends, prop_lists = None, molsPerRow = 4, mol_size = (600,600), out_path = None):
    num_mols = len(mols)
    n_rows = math.ceil(num_mols/molsPerRow)
    # Draw all molecule images
    mol_imgs = []
    for i in range(num_mols):
        props = prop_lists[i] if prop_lists else None
        img = draw_molecule_with_props(mols[i],legends[i],props,mol_size)
        mol_imgs.append(img)
    #determine max height per row
    row_heights = []
    for r in range(n_rows):
        imgs_in_row = mol_imgs[r*molsPerRow:(r + 1)*molsPerRow]
        row_heights.append(max(img.height for img in imgs_in_row))
    #build grid
    grid_width = molsPerRow*mol_size[0]
    grid_height = sum(row_heights)
    grid_img = Image.new("RGB",(grid_width,grid_height),"white")
    y = 0
    for r in range(n_rows):
        imgs_in_row = mol_imgs[r*molsPerRow:(r + 1)*molsPerRow]
        x = 0
        row_height = row_heights[r]
        for img in imgs_in_row:
            grid_img.paste(img, (x, y))
            x += mol_size[0]
        y += row_height
    if out_path:
        grid_img.save(out_path)
        print(f"High-resolution tautomer grid image saved to {out_path}")
    return grid_img

#full tautomer processing workflow with sorting and visualization
def process_tautomers(input_file,output_dir_base,max_tautomers = 16,amide_filter = True,image_mode = "lowres"):
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_dir = os.path.join(output_dir_base,f"{base_name}_tautomers")
    os.makedirs(output_dir,exist_ok = True)
    template_mol = load_template_mol(input_file)
    tautomers = generate_tautomers(input_file)
    raw_sdf_files,finalized_sdf_files,energies = [],[],[]
    #detect number of amide groups across all tautomers
    amide_smarts = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX3](=O)")
    max_amide_count = max(len(mol.GetSubstructMatches(amide_smarts)) for mol in tautomers)
    for i,tautomer in enumerate(tautomers):
        raw_output_file = os.path.join(output_dir,f"{base_name}_raw_tautomer_{i + 1}.sdf")
        finalized_output_file = os.path.join(output_dir,f"{base_name}_finalized_tautomer_{i + 1}.sdf")
        try:
            process_tautomer(template_mol,tautomer,raw_output_file,finalized_output_file)
            mol = Chem.MolFromMolFile(finalized_output_file,removeHs = False)            #before mol = Chem.MolFromMolFile(finalized_output_file)
            if mol is None:
                print(f"Skipping tautomer {i + 1} due to failed molecule loading.")
                os.remove(raw_output_file)
                continue
            amide_count = len(mol.GetSubstructMatches(amide_smarts))
            missing_amides = max_amide_count - amide_count
            energy = calculate_energy(mol,amide_filter = amide_filter,missing_amides = missing_amides)
            energies.append((finalized_output_file,mol,energy))
            raw_sdf_files.append(raw_output_file)
            finalized_sdf_files.append(finalized_output_file)
        except Exception as e:
            print(f"Skipping tautomer {i+1} due to error: {e}")
            if os.path.exists(raw_output_file):
                os.remove(raw_output_file)
            continue
    #sort tautomers by energy (ascending)
    sorted_tautomers = sorted(energies,key = lambda x: x[2])[:max_tautomers]
    molecules,legends = [],[]
    template_smiles = Chem.MolToSmiles(Chem.RemoveHs(template_mol))
    #separate out the original tautomer BEFORE sorting
    original_idx = None
    for i,(file,mol,energy) in enumerate(energies):
        smiles = Chem.MolToSmiles(Chem.RemoveHs(mol))
        if smiles == template_smiles:
            original_idx = i
            break
    original_tautomer = None
    if original_idx is not None:
        original_tautomer = energies.pop(original_idx)  #remove from list for special handling
    #now sort the rest by energy (ascending) and keep up to (max_tautomers-1) since original will be added as rank 1
    sorted_tautomers = sorted(energies,key = lambda x: x[2])[:max_tautomers - 1]
    #insert the original at the top
    if original_tautomer is not None:
        sorted_tautomers = [original_tautomer] + sorted_tautomers
    molecules,legends = [],[]
    tautomer_files = []
    tautomer_smiles = []
    for rank, (file, mol, energy) in enumerate(sorted_tautomers,start = 1):
        sorted_filename = os.path.join(output_dir, f"{base_name}_tautomer_{rank}.sdf")
        tautomer_name = f"{base_name}_tautomer_{rank}"
        os.rename(file,sorted_filename)
        mol.SetProp("_Name",tautomer_name)
        #add_properties(mol, sorted_filename)
        #get SMILES string for the molecule
        #protonate
        mol_with_h = Chem.AddHs(mol, addCoords=True)
        smiles = Chem.MolToSmiles(mol_with_h)
        write_sdf([mol_with_h],sorted_filename,filename = tautomer_name)
        mol_2d = Chem.RemoveHs(mol)
        AllChem.Compute2DCoords(mol_2d)
        molecules.append(mol_2d)
        legends.append(f"{os.path.basename(sorted_filename)}\nEnergy: {energy:.2f}")
        tautomer_files.append(sorted_filename)
        tautomer_smiles.append(smiles)
        #print(f"Sorted tautomer {rank} saved as {sorted_filename}, SMILES: \"{Chem.MolToSmiles(mol_2d)}\", Energy: {energy:.2f}")
    #image generation options
    if image_mode == "highres":
        grid_img = draw_molecule_grid(molecules, legends,
                                      prop_lists = None,
                                      molsPerRow = 4,
                                      mol_size = (800,800),
                                      out_path = os.path.join(output_dir,f"{base_name}_tautomers_highres.png"))
        #grid_img.show()
    elif image_mode == "lowres":
        img = Draw.MolsToGridImage(molecules,legends = legends,molsPerRow = 4,subImgSize = (300,300))
        image_file = os.path.join(output_dir,f"{base_name}_tautomers.png")
        img.save(image_file)
        #img.show()
    #    print(f"Tautomers image saved to {image_file}")
    elif image_mode == "none":
        #explicitly do nothing (no image generation)
        pass
    #remove intermediate raw SDF files
    for raw_file in raw_sdf_files:
        if os.path.exists(raw_file):
            os.remove(raw_file)
    #remove intermediate finalized SDF files
    for finalized_file in finalized_sdf_files:
        if os.path.exists(finalized_file):
            os.remove(finalized_file)
    return tautomer_files,tautomer_smiles

# Entry point for script execution
if __name__ == "__main__":
    #path to input SDF file
    input_sdf = "./axitinib.sdf"
    #directory for output files
    output_dir = "./"
    os.makedirs(output_dir,exist_ok = True)
    #limit the maximum number of generated tautomers
    max_tautomers = 16
    #if set to True, amides will be ranked higher if they are in the proper amide form, 
    #rather than the imidic or any other alternative tautomeric form. 
    #set to False to disable amide preference.
    amide_filter = True
    #set to "lowres" or "highres" to choose image quality, or to "none" if no png-image should be exported.
    image_mode = "highres"
    process_tautomers(input_sdf,output_dir,max_tautomers,amide_filter = amide_filter,image_mode = image_mode)
